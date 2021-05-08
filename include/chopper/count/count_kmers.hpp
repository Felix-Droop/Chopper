#include <future>
#include <fstream>
#include <thread>

#define SEQAN_HAS_ZLIB 1
#define SEQAN3_HAS_ZLIB 1

#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/utility/views/to.hpp>
#include <seqan3/std/filesystem>

#include <chopper/union/hyperloglog.hpp>

#include <robin_hood.h>

std::mutex mu;

void print_safely(std::pair<std::string, std::vector<std::string>> const & cluster,
                  uint64_t size,
                  std::ofstream & fout)
{
    std::lock_guard<std::mutex> lock(mu);  // Acquire the mutex
    assert(cluster.second.size() >= 1);

    fout << cluster.second[0]; // write first filename
    for (size_t i = 1; i < cluster.second.size(); ++i)
        fout << ";" << cluster.second[i];
    fout << '\t' << size << '\t' << cluster.first << std::endl;;

}// lock_guard object is destroyed and mutex mu is released

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

template <typename cluster_view_type, typename compute_view_type>
void compute_hashes(cluster_view_type && cluster_view, compute_view_type && compute_fn, std::ofstream & fout,
                    bool const exclusively_hlls, std::filesystem::path const & hll_dir, uint8_t sketch_bits)
{
    for (auto && [cluster, sequence_vector] : cluster_view)
    {
        robin_hood::unordered_set<uint64_t> result;
        hyperloglog sketch(sketch_bits);

        for (auto && seq : sequence_vector)
        {
            for (auto hash : seq | compute_fn)
            {
                if (!exclusively_hlls)
                {
                    result.insert(hash);
                }
                if (exclusively_hlls || !hll_dir.empty())
                {
                    sketch.add(reinterpret_cast<char*>(&hash), sizeof(hash));
                }
            }
        }

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t size = exclusively_hlls ? static_cast<uint64_t>(sketch.estimate()) : result.size();
        print_safely(cluster, size, fout);

        if (!hll_dir.empty())
        {
            // For more than one file in the cluster, Felix doesn't know how to name the file
            // and what exactly is supposed to happen.
            if (cluster.second.size() > 1)
            {
                throw std::runtime_error("This mode sadly was not implemented yet for multiple files grouped together.");
            }
            
            // For one file in the cluster, the file stem is used with the .hll ending
            std::filesystem::path path = hll_dir / std::filesystem::path(cluster.first).stem();
            path += ".hll";
            std::ofstream hll_fout(path, std::ios::binary);
            sketch.dump(hll_fout);
        }
    }
}

void count_kmers(std::unordered_map<std::string, std::vector<std::string>> const & filename_clusters,
                 count_config const & config)
{
    // output file
    std::ofstream fout{config.output_filename};

    // create the hll dir if it doesn't already exist
    if (!config.hll_dir.empty() && !std::filesystem::exists(config.hll_dir))
    {
        std::filesystem::create_directory(config.hll_dir);
    }

    size_t const counting_threads = (config.num_threads <= 1) ? 1 : config.num_threads - 1;

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k}, seqan3::window_size{config.w});
    auto compute_kmers = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    auto read_files = std::views::transform([] (auto const & cluster)
    {
        using result_t = std::vector<std::vector<seqan3::dna4>>;
        using pair_t = std::pair<std::pair<std::string, std::vector<std::string>>, result_t>;

        result_t result;
        for (auto const & filename : cluster.second)
        {
            sequence_file_type fin{filename};
            for (auto & [seq] : fin)
                result.push_back(seq);
        }

        return pair_t{cluster, result};
    });

    auto cluster_view = filename_clusters
                      | read_files
                      | seqan3::views::async_input_buffer(counting_threads);

    auto worker = [&config, &cluster_view, &compute_minimiser, &compute_kmers, &fout] ()
    {
        if (config.disable_minimizers)
            compute_hashes(cluster_view, compute_kmers, fout, config.exclusively_hlls, config.hll_dir, config.sketch_bits);
        else
            compute_hashes(cluster_view, compute_minimiser, fout, config.exclusively_hlls, config.hll_dir, config.sketch_bits);
    };

    // launch threads with worker
    std::vector<decltype(std::async(std::launch::async, worker))> handles;

    for (size_t i = 0; i < counting_threads; ++i)
        handles.emplace_back(std::async(std::launch::async, worker));

    // wait for the threads to finish to measure peak memory usage afterwards
    for (auto & handle : handles)
        handle.get();
}
