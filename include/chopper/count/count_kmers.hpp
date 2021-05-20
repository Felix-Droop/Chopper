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

void write_cluster_data(std::pair<std::string, std::vector<std::string>> const & cluster,
                        uint64_t size,
                        std::ofstream & fout)
{
    assert(cluster.second.size() >= 1);

    fout << cluster.second[0]; // write first filename
    for (size_t i = 1; i < cluster.second.size(); ++i)
        fout << ";" << cluster.second[i];
    fout << '\t' << size << '\t' << cluster.first << std::endl;
}

struct mytraits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

template <typename seq_type, typename compute_view_type>
void compute_hashes(seq_type && seq, compute_view_type && compute_fn, count_config const & config,
                    robin_hood::unordered_node_set<uint64_t> & result, hyperloglog & sketch)
{
    for (auto hash : seq | compute_fn)
    {
        if (!config.exclusively_hlls)
        {
            result.insert(hash);
        }
        if (config.exclusively_hlls || !config.hll_dir.empty())
        {
            sketch.add(reinterpret_cast<char*>(&hash), sizeof(hash));
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

    auto compute_minimiser = seqan3::views::minimiser_hash(seqan3::ungapped{config.k}, seqan3::window_size{config.w});
    auto compute_kmers = seqan3::views::kmer_hash(seqan3::ungapped{config.k});

    // copy filename clusters to vector
    std::vector<std::pair<std::string, std::vector<std::string>>> cluster_vector;
    for (auto const & cluster : filename_clusters)
    {
        cluster_vector.push_back(cluster);
    }

    #pragma omp parallel for schedule(static) num_threads(config.num_threads)
    for (size_t i = 0; i < cluster_vector.size(); ++i)
    {
        // read files
        std::vector<std::vector<seqan3::dna4>> sequence_vector;
        for (auto const & filename : cluster_vector[i].second)
        {
            sequence_file_type fin{filename};
            for (auto & [seq] : fin)
                sequence_vector.push_back(seq);
        }

        robin_hood::unordered_node_set<uint64_t> result;
        hyperloglog sketch(config.sketch_bits);

        for (auto && seq : sequence_vector)
        {
            if (config.disable_minimizers)
                compute_hashes(seq, compute_kmers, config, result, sketch);
            else
                compute_hashes(seq, compute_minimiser, config, result, sketch);
        }

        // print either the exact or the approximate count, depending on exclusively_hlls
        uint64_t size = config.exclusively_hlls ? static_cast<uint64_t>(sketch.estimate()) : result.size();

        #pragma omp critical
        write_cluster_data(cluster_vector[i], size, fout);

        if (!config.hll_dir.empty())
        {
            // For more than one file in the cluster, Felix doesn't know how to name the file
            // and what exactly is supposed to happen.
            if (cluster_vector[i].second.size() > 1)
            {
                throw std::runtime_error("This mode sadly was not implemented yet for multiple files grouped together.");
            }
            
            // For one file in the cluster, the file stem is used with the .hll ending
            std::filesystem::path path = config.hll_dir / std::filesystem::path(cluster.first).stem();
            path += ".hll";
            std::ofstream hll_fout(path, std::ios::binary);
            sketch.dump(hll_fout);
        }
    }
}
