#include <thread>
#include <future>
#include <vector>
#include <utility>
#include <tuple>
#include <iostream>
#include <chrono>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/views/async_input_buffer.hpp>

#include <chopper/count/count_config.hpp>
#include <chopper/count/read_data_file.hpp>
#include <chopper/print_peak_memory_usage.hpp>

int main(int argc, const char *argv [])
{
    auto start = std::chrono::high_resolution_clock::now();

    seqan3::argument_parser parser{"count_GC_content", argc, argv,
                                   seqan3::update_notifications::off};
    
    count_config config{};

    parser.add_option(config.data_file, 'f', "data_file", 
                      "Give me a filename to a seqinfo file.", seqan3::option_spec::required);
    parser.add_option(config.num_threads, 't', "threads", "Number of threads.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cerr << "[COMMAND LINE INPUT ERROR] " << ext.what() << std::endl;
        return -1;
    }

    size_t const counting_threads = (config.num_threads <= 1) ? 1 : config.num_threads - 1;
    auto filename_clusters = read_data_file(config);

    struct mytraits : public seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = seqan3::dna4;
    };

    using sequence_file_type = seqan3::sequence_file_input<mytraits,
                                                       seqan3::fields<seqan3::field::seq>,
                                                       seqan3::type_list<seqan3::format_fasta, seqan3::format_fastq>>;

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

    auto worker = [&cluster_view] ()
    { 
        uint64_t sum_GC = 0;
        uint64_t total = 0;

        seqan3::dna4 cytosin;
        cytosin.assign_char('C');
        seqan3::dna4 guanin;
        cytosin.assign_char('G');

        for (auto && [cluster, sequence_vector] : cluster_view)
        {
            for (auto && seq : sequence_vector)
            {
                for (seqan3::dna4 base : seq)
                {
                    sum_GC += static_cast<uint64_t>(base == cytosin || base == guanin);
                    ++total;
                }
            }
        }

        return std::make_tuple(sum_GC, total);
    };

    // launch threads with worker
    std::vector<decltype(std::async(std::launch::async, worker))> handles;

    for (size_t i = 0; i < counting_threads; ++i)
        handles.emplace_back(std::async(std::launch::async, worker));

    double sum_GC = 0;
    double total = 0;

    // wait for the threads to finish to measure peak memory usage afterwards
    for (auto & handle : handles)
    {
        auto [sum_GC_part, total_part] = handle.get();
        sum_GC += sum_GC_part;
        total += total_part;
    }

    std::cout << "GC content: " << sum_GC / total << std::endl;

    auto dur = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start);
    std::cerr << "Took " << dur.count() << " seconds.\n";

    print_peak_memory_usage();
}