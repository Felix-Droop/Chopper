#pragma once

#include <chopper/pack/aggregate_by.hpp>
#include <chopper/pack/hierarchical_binning.hpp>
#include <chopper/pack/filenames_data_input.hpp>
#include <chopper/print_peak_memory_usage.hpp>

int set_up_and_parse_subparser_split(seqan3::argument_parser & parser, pack_config & config)
{
    parser.info.version = "1.0.0";
    parser.info.author = "Svenja Mehringer";
    parser.info.email = "svenja.mehringer@fu-berlin.de";

    parser.info.description.emplace_back("The `pack` submodule will create a hierarchical binning that minimizes the "
                                         "space consumption of the resulting Interleaved Bloom Filter that you may "
                                         "build with the `build` submodule using the results.");

    parser.add_option(config.data_file, 'f', "filenames",
                      "A tab separated file that contains the filepaths of sequence data you want to analyse.\n"
                      "The first column must contain the paths to sequence files separated by ';'.\n"
                      "The second column must contain the (kmer) count you want you data to be packed into bins. "
                      " See the submodule count for more details on how to add kmer counts to your sequences\n."
                      "All other columns are optional and can be used to aggregate your data (e.g. taxonmic ids).",
                      seqan3::option_spec::required);

    parser.add_option(config.bins, 'b', "technical-bins",
                      "Into how many technical bins do you want your sequence data to be packed?");
    
    parser.add_option(config.hll_dir, 'd', "hll-dir", 
                      "If given, the hll sketches are dumped to this directory and restored when they already exist.");

    parser.add_option(config.alpha, 'a', "alpha",
                      "The scaling factor to influence the number of merged bins.");

    parser.add_option(config.max_ratio, 'm', "max-ratio",
                      "The maximal cardinality ratio in the clustering intervals (must be < 1).");

    parser.add_option(config.num_threads, 't', "num-threads",
                      "The number of threads to use to compute merged HLL sketches.");

    parser.add_option(config.output_filename, 'o', "outfile",
                      "An output file name for the binning results.");

    parser.add_option(config.aggregate_by_column, 'y', "aggregate-by",
                      "Which column do you want to aggregate your files by? Start counting your columns from 1!",
                      seqan3::option_spec::standard,
                      seqan3::arithmetic_range_validator{3, std::numeric_limits<int>::max()});

    parser.add_flag(config.union_estimate, 'u', "union-estimate",
                    "Estimate the union of kmer sets to possibly improve the binning.");
    
    parser.add_flag(config.rearrange_bins, 'r', "rearrange-bins",
                    "Do a rearrangement of the bins which takes into account similarity. Only works with -u.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
        std::cerr << "[CHOPPER PACK ERROR] " << ext.what() << '\n'; // customize your error message
        return -2;
    }

    return 0;
}

int chopper_pack(seqan3::argument_parser & parser)
{
    pack_config config;

    if (auto r = set_up_and_parse_subparser_split(parser, config); r != 0)
        return r;

    // Read in the data file containing file paths, kmer counts and additional information.
    pack_data data;
    read_filename_data_file(data, config);

    // Some sanity checks on the input file and user options.
    if (data.filenames.empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] File Error: you passed an empty file."};

    if (config.aggregate_by_column != -1 && data.extra_information[0].empty())
        throw std::runtime_error{"[CHOPPER PACK ERROR] Aggregate Error: You want to aggregate by something but your "
                                 "file does not contain any extra information columns."};
    if (config.aggregate_by_column > static_cast<int>(data.extra_information[0].size()))
        throw std::runtime_error{"[CHOPPER PACK ERROR] Aggregate Error: You want to aggregate by a column index that is "
                                 "larger than the number of extra information columns."};

    // If requested, aggregate the data before packing them
    if (config.aggregate_by_column != -1)
        aggregate_by(data, config.aggregate_by_column);

    std::stringstream output_buffer;
    std::stringstream header_buffer;
    data.output_buffer = &output_buffer;
    data.header_buffer = &header_buffer;

    // Execute the actual algorithm:
    hierarchical_binning algo{data, config};
    algo.execute();

    // brief Write the output to the result file.
    std::ofstream fout{config.output_filename};
    fout << header_buffer.rdbuf();
    fout << output_buffer.rdbuf();

    print_peak_memory_usage();

    return 0;
}
