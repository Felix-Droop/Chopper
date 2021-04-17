#include <seqan3/std/ranges>     // range comparisons
#include <string>                // strings
#include <vector>                // vectors

#include <fstream>

#include <seqan3/test/tmp_filename.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, small_example)
{
    seqan3::test::tmp_filename binning_filename{"test.binning"};
    seqan3::test::tmp_filename counts_filename{"kmer_counts.txt"};
    
    {
        std::ofstream fout{counts_filename.get_path()};
        fout << data("small.fa").string() << '\t' << "585" << '\t' 
             << data("small.fa").string() << '\n';
    }

    {
        std::ofstream fout{binning_filename.get_path()};
        fout << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n" // header is ignored anyway
             << "SPLIT_BIN_0\t" << data("small.fa").string() << "\t2\t500\n"
             << "MERGED_BIN_2_0\t" << data("small.fa").string() << "\t2\t2500\n"
             << "MERGED_BIN_2_1\t" << data("small.fa").string() << "\t2\t2500\n"
             << "SPLIT_BIN_3\t" << data("small.fa").string() + "\t3\t1000\n";
    }

    cli_test_result result = execute_app("count_HIBF_kmers_based_on_binning",
                                         "-c", counts_filename.get_path().c_str(),
                                         "-k", "25",
                                         "-t", "1",
                                         "-f", binning_filename.get_path().c_str());

    std::string expected_stdout
    {
        "SPLIT_BIN_0\t292\t0\n"
        "SPLIT_BIN_0\t292\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "SPLIT_BIN_3\t195\t0\n"
        "MERGED_BIN_2\t585\t1170\n"
    };

    // compare results
    EXPECT_EQ(result.out, expected_stdout);
    EXPECT_EQ(result.err, std::string{});
}
