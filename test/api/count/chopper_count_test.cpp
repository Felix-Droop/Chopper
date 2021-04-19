#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <unordered_map>
#include <vector>

#include <seqan3/test/tmp_filename.hpp>

#include <chopper/count/count_config.hpp>
#include <chopper/count/chopper_count.hpp>

TEST(chopper_count_test, small_example_parallel_2_threads)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};
    seqan3::test::tmp_filename output_filename{"kmer_counts.txt"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
             << input_filename << '\t' << "TAX2\n";
    }

    const char * argv[] = {"./chopper-count", "-k", "15", "-w", "25", "-t", "2", "-c", "2",
                           "-f", data_filename.get_path().c_str(), 
                           "-o", output_filename.get_path().c_str()};
    int argc = 13;
    seqan3::argument_parser count_parser{"chopper-count", argc, argv, seqan3::update_notifications::off};

    std::string expected
    {
        input_filename + ";" + input_filename + "\t95\tTAX2\n" +
        input_filename + "\t95\tTAX1\n"
    };

    chopper_count(count_parser);

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}

TEST(chopper_count_test, disable_minimizers)
{
    std::string input_filename = DATADIR"small.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};
    seqan3::test::tmp_filename output_filename{"kmer_counts.txt"};

    // generate data filename
    {
        std::ofstream fout{data_filename.get_path()};
        fout << input_filename << '\t' << "TAX1\n"
             << input_filename << '\t' << "TAX2\n"
             << input_filename << '\t' << "TAX2\n";
    }

    const char * argv[] = {"./chopper-count",
                           "-k", "25",
                           "-t", "1",
                           "-c", "2",
                           "--disable-minimizers",
                           "-f", data_filename.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 12;
    seqan3::argument_parser count_parser{"chopper-count", argc, argv, seqan3::update_notifications::off};

    std::string expected
    {
        input_filename + ";" + input_filename + "\t585\tTAX2\n" +
        input_filename + "\t585\tTAX1\n"
    };


    chopper_count(count_parser);

    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());

    EXPECT_EQ(expected, output_file_str);
}
