#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/std/filesystem>
#include <seqan3/test/tmp_filename.hpp>

#include <chopper/split/split_data.hpp>
#include <chopper/split/sequence_input.hpp>
#include <chopper/split/chopper_split.hpp>

TEST(chopper_split_test, simple_example)
{
    std::string input_filename1 = DATADIR"small.fa";
    std::string input_filename2 = DATADIR"small2.fa";
    seqan3::test::tmp_filename output_filename{"small.split"};
    const char * argv[] = {"./chopper-split", "-k", "15", "-w", "25", "-b", "3",
                           "-s", input_filename1.c_str(), "-s", input_filename2.c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 13;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, seqan3::update_notifications::off};

    chopper_split(split_parser);

    // compare results
    std::string expected_file_str
    {
        "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n" +
        input_filename1 + "\tseq1\t0\t163\t0\t-\n" +
        input_filename1 + "\tseq2\t0\t186\t0\t-\n" +
        input_filename1 + "\tseq3\t0\t163\t0\t-\n" +
        input_filename2 + "\tseq10\t0\t163\t0\t-\n" +
        input_filename2 + "\tseq20\t0\t186\t0\t-\n" +
        input_filename2 + "\tseq30\t0\t163\t0\t-\n" +
        input_filename1 + "\tseq1\t163\t247\t1\t-\n" +
        input_filename1 + "\tseq2\t186\t327\t1\t-\n" +
        input_filename1 + "\tseq3\t163\t284\t1\t-\n" +
        input_filename2 + "\tseq10\t163\t247\t1\t-\n" +
        input_filename2 + "\tseq20\t186\t327\t1\t-\n" +
        input_filename2 + "\tseq30\t163\t284\t1\t-\n" +
        input_filename1 + "\tseq1\t247\t400\t2\t-\n" +
        input_filename1 + "\tseq2\t327\t480\t2\t-\n" +
        input_filename1 + "\tseq3\t284\t481\t2\t-\n" +
        input_filename2 + "\tseq10\t247\t400\t2\t-\n" +
        input_filename2 + "\tseq20\t327\t480\t2\t-\n" +
        input_filename2 + "\tseq30\t284\t481\t2\t-\n"
    };

    // compare results
    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file_str);
}

TEST(chopper_split_test, no_s_or_f_option)
{
    std::string input_filename = DATADIR"small.fa";
    const char * argv[] = {"./chopper-split", "-k", "15"};
    int argc = 3;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, seqan3::update_notifications::off};

    EXPECT_THROW(chopper_split(split_parser), std::runtime_error);
}

TEST(chopper_split_test, data_file_as_input)
{
    std::string input_filename1 = DATADIR"small.fa";
    std::string input_filename2 = DATADIR"small2.fa";
    seqan3::test::tmp_filename data_filename{"data.tsv"};

    {
        std::ofstream fout{data_filename.get_path()};
        fout << "#BIN_ID\tSEQ_IDS\tNUM_TECHNICAL_BINS\tESTIMATED_MAX_TB_SIZE\n"
             << "SPLIT_BIN_0\t" << input_filename1 + "\t2\t500\n"
             << "SPLIT_BIN_2\t" << input_filename1 + "\t2\t500\n"
             << "MERGED_BIN_4_0\t" << input_filename1 << "\t2\t2500\n"
             << "MERGED_BIN_4_2\t" << input_filename1 << ";" << input_filename2 << "\t2\t2500\n"
             << "MERGED_BIN_4_4\t" << input_filename1 << "\t1\t500\n"
             << "SPLIT_BIN_5\t" << input_filename2 + "\t3\t1000\n";
    }

    seqan3::test::tmp_filename output_filename{"small.split"};

    const char * argv[] = {"./chopper-split", "-k", "15", "-w", "25",
                           "-f", data_filename.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 9;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, seqan3::update_notifications::off};

    EXPECT_EQ(chopper_split(split_parser), 0);

    std::string const expected_output_str
    {
        "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n" +
        /*SPLIT_BIN_0*/
        input_filename1 + "\tseq1\t0\t209\t0\t-\n" +
        input_filename1 + "\tseq2\t0\t289\t0\t-\n" +
        input_filename1 + "\tseq3\t0\t209\t0\t-\n" +
        input_filename1 + "\tseq1\t209\t400\t1\t-\n" +
        input_filename1 + "\tseq2\t289\t480\t1\t-\n" +
        input_filename1 + "\tseq3\t209\t481\t1\t-\n" +
        /*SPLIT_BIN_2*/
        input_filename1 + "\tseq1\t0\t209\t2\t-\n" +
        input_filename1 + "\tseq2\t0\t289\t2\t-\n" +
        input_filename1 + "\tseq3\t0\t209\t2\t-\n" +
        input_filename1 + "\tseq1\t209\t400\t3\t-\n" +
        input_filename1 + "\tseq2\t289\t480\t3\t-\n" +
        input_filename1 + "\tseq3\t209\t481\t3\t-\n" +
        /*MERGED_BIN_4_0*/
        input_filename1 + "\tseq1\t0\t209\t4\t0\n" +
        input_filename1 + "\tseq2\t0\t289\t4\t0\n" +
        input_filename1 + "\tseq3\t0\t209\t4\t0\n" +
        input_filename1 + "\tseq1\t209\t400\t4\t1\n" +
        input_filename1 + "\tseq2\t289\t480\t4\t1\n" +
        input_filename1 + "\tseq3\t209\t481\t4\t1\n" +
        /*MERGED_BIN_4_1*/
        input_filename1 + "\tseq1\t0\t209\t4\t2\n" +
        input_filename1 + "\tseq2\t0\t289\t4\t2\n" +
        input_filename1 + "\tseq3\t0\t209\t4\t2\n" +
        input_filename2 + "\tseq10\t0\t209\t4\t2\n" +
        input_filename2 + "\tseq20\t0\t289\t4\t2\n" +
        input_filename2 + "\tseq30\t0\t209\t4\t2\n" +
        input_filename1 + "\tseq1\t209\t400\t4\t3\n" +
        input_filename1 + "\tseq2\t289\t480\t4\t3\n" +
        input_filename1 + "\tseq3\t209\t481\t4\t3\n" +
        input_filename2 + "\tseq10\t209\t400\t4\t3\n" +
        input_filename2 + "\tseq20\t289\t480\t4\t3\n" +
        input_filename2 + "\tseq30\t209\t481\t4\t3\n" +
        /*MERGED_BIN_4_2*/
        input_filename1 + "\tseq1\t0\t400\t4\t4\n" +
        input_filename1 + "\tseq2\t0\t480\t4\t4\n" +
        input_filename1 + "\tseq3\t0\t481\t4\t4\n" +
        /*SPLIT_BIN_5*/
        input_filename2 + "\tseq10\t0\t163\t5\t-\n" +
        input_filename2 + "\tseq20\t0\t186\t5\t-\n" +
        input_filename2 + "\tseq30\t0\t163\t5\t-\n" +
        input_filename2 + "\tseq10\t163\t247\t6\t-\n" +
        input_filename2 + "\tseq20\t186\t327\t6\t-\n" +
        input_filename2 + "\tseq30\t163\t284\t6\t-\n" +
        input_filename2 + "\tseq10\t247\t400\t7\t-\n" +
        input_filename2 + "\tseq20\t327\t480\t7\t-\n" +
        input_filename2 + "\tseq30\t284\t481\t7\t-\n"
    };

    // compare results
    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_output_str);
}


TEST(chopper_split_test, big_fat_nodes)
{
    // big fat nodes (e.g. the extreme is a single node if only one sequence was given)
    // should be split several times, not only once.

    seqan3::test::tmp_filename seq_file{"one_seq.fasta"};

    {
        std::ofstream fout{seq_file.get_path()};
        fout << ">seq1\n"
             << "ACTGATCAGGGAGCTAGCAGGCAGGCAGCAGCTAGCGAGCGATCGAGCATCGAGCATCGAGCGATCGACGATCGACTAGC\n";
    }

    seqan3::test::tmp_filename output_filename{"small.split"};
    const char * argv[] = {"./chopper-split", "-k", "5", "-w", "7", "-b", "5",
                           "-s", seq_file.get_path().c_str(),
                           "-o", output_filename.get_path().c_str()};
    int argc = 11;
    seqan3::argument_parser split_parser{"chopper-split", argc, argv, seqan3::update_notifications::off};

    chopper_split(split_parser);

    // compare results
    std::string expected_file_str
    {
        "#FILE_ID\tSEQ_ID\tBEGIN\tEND\tHIBF_BIN_IDX\tLIBF_BIN_IDX\n" +
        seq_file.get_path().string() + "\tseq1\t0\t17\t0\t-\n" +
        seq_file.get_path().string() + "\tseq1\t17\t34\t1\t-\n" +
        seq_file.get_path().string() + "\tseq1\t34\t51\t2\t-\n" +
        seq_file.get_path().string() + "\tseq1\t51\t68\t3\t-\n" +
        seq_file.get_path().string() + "\tseq1\t68\t80\t4\t-\n"
    };

    // compare results
    std::ifstream output_file{output_filename.get_path()};
    std::string const output_file_str((std::istreambuf_iterator<char>(output_file)), std::istreambuf_iterator<char>());
    EXPECT_EQ(output_file_str, expected_file_str);
}
