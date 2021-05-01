#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <unordered_set>
#include <vector>
#include <fstream>
#include <iostream>

#include <chopper/union/hyperloglog.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/test/tmp_filename.hpp>

struct input_traits : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;
    using sequence_legal_alphabet = char;
    using sequence_contianer = std::vector<sequence_alphabet>;
};
using sequence_file_type = seqan3::sequence_file_input<input_traits, seqan3::fields<seqan3::field::seq>>;

TEST(hyperloglog_test, initialization)
{
    size_t const b = 4;
    size_t const m = 1 << b;

    hyperloglog sketch(b);
    
    EXPECT_EQ(sketch.registerSize(), m);

    // No elements were inserted, so the small values correction should be used.
    // Since there are only zeros in the register, the correction formula should be:
    // m * log(m / #zeros) = m * log(m/m) = m * log(1) = 0
    EXPECT_EQ(sketch.estimate(), 0.0);
}

TEST(hyperloglog_test, add_and_estimate_small)
{
    size_t const b = 4;
    size_t const m = 1 << b;

    hyperloglog sketch(b);

    // XXH3_64bits hash -> first 4 bits: 0000, rank: 3
    sketch.add("bla", 3); 
    // XXH3_64bits hash -> first 4 bits: 1011, rank: 2
    sketch.add("bli", 3);
    // XXH3_64bits hash -> first 4 bits: 0100, rank: 2
    sketch.add("blub", 4);
    // XXH3_64bits hash -> first 4 bits: 1110, rank: 1
    sketch.add("bloink", 6);
    // XXH3_64bits hash -> first 4 bits: 0000, rank: 3
    sketch.add("blubba", 6); 
    // XXH3_64bits hash -> first 4 bits: 1101, rank: 1
    sketch.add("blumpf", 6);
    // XXH3_64bits hash -> first 4 bits: 1001, rank: 2
    sketch.add("blarkse", 7);
    // XXH3_64bits hash -> first 4 bits: 1000, rank: 2
    sketch.add("bladuzel", 8);

    // estimate = alpha * m  * m  / sum(2^(-M_[j]))
    //          = 0.673 * 16 * 16 / (89/8) = 15,48...

    // this still is in the range of small value corrections (< 2.5 * 16)
    // m * log(#zeros / m) = 16 * log(16/9) = 9.205826318 (with calculator)

    EXPECT_NEAR(sketch.estimate(), 9.205826318, 0.0000001);
}

TEST(hyperloglog_test, add_and_estimate_large)
{
    std::string input_file = DATADIR"small.fa";
    size_t const k = 16;

    size_t const b = 4;
    size_t const m = 1 << b;
    hyperloglog sketch(b);

    std::unordered_set<std::string> control;

    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the sketch
    for (auto && [seq] : seq_file)
    {
        // we have to go C-style here for the HyperLogLog Interface
        char* c_seq_it = &*seq.begin();
        char* end = c_seq_it + seq.size();

        while(c_seq_it + k <= end)
        {
            control.insert(std::string(c_seq_it, k));
            sketch.add(c_seq_it, k);
            ++c_seq_it;
        }
    }

    // the estimate is greater than 2.5 * 16, therefore this is the raw estimate
    // 1.04 / sqrt(m) is the usual relative error bound proved in the paper 
    EXPECT_NEAR(sketch.estimate(), control.size(), control.size() * 1.04 / 4);
}

TEST(hyperloglog_test, merge_and_merge_SSE)
{
    std::string input_file = DATADIR"small.fa";
    size_t const k = 16;

    size_t const b = 4;
    size_t const m = 1 << b;
    hyperloglog full_sketch(b);
    hyperloglog merge_sketch(b);
    hyperloglog merge_SSE_sketch(b);

    std::vector<hyperloglog> partial_sketches;

    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the full_sketch
    // and add a disjointed sketch for every sequence to partial_sketches
    for (auto && [seq] : seq_file)
    {
        partial_sketches.emplace_back(b);

        // we have to go C-style here for the HyperLogLog Interface
        char* c_seq_it = &*seq.begin();
        char* end = c_seq_it + seq.size();

        while(c_seq_it + k <= end)
        {
            partial_sketches.back().add(c_seq_it, k);
            full_sketch.add(c_seq_it, k);
            ++c_seq_it;
        }
    }   

    double merge_SSE_estimate;
    // merge all partial sketches into merge_sketch
    for (auto & partial_sketch : partial_sketches)
    {
        merge_sketch.merge(partial_sketch);
        merge_SSE_estimate = merge_SSE_sketch.merge_and_estimate_SSE(partial_sketch);
    }

    // now full_sketch and merged_sketch should be equal
    EXPECT_EQ(full_sketch.estimate(), merge_sketch.estimate());
    EXPECT_EQ(full_sketch.estimate(), merge_SSE_sketch.estimate());
}

TEST(hyperloglog_test, dump_and_restore)
{
    std::string input_file = DATADIR"small.fa";
    size_t const k = 16;

    size_t const b = 4;
    size_t const m = 1 << b;
    hyperloglog dump_sketch(b);
    hyperloglog restore_sketch(b);

    sequence_file_type seq_file{input_file};

    // put every sequence in this file into the dump_sketch
    for (auto && [seq] : seq_file)
    {
        // we have to go C-style here for the HyperLogLog Interface
        char* c_seq_it = &*seq.begin();
        char* end = c_seq_it + seq.size();

        while(c_seq_it + k <= end)
        {
            dump_sketch.add(c_seq_it, k);
            ++c_seq_it;
        }
    }

    // create temp file
    seqan3::test::tmp_filename dump_filename{"dump.hll"};
    
    // dump sketch
    std::ofstream ostrm(dump_filename.get_path(), std::ios::binary);
    dump_sketch.dump(ostrm);

    // restore sketch
    std::ifstream istrm(dump_filename.get_path(), std::ios::binary);
    restore_sketch.restore(istrm);

    // now dump_sketch and restore_sketch should be equal
    EXPECT_EQ(dump_sketch.estimate(), restore_sketch.estimate());
}