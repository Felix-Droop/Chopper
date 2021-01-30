#pragma once

#include <seqan3/std/filesystem>

struct pack_config
{
    std::filesystem::path data_file;
    std::filesystem::path output_filename{"binning.out"};
    uint16_t bins{64};
    uint8_t k{25};
    uint8_t sketch_bits{12};
    size_t num_threads{std::thread::hardware_concurrency()};
    double alpha{10};
    int aggregate_by_column{-1};
    bool union_estimate{false};
    bool resort_bins{false};
};
