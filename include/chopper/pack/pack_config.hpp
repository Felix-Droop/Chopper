#pragma once

#include <seqan3/std/filesystem>

#include <thread>

struct pack_config
{
    std::filesystem::path data_file;
    std::filesystem::path output_filename{"binning.out"};
    uint16_t bins{64};
    std::filesystem::path hll_dir{};
    double alpha{10};
    double max_ratio{0.5};
    int aggregate_by_column{-1};
    bool union_estimate{false};
    bool resort_bins{false};
};
