#pragma once

#include <vector>
#include <cmath>
#include <fstream>
#include <future>
#include <tuple>
#include <limits>

#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/range/views/async_input_buffer.hpp>

#include <chopper/union/hyperloglog.hpp>
#include <chopper/union/distance_matrix.hpp>
#include <chopper/union/clustering_node.hpp>

#include <robin_hood.h>

struct user_bin_sequence
{
private:
    //!\brief A reference to the filenames of the user input sequences.
    std::vector<std::string> & filenames;
    
    //!\brief A referece to kmer counts associated with the above files used to pack user bin into technical bins.
    std::vector<size_t> & user_bin_kmer_counts;

    //!\brief HyperLogLog sketches on the k-mer sets of the sequences from the files of filenames.
    std::vector<hyperloglog> sketches;

public:
    /*!\brief A sequence of user bins for which filenames and counts are given. This constructor reads 
              the HyperLogLog sketches from the hll_dir
     * \param[in] filenames_ filenames of the sequence files for the user bins
     * \param[in] user_bin_kmer_counts_ counts of the k-mer sets of the bins corresponding to filenames
     * \param[in] hll_dir path to the directory where hll caches will be found
     */
    user_bin_sequence(std::vector<std::string> & filenames_,
                      std::vector<size_t> & user_bin_kmer_counts_,
                      std::filesystem::path const & hll_dir) :
        filenames{filenames_},
        user_bin_kmer_counts{user_bin_kmer_counts_}
    {
        if (hll_dir.empty())
        {
            throw std::runtime_error("A directory where the HyperLogLog sketches are stored must be given "
                                     "when union estimates are enabled");
        }

        sketches.reserve(filenames.size());

        try 
        {
            for (auto & filename : filenames)
            {
                std::filesystem::path path = hll_dir / std::filesystem::path(filename).stem();
                path += ".hll";

                // the sketch bits will be automatically read from the files
                sketches.emplace_back();

                std::ifstream hll_fin(path, std::ios::binary);
                sketches.back().restore(hll_fin);
            }
        }
        catch (std::ios::failure const & fail)
        {
            std::cerr << "[CHOPPER PACK ERROR] Something went wrong trying to read the HyperLogLog sketches from files:\n"
                      << fail.what() << '\n';

        }
        catch (std::runtime_error const & err)
        {
            std::cerr << "[CHOPPER PACK ERROR] Something went wrong trying to read the HyperLogLog sketches from files:\n"
                      << err.what() << '\n';
        }
        catch (std::invalid_argument const & err)
        {
            std::cerr << "[CHOPPER PACK ERROR] Something went wrong trying to read the HyperLogLog sketches from files:\n"
                      << err.what() << '\n';
        }
    }

    /*!\brief For all intervals of filenames: estimate the cardinality of the union
     * of k-mer sets of all sequences in the files of the interval.
     * estimates[i][j] will be the union of the interval i, ..., i+j
     * \param[in] num_threads the number of threads to use
     * \param[out] estimates output table
     */
    void estimate_interval_unions(std::vector<std::vector<uint64_t>> & estimates, size_t num_threads)
    {
        estimates.clear();
        size_t const n = filenames.size();
        estimates.resize(n);

        auto indices = std::views::iota(0u, n)
                     | seqan3::views::async_input_buffer(num_threads * 5);
        
        auto worker = [&] ()
        {
            for (auto i : indices)
            {
                estimates[i].resize(n - i);
                estimates[i][0] = user_bin_kmer_counts[i];
                hyperloglog temp_hll = sketches[i];

                for (size_t j = i + 1; j < n; ++j)
                {
                    // merge next sketch into the current union
                    estimates[i][j - i] = static_cast<uint64_t>(temp_hll.merge_and_estimate_SIMD(sketches[j]));
                }
            }
        };

        // launch threads with worker
        std::vector<decltype(std::async(std::launch::async, worker))> handles;

        for (size_t i = 0; i < num_threads; ++i)
            handles.emplace_back(std::async(std::launch::async, worker));

        // wait for the threads to finish
        for (auto & handle : handles)
            handle.get();
    }

    /*!\brief Rearrange filenames, sketches and counts such that similar bins are close to each other
     * \param[in] max_ratio the maximal cardinality ratio in the clustering intervals (must be <= 1 and >= 0)
     */
    void rearrange_bins(double const max_ratio)
    {
        std::vector<size_t> permutation;

        size_t first = 0;
        size_t last = 1;

        while (first < filenames.size())
        {
            // size difference is too large or sequence is over -> do the clustering
            if (last == filenames.size() || user_bin_kmer_counts[first] * max_ratio > user_bin_kmer_counts[last])
            {
                // if this is not the first group, we want one bin overlap
                cluster_bins(permutation, first, last);
                first = last;
            }
            ++last;
        }

        // apply permutation to filenames and sketches
        for (size_t i = 0; i < permutation.size(); ++i)
        {
            size_t current = i;
            while (i != permutation[current])
            {
                size_t next = permutation[current];
                std::swap(filenames[current], filenames[next]);
                std::swap(user_bin_kmer_counts[current], user_bin_kmer_counts[next]);
                std::swap(sketches[current], sketches[next]);
                permutation[current] = current;
                current = next;
            }
            permutation[current] = current;
        }
    }

private:
    /*!\brief Perform an agglomerative clustering variant on the index range [first:last)
     * \param[in] first id of the first cluster of the interval
     * \param[in] last id of the last cluster of the interval plus one
     * \param[out] permutation append the new order to this
     */
    void cluster_bins(std::vector<size_t> & permutation,
                      size_t first,
                      size_t last)
    {
        robin_hood::unordered_map<size_t, clustering_node> clustering;
        robin_hood::unordered_map<size_t, double> estimates;
        distance_matrix dist(clustering, estimates);

        size_t const none = std::numeric_limits<size_t>::max();

        // initialize clustering and estimates
        for (size_t i = first; i < last; ++i)
        {
            clustering[i] = {none, none, sketches[i]};
            estimates[i] = sketches[i].estimate();
        }

        // if this is not the first group, we want to have one overlapping bin
        size_t previous_rightmost = none;
        if (first != 0)
        {
            previous_rightmost = permutation.back();
            clustering[previous_rightmost] = {none, none, sketches[previous_rightmost]};
            estimates[previous_rightmost] = sketches[previous_rightmost].estimate();
        }
        
        // initalize distance matrix
        dist.initialize(sketches);

        // ids of new nodes start at last's value
        size_t id = last;

        // keep merging nodes until we have a complete tree
        while (dist.size() > 1)
        {
            auto [min_id, neighbor_id] = dist.get_min_pair();

            // merge the two nodes with minimal distance together insert the new node into the clustering
            clustering[id] = {min_id, neighbor_id, std::move(clustering[min_id].hll)};
            estimates[id] = clustering[id].hll.merge_and_estimate_SIMD(clustering[neighbor_id].hll);

            // insert new cluster and update distances 
            dist.update(id, min_id, neighbor_id);

            ++id;
        }

        size_t final_root = dist.get_remaining_cluster_id();

        // rotate the previous rightmost to the left so that it has the correct place in the permutation
        if (first != 0)
        {
            rotate(clustering, previous_rightmost, final_root);
        }

        // traceback into permutation and ignore the previous rightmost
        trace(clustering, permutation, final_root, previous_rightmost);
    }

    /*!\brief Rotate the previous rightmost bin to the left of the clustering tree
     * \param[in, out] clustering the tree to do the rotation on
     * \param[in] previous_rightmost the id of the node to be rotated to the left
     * \param[in] id the id of the current node
     * 
     * \return whether previous rightmost was in the subtree rooted at id
     */
    bool rotate(robin_hood::unordered_map<size_t, clustering_node> & clustering,
                size_t const previous_rightmost,
                size_t const id)
    {
        if (id == previous_rightmost) return true;
        
        clustering_node & curr = clustering.at(id);
        if (curr.left == std::numeric_limits<size_t>::max())
        {
            return false;
        }

        // nothing to do if previous_rightmost is in the left subtree
        if(rotate(clustering, previous_rightmost, curr.left)) return true;
        
        // rotate if previous_rightmost is in the right subtree
        if(rotate(clustering, previous_rightmost, curr.right))
        {
            size_t const temp = curr.right;
            curr.right = curr.left;
            curr.left = temp;
            return true; 
        }

        // previous_rightmost is not in this subtree
        return false;
    }

    /*!\brief Do a recursive traceback to find the order of leaves in the clustering tree 
     * \param[in] clustering the tree to do the traceback on
     * \param[out] permutation append the new order to this
     * \param[in] id the id of the current node
     * \param[in] previous_rightmost the id of the node on the left which should be ignored
     */
    void trace(robin_hood::unordered_map<size_t, clustering_node> const & clustering,
               std::vector<size_t> & permutation,
               size_t const id,
               size_t const previous_rightmost)
    {
        clustering_node const & curr = clustering.at(id);
        
        if (curr.left == std::numeric_limits<size_t>::max())
        {
            if (id != previous_rightmost) permutation.push_back(id);
            return;
        }

        trace(clustering, permutation, curr.left, previous_rightmost);
        trace(clustering, permutation, curr.right, previous_rightmost);
    }
};