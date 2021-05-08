#pragma once

#include <vector>
#include <queue>
#include <functional>
#include <limits>
#include <fstream>

#include <future>
#include <barrier>

#include <cmath>
#include <cstddef>
#include <cassert>

#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/io/views/async_input_buffer.hpp>

#include <chopper/union/hyperloglog.hpp>

#include <robin_hood.h>

struct user_bin_sequence
{
private:
    //!\brief type for a node in the clustering tree when for the rearrangement
    struct clustering_node 
    {
        // children in the tree
        size_t left;
        size_t right;
        // hll sketch of the union if the node is still a root
        hyperloglog hll;
    };

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
    void estimate_interval_unions(std::vector<std::vector<uint64_t>> & estimates, size_t const num_threads)
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
     * \param[in] num_threads the number of threads to use
     */
    void rearrange_bins(double const max_ratio, size_t const num_threads)
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
                cluster_bins(permutation, first, last, num_threads);
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
     * \param[in] num_threads the number of threads to use
     * \param[out] permutation append the new order to this
     */
    void cluster_bins(std::vector<size_t> & permutation,
                      size_t const first,
                      size_t const last,
                      size_t const num_threads)
    {
        assert(num_threads >= 1);

        //!\brief element of the second priority queue layer of the distance matrix
        struct neighbor
        {
            double dist;
            size_t id;

            bool operator>(neighbor const & other) const
            {
                return dist > other.dist;
            }
        };

        //!\brief type of a min heap based priority queue
        using prio_queue = std::priority_queue<neighbor, std::vector<neighbor>, std::greater<neighbor>>;

        /*!\brief internal map that stores the distances
        *
        * The first layer is a hash map with the ids of active clusters as keys.
        * The values (second layer) are priority queues with neighbors of the cluster 
        * with the respective key in the first layer.
        * These neighbors are themselves clusters with an id and store a distance to the 
        * cluster of the first layer.
        */
        robin_hood::unordered_map<size_t, prio_queue> dist;
        
        // clustering tree
        robin_hood::unordered_map<size_t, clustering_node> clustering;

        // cache for hll cardinality estimates
        robin_hood::unordered_map<size_t, double> estimates;

        size_t const none = std::numeric_limits<size_t>::max();

        // every thread will write its observed id with minimal distance to some other here
        // id == none means that the thread observed only empty or no priority queues
        std::vector<size_t> min_ids(num_threads, none);
        
        // these will be the new ids for new clusters
        // the first one is invalid, but it will be incremented before it is used for the first time
        size_t id = last - 1;

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
        
        // only the keys are needed when iterating through the clustering later
        auto to_only_key = [] (robin_hood::unordered_map<size_t, clustering_node>::value_type const & v)
        {
            return v.first;
        };

        // async buffer to iterate through the clustering in parallel during the initialization
        auto clustering_buffer = clustering 
                               | std::views::transform(to_only_key)
                               | seqan3::views::async_input_buffer(num_threads * 5);
        
        // a key and pointer is needed when iterating through dist later, because the 
        // async_input_buffer moves the values from the range
        auto to_key_and_pointer = [] (robin_hood::unordered_map<size_t, prio_queue>::value_type & v)
        {
            return std::make_tuple(v.first, &v.second);
        };

        // async buffer to iterate through the distance matrix in parallel when updating it
        // this assignment here is just to declare the variable with the correct type
        // it will be reset in every iteration in the end of the critical clustering update
        auto dist_buffer = dist 
                         | std::views::transform(to_key_and_pointer)
                         | seqan3::views::async_input_buffer(num_threads * 5);
        

        // merging two clusters, deleting the old ones and inserting the new one is the synchronisation step
        std::barrier critical_clustering_update(static_cast<std::ptrdiff_t>(num_threads), [&] ()
        {
            // increment id for the new cluster (must be done here at the beginning)
            ++id;

            // compute the final min_id from the min_ids of the worker threads
            size_t min_id = min_ids[0];
            double min_dist = std::numeric_limits<double>::max();
            for (auto candidate_id : min_ids)
            {
                // check if the thread saw any id
                if (candidate_id == none) continue;

                neighbor const & curr = dist.at(candidate_id).top();
                if (curr.dist < min_dist)
                {
                    min_dist = curr.dist;
                    min_id = candidate_id;
                }
            }

            size_t neighbor_id = dist.at(min_id).top().id;

            // merge the two nodes with minimal distance together insert the new node into the clustering
            clustering[id] = {min_id, neighbor_id, std::move(clustering.at(min_id).hll)};
            estimates[id] = clustering.at(id).hll.merge_and_estimate_SIMD(clustering.at(neighbor_id).hll);
            
            // remove old clusters
            dist.erase(min_id);
            dist.erase(neighbor_id);

            // initialize priority queue for the new cluster 
            dist[id];

            // reset the async buffer of the distance matrix for the following updating step
            dist_buffer = dist 
                         | std::views::transform(to_key_and_pointer)
                         | seqan3::views::async_input_buffer(num_threads * 5);
        });

        // initialize priority queues in the distance matrix (sequentially)
        for (auto & [i, i_node] : clustering)
        {
            // empty priority queue for every item in clustering
            dist[i];
        }

        auto worker = [&] (size_t thread_id)
        {
            // minimum distance exclusively for this thread
            double min_dist = std::numeric_limits<double>::max();

            // initialize all the priority queues of the distance matrix
            // while doing that, compute the first min_id
            for (auto i : clustering_buffer)
            {
                for (auto const & [j, j_node] : clustering)
                {
                    // we only want one diagonal of the distance matrix
                    if (i < j)
                    {
                        // this must be a copy, because merging changes the hll sketch
                        hyperloglog temp_hll = sketches.at(i);
                        double const estimate_ij = temp_hll.merge_and_estimate_SIMD(j_node.hll);
                        // Jaccard distance estimate
                        double const distance = 2 - (estimates.at(i) + estimates.at(j)) / estimate_ij;
                        dist.at(i).push({distance, j});
                    }
                }
                if (dist.at(i).empty()) continue;

                // check if the just initialized priority queue contains the minimum value for this thread
                neighbor const & curr = dist.at(i).top();
                if (curr.dist < min_dist)
                {
                    min_dist = curr.dist;
                    min_ids[thread_id] = i;
                }
            }

            // main loop of the clustering
            // keep merging nodes until we have a complete tree
            while (dist.size() > 1)
            {
                // synchronize and perform critical update                
                critical_clustering_update.arrive_and_wait();

                // reset values for the computation of the new minimum
                min_ids[thread_id] = none;
                min_dist = std::numeric_limits<double>::max();
 
                hyperloglog const new_hll = clustering.at(id).hll;
                
                // update distances in dist 
                // while doing that, compute the new min_id
                for (auto [i, prio_q_ptr] : dist_buffer)
                {
                    if (i == id) continue;

                    // this must be a copy, because merge_and_estimate_SIMD() changes the hll
                    hyperloglog temp_hll = new_hll;
                    double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering.at(i).hll);
                    // Jaccard distance estimate
                    double const distance = 2 - (estimates.at(i) + estimates.at(id)) / estimate_ij;
                    prio_q_ptr->push({distance, id});

                    // make sure the closest neighbor is not yet deleted (this is a lazy update)
                    while (!prio_q_ptr->empty() && !dist.contains(prio_q_ptr->top().id))
                    {
                        prio_q_ptr->pop();
                    }

                    if (prio_q_ptr->empty()) continue;

                    // check if the just updated priority queue contains the minimum value for this thread
                    neighbor const & curr = prio_q_ptr->top();
                    if (curr.dist < min_dist)
                    {
                        min_dist = curr.dist;
                        min_ids[thread_id] = i;
                    }
                }
                
            }
        };

        // launch threads with worker
        std::vector<decltype(std::async(std::launch::async, worker, size_t{}))> handles;

        for (size_t i = 0; i < num_threads; ++i)
            handles.emplace_back(std::async(std::launch::async, worker, i));

        // wait for the threads to finish
        for (auto & handle : handles)
            handle.get();

        size_t final_root = dist.begin()->first;

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