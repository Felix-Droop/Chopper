#pragma once

#include <vector>
#include <queue>
#include <functional>
#include <tuple>
#include <limits>

#include <seqan3/range/views/async_input_buffer.hpp>
#include <seqan3/std/ranges>

#include <chopper/union/hyperloglog.hpp>
#include <chopper/union/clustering_node.hpp>

#include <robin_hood.h>

class distance_matrix
{
private:
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

    //!\brief read-only reference to the clustering tree
    robin_hood::unordered_map<size_t, clustering_node> const & clustering;
    
    //!\brief read-only reference to already computes estimates of cardinalities of clusters
    robin_hood::unordered_map<size_t, double> const & estimates;

    //!\brief The number of threads to use to compute merged HLL sketches.
    size_t const num_threads;

public:
    /*!\brief Distance matrix for a hierarchical clustering algorithm
     *
     * \param[in] clustering_ read-only reference to the clustering tree
     * \param[in] estimates_ read-only reference to already computes estimates of cardinalities
     */
    distance_matrix(robin_hood::unordered_map<size_t, clustering_node> const & clustering_,
                    robin_hood::unordered_map<size_t, double> const & estimates_,
                    size_t const num_threads_) :
        clustering{clustering_},
        estimates{estimates_},
        num_threads{num_threads_}
    {}

    /*!\brief Initialize the pairwise jaccard distances in the distance matrix 
     * \param[in] sketches HyperLogLog sketches used to estimate the distances
     */
    void initialize()
    {
        for (auto & [i, i_node] : clustering)
        {
            // we want an empty priority queue for every item
            dist[i];
        }

        // only the keys are needed when iterating through the clustering
        auto to_key_and_pointer = [] (robin_hood::unordered_map<size_t, clustering_node>::value_type const & v)
        {
            return std::make_tuple(v.first, &v.second);
        };

        auto clustering_buffer = clustering 
                               | std::views::transform(to_key_and_pointer)
                               | seqan3::views::async_input_buffer(num_threads * 5);

        auto worker = [this, &clustering_buffer] ()
        {
            for (auto [i, i_node_ptr] : clustering_buffer)
            {
                for (auto & [j, j_node] : clustering)
                {
                    // we only want one diagonal of the distance matrix
                    if (i < j)
                    {
                        // this must be a copy, because merging changes the hll sketch
                        hyperloglog temp_hll = i_node_ptr->hll;
                        double const estimate_ij = temp_hll.merge_and_estimate_SIMD(j_node.hll);
                        // Jaccard distance estimate
                        double const distance = 2 - (estimates.at(i) + estimates.at(j)) / estimate_ij;
                        dist.at(i).push({distance, j});
                    }
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

    //!\brief Get the pair of ids that belong to the clusters with the lowest distance to each other
    std::tuple<size_t, size_t> get_min_pair()
    {
        size_t min_id = get_remaining_cluster_id();
        double min_dist = std::numeric_limits<double>::max();

        // maybe do at the end of update()
        // find the two nodes with the minimal distance
        for (auto & [i, prio_q] : dist)
        {
            if (prio_q.empty()) continue;

            neighbor const & curr = prio_q.top();
            if (curr.dist < min_dist)
            {
                min_dist = curr.dist;
                min_id = i;
            }
        }

        return std::make_tuple(min_id, dist.at(min_id).top().id);
    }

    /*!\brief Update the distance matrix.
     * 
     * Delete the two old clusters with the given ids.
     * Insert a new cluster with the given ids.
     * Compute distances from previously existing clusters to the new one.
     * 
     * \param[in] new_id id of the new cluster
     * \param[in] old_id_0 id of the first old cluster
     * \param[in] old_id_1 id of the second old cluster
     */
    void update(size_t const new_id, size_t const old_id_0, size_t const old_id_1)
    {
        // remove old clusters
        dist.erase(old_id_0);
        dist.erase(old_id_1);

        // initialize priority queue for the new cluster 
        dist[new_id];

        // only the keys are needed when iterating through the clustering
        auto to_key_and_pointer = [] (robin_hood::unordered_map<size_t, prio_queue>::value_type & v)
        {
            return std::make_tuple(v.first, &v.second);
        };

        auto dist_buffer = dist 
                         | std::views::transform(to_key_and_pointer)
                         | seqan3::views::async_input_buffer(num_threads * 5);
                         
        auto worker = [this, new_id, &dist_buffer] () {
            hyperloglog const new_hll = clustering.at(new_id).hll;
            
            for (auto [i, prio_q_ptr] : dist_buffer)
            {
                if (i == new_id) continue;

                // this must be a copy, because merge() changes the hll
                hyperloglog temp_hll = new_hll;
                double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering.at(i).hll);
                // Jaccard distance estimate
                double const distance = 2 - (estimates.at(i) + estimates.at(new_id)) / estimate_ij;
                prio_q_ptr->push({distance, new_id});

                // make sure the closest neighbor is not yet deleted (this is a lazy update)
                while (!prio_q_ptr->empty() && !dist.contains(prio_q_ptr->top().id))
                {
                    prio_q_ptr->pop();
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

    //!\brief Number of clusters in the matrix (NOT number of distances)
    size_t size()
    {
        return dist.size();
    }

    /*!\brief Get the id of some cluster in the matrix. 
     * If there is only one cluster left, get that clusters' id.
     */
    size_t get_remaining_cluster_id()
    {
        return dist.begin()->first;
    }
};