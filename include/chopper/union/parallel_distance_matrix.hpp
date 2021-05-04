#pragma once

#include <vector>
#include <unordered_map>
#include <queue>
#include <functional>
#include <tuple>
#include <limits>

#include <chopper/union/hyperloglog.hpp>
#include <chopper/union/clustering_node.hpp>

class parallel_distance_matrix
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
    std::unordered_map<size_t, prio_queue> dist;

    //!\brief read-only reference to the clustering tree
    std::unordered_map<size_t, clustering_node> const & clustering;
    
    //!\brief read-only reference to already computes estimates of cardinalities of clusters
    std::unordered_map<size_t, double> const & estimates;

    //!\brief number of threads to use in parallel operations
    size_t num_threads;

    //!\brief execution handler and thread pool for parallel operations
    // parallel_execution_handler exec_handler 

public:
    /*!\brief Distance matrix for a hierarchical clustering algorithm
     * that computes the expensive operations in parallel
     *
     * \param[in] clustering_ read-only reference to the clustering tree
     * \param[in] estimates_ read-only reference to already computes estimates of cardinalities
     * \param[in] num_threads_ number of threads to use in parallel operations
     */
    parallel_distance_matrix(std::unordered_map<size_t, clustering_node> const & clustering_,
                             std::unordered_map<size_t, double> const & estimates_,
                             size_t num_threads_) :
        clustering{clustering_},
        estimates{estimates_},
        num_threads{num_threads_}
        // exec_handler{num_threads}
    {}

    /*!\brief Initialize the pairwise jaccard distances in the distance matrix 
     * \param[in] sketches HyperLogLog sketches used to estimate the distances
     */
    void initialize(std::vector<hyperloglog> const & sketches)
    {
        // initialize empty priority queues on main thread
        for (auto & [i, i_node] : clustering)
        {
            dist[i];
        }

        // OUTER LOOP COULD BE PARALLELIZED WITH THREADS (LIKELY PART OF THE MAIN BOTTLENECK)
        // ELEMENTS ARE CHANGED -> IS IT SAFE???
        // DEPENDING ON i MORE OR LESS WORK IN THE INNER LOOP -> DYNAMIC SCHEDULING NEEDED

        for (auto & [i, i_node] : clustering)
        {
            // we want an empty priority queue for every item
            for (auto & [j, j_node] : clustering)
            {
                // we only want one diagonal of the distance matrix
                if (i < j)
                {
                    // this must be a copy, because merging changes the hll sketch
                    hyperloglog temp_hll = sketches.at(i);
                    double const estimate_ij = temp_hll.merge_and_estimate_SIMD(sketches.at(j));
                    // Jaccard distance estimate
                    double const distance = 2 - (estimates.at(i) + estimates.at(j)) / estimate_ij;
                    dist.at(i).push({distance, j});
                }
            }
        }
    }

    //!\brief Get the pair of ids that belong to the clusters with the lowest distance to each other
    std::tuple<size_t, size_t> get_min_pair()
    {
        size_t min_id = get_remaining_cluster_id();
        double min_dist = std::numeric_limits<double>::max();

        // THIS LOOP COULD BE PARALLELIZED WITH THREADS (LIKELY NOT THE MAIN BOTTLENECK)
        // ONLY CONST OPERATIONS > SAFE
        // EVERY ITERATION HAS CONSTANT AMOUNT OF WORK -> STATIC (CHUNKED) SCHEDULING POSSIBLE

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

        return std::make_tuple(min_id, dist[min_id].top().id);
    }

    //!\brief Delete the two clusters with the given ids
    void delete_pair(size_t id_0, size_t id_1)
    {
        dist.extract(id_0);
        dist.extract(id_1);
    }

    /*!\brief Initialize the pairwise jaccard distances in the distance matrix 
     * \param[in] sketches HyperLogLog sketches used to estimate the distances
     */
    void update(size_t new_id)
    {
        // initialize priority queue for the new cluster 
        dist[new_id];

        // THIS COULD BE PARALLELIZED WITH THREADS (LIKELY PART OF THE MAIN BOTTLENECK)
        // ELEMENTS ARE CHANGED -> IS IT SAFE???
        // DEPENDING ON i MORE OR LESS WORK IN THE INNER LOOP -> DYNAMIC SCHEDULING NEEDED
        
        // update distances
        for (auto & [i, prio_q] : dist)
        {
            if (i == new_id) continue;

            // this must be a copy, because merge() changes the hll
            hyperloglog temp_hll = clustering.at(new_id).hll;
            double const estimate_ij = temp_hll.merge_and_estimate_SIMD(clustering.at(i).hll);
            // Jaccard distance estimate
            double const distance = 2 - (estimates.at(i) + estimates.at(new_id)) / estimate_ij;
            prio_q.push({distance, new_id});

            // make sure the closest neighbor is not yet deleted (this is a lazy update)
            while (dist.find(prio_q.top().id) == dist.end() && !prio_q.empty())
            {
                prio_q.pop();
            }
        }
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