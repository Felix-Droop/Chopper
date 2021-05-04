#pragma once

#include <vector>
#include <unordered_map>
#include <queue>
#include <functional>
#include <tuple>
#include <limits>

#include <chopper/union/hyperloglog.hpp>
#include <chopper/union/clustering_node.hpp>

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
    std::unordered_map<size_t, prio_queue> dist;

    //!\brief read-only reference to the clustering tree
    std::unordered_map<size_t, clustering_node> const & clustering;
    
    //!\brief read-only reference to already computes estimates of cardinalities of clusters
    std::unordered_map<size_t, double> const & estimates;

public:
    /*!\brief Distance matrix for a hierarchical clustering algorithm
     *
     * \param[in] clustering_ read-only reference to the clustering tree
     * \param[in] estimates_ read-only reference to already computes estimates of cardinalities
     */
    distance_matrix(std::unordered_map<size_t, clustering_node> const & clustering_,
                             std::unordered_map<size_t, double> const & estimates_) :
        clustering{clustering_},
        estimates{estimates_}
    {}

    /*!\brief Initialize the pairwise jaccard distances in the distance matrix 
     * \param[in] sketches HyperLogLog sketches used to estimate the distances
     */
    void initialize(std::vector<hyperloglog> const & sketches)
    {
        for (auto & [i, i_node] : clustering)
        {
            // we want an empty priority queue for every item
            dist[i];
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
    void update(size_t new_id, size_t old_id_0, size_t old_id_1)
    {
        // remove old clusters
        dist.extract(old_id_0);
        dist.extract(old_id_1);

        // initialize priority queue for the new cluster 
        dist[new_id];

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