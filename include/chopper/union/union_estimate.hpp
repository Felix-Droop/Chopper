#include <vector>
#include <string>
#include <unordered_map>
#include <queue>
#include <functional>

#include <seqan3/io/sequence_file/input.hpp>

#include <chopper/union/hyperloglog.hpp>

/*!\brief For all intervals of names: estimate the cardinality of the union of all sequences in the files of the interval
    * \param[in] names names of the sequence files
    * \param[in] user_bin_kmer_counts cardinalities of the single sets for the diagonal of the output table
    * \param[in] kmer_size size of k-mers
    * \param[in] sketch_bits The number of bits the HyperLogLog sketch should use to distribute the values into bins
    * \param[out] union_estimates output table
    */
struct union_estimate
{
private:
    //!\brief  traits for data input
    struct input_traits : public seqan3::sequence_file_input_default_traits_dna
    {
        using sequence_alphabet = char;

        using sequence_legal_alphabet = char;

        using sequence_contianer = std::vector<sequence_alphabet>;
    };

    //!\brief node for the clustering
    struct node 
    {
        // children in the tree
        size_t left;
        size_t right;
        // hll sketch of the union if the node is still a root
        hll::HyperLogLog hll;
    };

    //!\brief The file names of the user input. They might be resorted.
    std::vector<std::string> & names;
    
    //!\brief The kmer counts associated with the above files used to pack user bin into technical bins.
    std::vector<size_t> & user_bin_kmer_counts;

    //!\brief HyperLogLog sketches on the sequences from the above files.
    std::vector<hll::HyperLogLog> hlls;

    //!\brief The number of bits the HyperLogLog sketch should use to distribute the values into bins.
    uint8_t sketch_bits;

    //!\brief The k of the k-mers.
    size_t kmer_size;

    //!\brief Whether build_hlls() has been called
    bool built_hlls;

public:
    /*!\brief Constructor of the union_estimate that supplies it with all the important data
     * \param[in] names names of the sequence files
     * \param[in] user_bin_kmer_counts cardinalities of the single sets for the diagonal of the output table
     * \param[in] kmer_size size of k-mers
     * \param[in] sketch_bits The number of bits the HyperLogLog sketch should use to distribute the values into bins
     */
    union_estimate(std::vector<std::string> & names_,
                   std::vector<size_t> & user_bin_kmer_counts_,
                   size_t kmer_size_, 
                   uint8_t sketch_bits_) :
        names{names_},
        user_bin_kmer_counts{user_bin_kmer_counts_},
        hlls{},
        sketch_bits{sketch_bits_},
        kmer_size{kmer_size_},
        built_hlls{false}
    {
    }

    //!\brief Construct HyperLogLog sketches from the sequences given by the member names
    void build_hlls()
    {
        using sequence_file_type = seqan3::sequence_file_input<input_traits, seqan3::fields<seqan3::field::seq>>;

        // build hll sketches from the sequences of all files
        for (auto && filename : names)
        {
            sequence_file_type seq_file{filename};
            hll::HyperLogLog & curr_hll = hlls.emplace_back(sketch_bits);

            // put every sequence in this file into the sketch
            for (auto && [seq] : seq_file)
            {
                // we have to go C-style here for the HyperLogLog Interface
                char* c_seq_it = &*seq.begin();
                char* end = c_seq_it + seq.size();

                while(c_seq_it + kmer_size <= end)
                {
                    curr_hll.add(c_seq_it, kmer_size);
                    ++c_seq_it;
                }
            }
        }

        built_hlls = true;
    }

    //!\brief Reorder names, user_bin_kmer_counts and hlls such that similar bins are close to each other
    void resort_bins()
    {
        std::vector<size_t> permutation;

        // TODO do some small overlap, e.g. the last cluster from before

        size_t first = 0;
        size_t last = 1;
        double const max_ratio = 0.5;
        while(first < names.size())
        {
            // size difference is too large or sequence is over -> do the clustering
            if (user_bin_kmer_counts[first] * max_ratio > user_bin_kmer_counts[last] || last == names.size())
            {
                cluster_bins(permutation, first, last);
                first = last;
            }
            ++last;
        }

        // apply permutation to names, user_bin_kmar_counts and hlls
        for (size_t i = 0; i < permutation.size(); ++i)
        {
            size_t current = i;
            while (i != permutation[current])
            {
                size_t next = permutation[current];
                std::swap(names[current], names[next]);
                std::swap(user_bin_kmer_counts[current], user_bin_kmer_counts[next]);
                std::swap(hlls[current], hlls[next]);
                permutation[current] = current;
                current = next;
            }
            permutation[current] = current;
        }
    }

    /*!\brief For all intervals of names: estimate the cardinality of the union of all sequences in the files of the interval
     * \param[out] union_estimates output table
     */
    void estimate_unions(std::vector<std::vector<size_t>> & union_estimates)
    {
        if (!built_hlls) throw std::runtime_error{"Must build hlls first."};

        // union_estimates[i][j] will be the union of the interval i, ..., i+j
        union_estimates.clear();
        size_t const n = user_bin_kmer_counts.size();

        for (size_t i = 0; i < n; ++i)
        {
            std::vector<size_t> & curr_vec = union_estimates.emplace_back();

            // for the single set we have the exact value - no need for the hll estimate here
            curr_vec.push_back(user_bin_kmer_counts[i]);

            for (size_t j = i + 1; j < n; ++j)
            {
                // merge next sketch into the current union
                hlls[i].merge(hlls[j]);

                curr_vec.push_back(static_cast<size_t>(hlls[i].estimate()));
            }
        }
    }

private:
    /*!\brief Perform an agglomerative clustering variant on the index range [first:last)
     * \param[in] start
     * \param[in] last
     * \param[out] permutation append the new order to this
     */
    void cluster_bins(std::vector<size_t> & permutation,
                      size_t first,
                      size_t last)
    {
        struct neighbor
        {
            double dist;
            size_t id;

            bool operator>(neighbor const & other) const
            {
                return dist > other.dist;
            }
        };

        // (possible improvement: these two could be just std::vectors with an index shift by 'start')
        std::unordered_map<size_t, node> clustering;
        std::unordered_map<size_t, double> estimates;

        // store distances in a heap to improve running time
        using prio_queue = std::priority_queue<neighbor, std::vector<neighbor>, std::greater<neighbor>>;
        std::unordered_map<size_t, prio_queue> dist;

        size_t const none = std::numeric_limits<size_t>::max();

        // initialize clustering and estimates
        for (size_t i = first; i < last; ++i)
        {
            clustering[i] = {none, none, hlls[i]};
            estimates[i] = hlls[i].estimate();
        }

        // initialize dist
        for (auto & [i, i_node] : clustering)
        {
            for (auto & [j, j_node] : clustering)
            {
                // we only want one diagonal of the distance matrix
                if (i < j)
                {
                    // this must be a copy, because merge() changes the hll
                    hll::HyperLogLog temp_hll = hlls[i];
                    temp_hll.merge(hlls[j]);
                    double const estimate_ij = temp_hll.estimate();
                    // Jaccard distance estimate
                    double const distance = 2 - (estimates[i] + estimates[j]) / estimate_ij;
                    dist[i].push({distance, j});
                }
                // we want an empty priority queue for every item
                dist[i];
            }
        }

        // ids of new nodes start at last's value
        size_t id = last;

        // keep merging nodes until we have a complete tree
        while (dist.size() > 1)
        {
            size_t min_id = none;
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

            size_t neighbor_id = dist[min_id].top().id;

            // merge the two nodes with minimal distance together and insert the new node into the clustering
            clustering[id] = {min_id, neighbor_id, std::move(clustering[min_id].hll)};
            node & new_root = clustering[id];
            new_root.hll.merge(clustering[neighbor_id].hll);

            // insert the new node into dist and update estimates
            estimates[id] = new_root.hll.estimate();
            dist[id];

            // delete them from dist
            dist.extract(min_id);
            dist.extract(neighbor_id);

            // update distances
            for (auto & [i, prio_q] : dist)
            {
                if (i == id) continue;

                // this must be a copy, because merge() changes the hll
                hll::HyperLogLog temp_hll = new_root.hll;
                temp_hll.merge(clustering[i].hll);
                double const estimate_ij = temp_hll.estimate();
                // Jaccard distance estimate
                double const distance = 2 - (estimates[i] + estimates[id]) / estimate_ij;
                prio_q.push({distance, id});

                // make sure the closest neighbor is not yet deleted (this is kind of a lazy update)
                while (dist.find(prio_q.top().id) == dist.end() && !prio_q.empty())
                {
                    prio_q.pop();
                }
            }

            ++id;
        }

        // traceback into permutation
        size_t final_root = dist.begin()->first;
        trace(clustering, permutation, final_root);
    }

    /*!\brief Do a recursive traceback to find the order of leaves in the clustering tree 
     * \param[in] clustering the tree to do the traceback on
     * \param[out] permutation append the new order to this
     * \param[in] id the id of the current node
     */
    void trace(std::unordered_map<size_t, node> const & clustering,
               std::vector<size_t> & permutation,
               size_t id)
    {
        if (clustering.at(id).left == std::numeric_limits<size_t>::max())
        {
            permutation.push_back(id);
            return;
        }

        trace(clustering, permutation, clustering.at(id).left);
        trace(clustering, permutation, clustering.at(id).right);
    }
};