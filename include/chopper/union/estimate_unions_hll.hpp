#include <vector>
#include <string>
#include <seqan3/io/sequence_file/input.hpp>

#include <chopper/union/hyperloglog.hpp>

struct mytraits_union_hll : public seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = char;

    using sequence_legal_alphabet = char;

    using sequence_contianer = std::vector<sequence_alphabet>;
};

/*!\brief For all intervals of names: estimate the cardinality of the union of all sequences in the files of the interval
    * \param[in] names names of the sequence files
    * \param[in] user_bin_kmer_counts cardinalities of the single sets for the diagonal of the output table
    * \param[in] kmer_size size of k-mers
    * \param[in] sketch_bits The number of bits the HyperLogLog sketch should use to distribute the values into bins
    * \param[out] union_estimates output table
    */
void estimate_unions_hll(std::vector<std::string> const & names,
                         std::vector<size_t> const & user_bin_kmer_counts,
                         size_t kmer_size, 
                         uint8_t sketch_bits,
                         std::vector<std::vector<size_t>> & union_estimates)
{
    using sequence_file_type = seqan3::sequence_file_input<mytraits_union_hll, seqan3::fields<seqan3::field::seq>>;

    std::vector<hll::HyperLogLog> hlls;

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

    // union_estimates[i][j] will be the union of the interval i, ..., i+j
    union_estimates.clear();
    size_t const n = names.size();

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

