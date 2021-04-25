#pragma once

/**
 * @file hyperloglog.hpp
 * @brief HyperLogLog cardinality estimator
 * @date Created 2013/3/20
 * @author Hideaki Ohno
 * 
 * Copied to this location from github by Felix Droop - Jan 5 2021
 * Modified a lot for a bugfix, improvements and functional changes (64 bit hashes)
 */

#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <iostream> 
#include <xxh3.h>
#include <numeric>

/** @class HyperLogLog
 *  @brief Implement of 'HyperLogLog' estimate cardinality algorithm
 */
class hyperloglog 
{
public:

    /**
     * Constructor
     *
     * @param[in] b bit width (register size will be 2 to the b power).
     *            This value must be in the range[4,30].Default value is 4.
     *
     * @exception std::invalid_argument the argument is out of range.
     */
    hyperloglog(uint8_t b = 4) :
            m_(1 << b), b_(b), M_(m_, 0) {

        if (b < 4 || 30 < b) 
        {
            throw std::invalid_argument("bit width must be in the range [4,30]");
        }
        M_.shrink_to_fit();
        double alpha;
        switch (m_) {
            case 16:
                alpha = 0.673;
                break;
            case 32:
                alpha = 0.697;
                break;
            case 64:
                alpha = 0.709;
                break;
            default:
                alpha = 0.7213 / (1.0 + 1.079 / m_);
                break;
        }
        alphaMM_ = alpha * m_ * m_;
        // 64 bits where the last b are ones and the rest zeroes
        mask_ = (1 << b) - 1;
    }

    /**
     * Adds element to the estimator
     *
     * @param[in] str string to add
     * @param[in] len length of string
     */
    void add(const char* str, uint64_t len) 
    {
        uint64_t hash = XXH3_64bits(str, len);
        // the first b_ bits are used to distribute the leading zero counts along M_
        uint64_t index = hash >> (64 - b_);
        // WARNING: __builtin_clzl() only works with g++ and clang
        // the bitwise-or with mask_ assures that we get at most 64 - b_ as value. 
        // Otherwise the count for hash = 0 would be 64
        uint8_t rank = __builtin_clzl((hash << b_) | mask_) + 1;
        M_[index] = std::max(rank, M_[index]);
    }

    /**
     * Estimates cardinality value.
     *
     * @return Estimated cardinality value.
     */
    double estimate() const 
    {
        // compute indicator formula
        double sum = 0.0;
        for (uint8_t c : M_)
        {
            sum += 1.0 / (1ull << c);
        }
        double estimate = alphaMM_ / sum; 

        // use linear counting of zeros for small values
        if (estimate <= 2.5 * m_) 
        {
            uint64_t zeros = std::count(M_.cbegin(), M_.cend(), 0);
            if (zeros != 0ull) 
            {
                estimate = m_ * std::log(static_cast<double>(m_) / zeros);
            }
        } 
        return estimate;
    }

    /**
     * Merges the estimate from 'other' into this object
     * The number of registers in each must be the same.
     *
     * @param[in] other HyperLogLog instance to be merged
     * 
     * @exception std::invalid_argument number of registers doesn't match.
     */
    void merge(hyperloglog const & other) 
    {
        if (m_ != other.m_) 
        {
            std::stringstream ss;
            ss << "number of registers doesn't match: " << m_ << " != " << other.m_;
            throw std::invalid_argument(ss.str().c_str());
        }

        std::transform(M_.begin(), M_.end(), other.M_.begin(), 
            M_.begin(), [] (uint8_t const x, uint8_t const y) 
            {
                return std::max(x, y);
            });
    }

    /**
     * Clears all internal registers.
     */
    void clear() 
    {
        std::fill(M_.begin(), M_.end(), 0);
    }

    /**
     * Returns size of register.
     *
     * @return Register size
     */
    uint64_t registerSize() const 
    {
        return m_;
    }

    /**
     * Exchanges the content of the instance
     *
     * @param[in,out] rhs Another HyperLogLog instance
     */
    void swap(hyperloglog& rhs) 
    {
        std::swap(mask_, rhs.mask_);
        std::swap(alphaMM_, rhs.alphaMM_);
        std::swap(m_, rhs.m_);
        std::swap(b_, rhs.b_);
        M_.swap(rhs.M_);       
    }

    /**
     * Dump the current status to a stream
     *
     * @param[out] os The output stream where the data is saved
     *
     * @exception std::runtime_error When failed to dump.
     */
    void dump(std::ostream& os) const 
    {
        os.write((char*)&b_, sizeof(b_));
        os.write((char*)&M_[0], sizeof(M_[0]) * M_.size());
        os.flush();
        if (os.fail())
        {
            throw std::runtime_error("Failed to dump");
        }
    }

    /**
     * Restore the status from a stream
     * 
     * @param[in] is The input stream where the status is saved
     *
     * @exception std::runtime_error When failed to restore.
     */
    void restore(std::istream& is) 
    {
        uint8_t b = 0;
        is.read((char*)&b, sizeof(b));
        hyperloglog tempHLL(b);
        is.read((char*)&(tempHLL.M_[0]), sizeof(M_[0]) * tempHLL.m_);
        if (is.fail())
        {
           throw std::runtime_error("Failed to restore");
        }       
        swap(tempHLL);
    }

private:
    uint64_t mask_; ///< mask for the rank bits
    double alphaMM_; ///< alpha * m^2
    uint64_t m_; ///< register size
    uint8_t b_; ///< register bit width
    std::vector<uint8_t> M_; ///< registers
};
