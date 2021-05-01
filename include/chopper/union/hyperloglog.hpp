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
#include <numeric>
#include <cassert>

#include <xxh3.h>
#include <emmintrin.h>
#include <xmmintrin.h>

/** @class hyperloglog
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
            std::stringstream ss;
            ss << "bit width must be in the range [4,30] and it is " << (int)b;
            throw std::invalid_argument(ss.str().c_str());
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
        alphaMM_float_ = static_cast<float>(alphaMM_);
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
            sum += 1.0 / static_cast<double>(1ull << c);
        }
        double estimate = alphaMM_ / sum; 

        // use linear counting of zeros for small values
        if (estimate <= 2.5 * m_) 
        {
            uint32_t zeros = 0;

            for(size_t i = 0; i < m_; ++i) {
                if (!M_[i]) ++zeros;
            }

            if (zeros != 0u) 
            {
                estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
            }
        } 
        return estimate;
    }

    /**
     * Merges the estimate from 'other' into this object
     * The number of registers in each must be the same.
     * 
     * @param[in] other HyperLogLog instance to be merged
     */    
    void merge(hyperloglog const & other)
    {
        assert(m_ == other.m_);

        for (size_t i = 0; i < m_; ++i)
        {
            if (M_[i] < other.M_[i])
            {
                M_[i] = other.M_[i];
            }
        }
    }

    /**
     * Merges the estimate from 'other' into this object
     * The number of registers in each must be the same.
     * This function is implemented using SSE instructions.
     * 
     * @param[in] other HyperLogLog instance to be merged
     * 
     * @return estimated cardinality of the new merged sketch
     */
    double merge_and_estimate_SSE(hyperloglog const & other) 
    {
        assert(m_ == other.m_);

        // this is safe, because b_ is at least 4 and therefore M_'s size in bits is a power of 2^4 * 8 = 128
        __m128i* it = reinterpret_cast<__m128i*>(&*(M_.begin()));
        const __m128i* other_it = reinterpret_cast< const __m128i*>(&*(other.M_.begin()));
        __m128i* end = reinterpret_cast<__m128i*>(&*(M_.end()));

        __m128 packed_sum = _mm_set_ps1(0.0);

        for (; it != end; ++it, ++other_it)
        {
            // this merges the registers by computing the byte-wise maximum
            *it = _mm_max_epu8(*it, *other_it);
            
            // get pointer to iterate over the single merged registers
            uint8_t* reg_it = reinterpret_cast<uint8_t*>(it);

            // get floats with two to the power of the value in the merged registers
            __m128 values = _mm_set_ps(
                static_cast<float>(((uint64_t)1) << *reg_it),
                static_cast<float>(((uint64_t)1) << *(reg_it + 1)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 2)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 3))
            );

            // compute their reciprocal
            values = _mm_rcp_ps(values); 

            // sum up
            packed_sum = _mm_add_ps(packed_sum, values);

            // and repeat 3 times ...
            reg_it += 4;

            values = _mm_set_ps(
                static_cast<float>(((uint64_t)1) << *reg_it),
                static_cast<float>(((uint64_t)1) << *(reg_it + 1)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 2)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 3))
            );

            values = _mm_rcp_ps(values); 
            packed_sum = _mm_add_ps(packed_sum, values);
            reg_it += 4;

            values = _mm_set_ps(
                static_cast<float>(((uint64_t)1) << *reg_it),
                static_cast<float>(((uint64_t)1) << *(reg_it + 1)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 2)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 3))
            );

            values = _mm_rcp_ps(values); 
            packed_sum = _mm_add_ps(packed_sum, values);
            reg_it += 4;

            values = _mm_set_ps(
                static_cast<float>(((uint64_t)1) << *reg_it),
                static_cast<float>(((uint64_t)1) << *(reg_it + 1)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 2)),
                static_cast<float>(((uint64_t)1) << *(reg_it + 3))
            );

            values = _mm_rcp_ps(values); 
            packed_sum = _mm_add_ps(packed_sum, values);
        }

        // sum up the 4 values in the packed SSE variable
        float sum = 0.0;
        float* sum_it = reinterpret_cast<float*>(&packed_sum);
        sum += *sum_it;
        sum += *(sum_it + 1);
        sum += *(sum_it + 2);
        sum += *(sum_it + 3);

        // compute first estimate
        double estimate = alphaMM_float_ / sum; 

        // use linear counting of zeros for small values
        if (estimate <= 2.5 * m_) 
        {
            uint32_t zeros = 0u;

            for(size_t i = 0; i < m_; ++i) {
                if (!M_[i]) ++zeros;
            }

            if (zeros != 0u) 
            {
                estimate = m_ * std::log(static_cast<double>(m_) / static_cast<double>(zeros));
            }
        }

        return estimate;
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
        std::swap(alphaMM_float_, rhs.alphaMM_float_);
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
            throw std::runtime_error("Failed to dump a HyperLogLog sketch to a file.");
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
           throw std::runtime_error("Failed to restore a HyperLogLog sketch from a file.");
        }       
        swap(tempHLL);
    }

private:
    uint64_t mask_; ///< mask for the rank bits
    double alphaMM_; ///< alpha * m^2
    float alphaMM_float_; ///< alpha * m^2
    uint64_t m_; ///< register size
    uint8_t b_; ///< register bit width
    std::vector<uint8_t> M_; ///< registers
};
