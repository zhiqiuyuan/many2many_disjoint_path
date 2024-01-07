#ifndef _BATCHSET_H
#define _BATCHSET_H

#include "tools.h"

inline uint16_t _mm_extract_epi16_var(__m128i A, uint16_t imm)
{
    switch (imm)
    {
    case 0:
        return _mm_extract_epi16(A, 0);
    case 1:
        return _mm_extract_epi16(A, 1);
    case 2:
        return _mm_extract_epi16(A, 2);
    case 3:
        return _mm_extract_epi16(A, 3);
    case 4:
        return _mm_extract_epi16(A, 4);
    case 5:
        return _mm_extract_epi16(A, 5);
    case 6:
        return _mm_extract_epi16(A, 6);
    case 7:
        return _mm_extract_epi16(A, 7);
    default:
        return 0;
    }
}

class M128iOnesIdx;
class M128iWrapper
{
    static __m128i masks[129];

public:
    __m128i reg;
    typedef M128iOnesIdx OnesIdx;

    static void init_masks()
    {
        M128iWrapper tmp;
        tmp.set_all();
        for (int i = 128, next_i; i > 0; i = next_i)
        {
            masks[i] = tmp.reg;
            next_i = i - 1;
            tmp.reset_idx(next_i);
        }
        masks[0] = _mm_setzero_si128();
    }

    M128iWrapper() : reg(_mm_setzero_si128())
    {
    }
    M128iWrapper(const M128iWrapper &other) : reg(other.reg)
    {
    }
    M128iWrapper(M128iWrapper &&other) : reg(other.reg)
    {
    }
    M128iWrapper(__m128i other) : reg(other)
    {
    }
    // (a&(~b))&c
    M128iWrapper(const M128iWrapper &a, const M128iWrapper &b, const M128iWrapper &c)
    {
        reg = _mm_and_si128(_mm_andnot_si128(b.reg, a.reg), c.reg);
    }
    // a&b
    M128iWrapper(const M128iWrapper &a, const M128iWrapper &b)
    {
        reg = _mm_and_si128(b.reg, a.reg);
    }
    // a&(~b)
    M128iWrapper(const M128iWrapper &a, const M128iWrapper &b, int)
    {
        reg = _mm_andnot_si128(b.reg, a.reg);
    }
    void operator=(const M128iWrapper &other)
    {
        reg = other.reg;
    }
    void operator=(M128iWrapper &&other)
    {
        reg = other.reg;
    }

    void set_idx(BatchSizeType idx)
    {
        assert(idx < 128 && idx >= 0);
        uint16_t j = idx % 16;
        ((uint16_t *)&reg)[idx / 16] |= (1u << j);
    }
    void reset_idx(BatchSizeType idx)
    {
        assert(idx < 128 && idx >= 0);
        uint16_t j = idx % 16;
        ((uint16_t *)&reg)[idx / 16] &= (~(1u << j));
    }
    bool idx_is_one(BatchSizeType idx) const
    {
        assert(idx < 128 && idx >= 0);
        uint16_t j = idx % 16;
        uint16_t seg = _mm_extract_epi16_var(reg, idx / 16);
        return (seg & (1u << j));
    }
    void set_all() { reg = _mm_set1_epi32(-1u); }
    void zero_all() { reg = _mm_setzero_si128(); }
    bool __attribute__((hot)) is_all_zero(BatchSizeType len) const
    {
        __m128i zero = _mm_setzero_si128();
        __m128i masked_reg;
        if (len > 0 && len < 128)
            masked_reg = _mm_and_si128(reg, masks[len]);
        else
            masked_reg = reg;

        if (_mm_mask_cmp_epi32_mask((__mmask8)-1u, masked_reg, zero, 4) == 0)
            return 1;
        return 0;
    }

    M128iWrapper operator&(const M128iWrapper &r) const { return _mm_and_si128(r.reg, reg); }
    void operator&=(const M128iWrapper &r) { reg = _mm_and_si128(r.reg, reg); }
    M128iWrapper operator|(const M128iWrapper &r) const { return _mm_or_si128(r.reg, reg); }
    void operator|=(const M128iWrapper &r) { reg = _mm_or_si128(r.reg, reg); }
    M128iWrapper operator~() const { return _mm_xor_si128(reg, _mm_set1_epi32(-1u)); }
    // this & ~r
    void andnot(const M128iWrapper &r)
    {
        reg = _mm_andnot_si128(r.reg, reg);
    }
    void notall() { reg = _mm_xor_si128(reg, _mm_set1_epi32(-1u)); }

    void print_bit() const
    {
        std::cout << "left->right,up->bottom: low addr->high addr" << std::endl;
        // 128/16=8
        for (int8_t i = 0; i < 8; ++i)
        {
            uint16_t seg = _mm_extract_epi16_var(reg, i);
            std::cout << std::setw(3) << i * 16 << " ";
            for (int8_t j = 0; j < 16; ++j)
            {
                std::cout << (seg & (uint16_t)1);
                seg = (seg >> 1);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    static BatchSizeType batch_sz() { return 128; }
};
class M128iOnesIdx
{
    __m128i reg;
    uint32_t bit1pos[4][32];
    uint32_t cnt[4] = {0};
    uint32_t __i, __j;
    uint32_t sz;

public:
    M128iOnesIdx(const M128iWrapper &other, BatchSizeType len) : reg(other.reg)
    {
        __mmask8 len_msk;
        {
            __mmask8 msk1cnt = len / 32;
            if (len % 32)
                ++msk1cnt;
            len_msk = (((__mmask8)1) << msk1cnt) - (__mmask8)1;
        }

        __m128i step = _mm_set1_epi32(1);
        __m128i all32 = _mm_set1_epi32(32);
        __m128i idx_upper = _mm_set_epi32(127, 95, 63, 31);

        __m128i copy_reg = reg;

        while (1)
        {
            __m128i lz = _mm_lzcnt_epi32(copy_reg);

            __mmask8 msk = _mm_mask_cmp_epi32_mask((__mmask8)-1u, lz, all32, 4);
            msk &= len_msk;
            if (msk == 0)
                break;
            __m128i four_1s = _mm_sub_epi32(idx_upper, lz);

            __mmask8 k = 1;
            uint32_t pos;
            for (int i = 0; i < 4; ++i, k = (k << 1))
            {
                pos = _mm_extract_epi16_var(four_1s, i << 1);
                if ((msk & k) && pos < len)
                {
                    bit1pos[i][cnt[i]++] = pos;
                }
            }

            __m128i lshift = _mm_add_epi32(lz, step);
            idx_upper = _mm_sub_epi32(idx_upper, lshift);
            for (k = 1; k < 9; k = (k << 1)) // 1 2(10) 4(100) 8(1000)
            {
                if (k & msk)
                {
                    copy_reg = _mm_mask_sllv_epi32(copy_reg, k, copy_reg, lshift);
                }
            }
        }

        init();
        sz = 0;
        for (int i = 0; i < 4; ++i)
        {
            sz += cnt[i];
        }
    }

    void init()
    {
        __i = 0;
        while (__i < 4 && cnt[__i] == 0)
            ++__i;
        if (__i < 4)
            __j = cnt[__i] - 1;
    }

    BatchSizeType next()
    {
        if (__i >= 4 || (__i < 4 && __j >= cnt[__i]))
            return (BatchSizeType)-1ull;
        uint32_t i0 = __i, j0 = __j;
        if (__j == 0)
        {
            ++__i;
            while (__i < 4 && cnt[__i] == 0)
                ++__i;
            if (__i < 4)
                __j = cnt[__i] - 1;
        }
        else
            --__j;
        return bit1pos[i0][j0];
    }
    BatchSizeType size() { return sz; }
};

template <typename UInt>
class UIntOnesIdx;
template <typename UInt>
class BuiltinUIntWrapper
{
    static UInt masks[sizeof(UInt) * 8 + 1];

public:
    UInt reg;
    typedef UIntOnesIdx<UInt> OnesIdx;
    friend class UIntOnesIdx<UInt>;

    static void init_masks()
    {
        UInt msk = (UInt)-1ull;
        for (int i = sizeof(UInt) * 8; i > 0; --i)
        {
            masks[i] = msk;
            msk = (msk >> 1);
        }
        masks[0] = 0;
    }

    BuiltinUIntWrapper() : reg((UInt)0){};
    BuiltinUIntWrapper(const BuiltinUIntWrapper &other) : reg(other.reg)
    {
    }
    BuiltinUIntWrapper(BuiltinUIntWrapper &&other) : reg(other.reg)
    {
    }
    BuiltinUIntWrapper(UInt other) : reg(other)
    {
    }
    // (a&(~b))&c
    BuiltinUIntWrapper(const BuiltinUIntWrapper &a, const BuiltinUIntWrapper &b, const BuiltinUIntWrapper &c)
    {
        reg = ((a.reg) & (~(b.reg))) & (c.reg);
    }
    // a&b
    BuiltinUIntWrapper(const BuiltinUIntWrapper &a, const BuiltinUIntWrapper &b)
    {
        reg = (a.reg) & (b.reg);
    }
    // a&(~b)
    BuiltinUIntWrapper(const BuiltinUIntWrapper &a, const BuiltinUIntWrapper &b, int)
    {
        reg = (a.reg) & (~b.reg);
    }
    ~BuiltinUIntWrapper()
    {
    }
    void operator=(const BuiltinUIntWrapper &other)
    {
        reg = other.reg;
    }
    void operator=(BuiltinUIntWrapper &&other)
    {
        reg = other.reg;
    }

    void set_idx(BatchSizeType idx)
    {
        assert(idx < sizeof(UInt) * 8 && idx >= 0);
        reg |= (((UInt)1) << idx);
    }
    void reset_idx(BatchSizeType idx)
    {
        assert(idx < sizeof(UInt) * 8 && idx >= 0);
        reg &= (~(((UInt)1) << idx));
    }
    bool idx_is_one(BatchSizeType idx) const
    {
        assert(idx < sizeof(UInt) * 8 && idx >= 0);
        return (reg & (((UInt)1) << idx));
    }
    void set_all()
    {
        reg = (UInt)0;
        reg = ~reg;
    }
    void zero_all() { reg = (UInt)0; }
    bool __attribute__((hot)) is_all_zero(BatchSizeType len) const
    {
        if (len > 0 && len < sizeof(UInt) * 8)
            return (reg & masks[len]) == 0;
        return reg == (UInt)0;
    }

    BuiltinUIntWrapper operator&(const BuiltinUIntWrapper &r) const { return reg & r.reg; }
    void operator&=(const BuiltinUIntWrapper &r) { reg &= r.reg; }
    BuiltinUIntWrapper operator|(const BuiltinUIntWrapper &r) const { return (r.reg) | reg; }
    void operator|=(const BuiltinUIntWrapper &r) { reg |= (r.reg); }
    BuiltinUIntWrapper operator~() const { return ~reg; }
    void andnot(const BuiltinUIntWrapper &r) { reg = reg & (~(r.reg)); } // this & ~r
    void notall() { reg = ~reg; }

    void print_bit() const
    {
        std::cout << "left->right,up->bottom: low addr->high addr" << std::endl;
        BatchSizeType ed = sizeof(UInt) * 8;
        for (BatchSizeType i = 0; i < ed; i += 16)
        {
            uint16_t seg = (reg >> i);
            std::cout << std::setw(3) << i << " ";
            for (int8_t j = 0; j < 16; ++j)
            {
                std::cout << (seg & (uint16_t)1);
                seg = (seg >> 1);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    static BatchSizeType batch_sz() { return sizeof(UInt) * 8; }
};
template <typename UInt>
class UIntOnesIdx
{
    UInt ori_reg, reg;
    uint32_t has_shift_len;

public:
    UIntOnesIdx(const BuiltinUIntWrapper<UInt> &other, BatchSizeType len)
    {
        if (len > 0 && len < sizeof(UInt) * 8)
            ori_reg = other.reg & (BuiltinUIntWrapper<UInt>::masks[len]);
        else
            ori_reg = other.reg;
        init();
    }
    void init()
    {
        has_shift_len = 0;
        reg = ori_reg;
    }
    BatchSizeType next()
    {
        if (reg == 0)
            return (BatchSizeType)-1ull;
        BatchSizeType posp1 = __builtin_ffsll(reg);
        if (posp1 >= sizeof(UInt) * 8)
            reg = 0;
        else
            reg = (reg >> posp1);
        has_shift_len += posp1;
        return has_shift_len - 1;
    }
    BatchSizeType size()
    {
        return __builtin_popcountll(ori_reg);
    }
};

#endif //_BATCHSET_H