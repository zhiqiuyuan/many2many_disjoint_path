#ifndef _PARTITION_H
#define _PARTITION_H

#include "BatchSet.h"

template <typename BatchSetType>
BatchSizeType *uniform_split(BatchSizeType total_vp_num, BatchSizeType &batch_num)
{
    BatchSizeType bsz = BatchSetType::batch_sz();
    batch_num = total_vp_num / bsz;
    if (total_vp_num % bsz)
        ++batch_num;
    BatchSizeType *ranges = new BatchSizeType[batch_num * 2];
    BatchSizeType idx = 0;
    BatchSizeType nextb;
    for (BatchSizeType b = 0; b < total_vp_num; b = nextb)
    {
        ranges[idx++] = b;
        nextb = b + bsz;
        ranges[idx++] = std::min(total_vp_num, nextb);
    }
    return ranges;
}

// random
template <typename BatchSetType>
BatchSizeType *uniform(Graph &g, VertexIdType *&vps, BatchSizeType total_vp_num, BatchSizeType &batch_num)
{
    std::vector<BatchSizeType> vp_idxs(total_vp_num);
    for (BatchSizeType i = 0; i < total_vp_num; ++i)
        vp_idxs[i] = i;
    std::random_shuffle(vp_idxs.begin(), vp_idxs.end());
    VertexIdType *new_vps = new VertexIdType[2 * total_vp_num];
    for (BatchSizeType i = 0, pos; i < total_vp_num; ++i)
    {
        pos = vp_idxs[i];
        new_vps[i * 2] = vps[2 * pos];
        new_vps[i * 2 + 1] = vps[2 * pos + 1];
    }
    delete[] vps;
    vps = new_vps;

    return uniform_split<BatchSetType>(total_vp_num, batch_num);
}

// deg(s)+deg(t) descending, then diff(deg(s),deg(t)) ascending
template <typename BatchSetType>
BatchSizeType *degreeDescending(Graph &g, VertexIdType *&vps, BatchSizeType total_vp_num, VertexIdType &batch_num)
{
    struct SortVertexPair
    {
        VertexIdType v1, v2;
        EdgeNumType key1, key2;
        SortVertexPair(VertexIdType v1, VertexIdType v2, EdgeNumType key1, EdgeNumType key2) : v1(v1), v2(v2), key1(key1), key2(key2)
        {
        }
    };
    std::vector<SortVertexPair> arr;
    EdgeNumType *offset = g.get_offset(), key1, key2, diff;
    for (BatchSizeType i = 0; i < total_vp_num; ++i)
    {
        VertexIdType v1 = vps[2 * i], v2 = vps[2 * i + 1];
        key1 = offset[v1 + 1] - offset[v1];
        key2 = offset[v2 + 1] - offset[v2];
        if (key1 < key2)
            diff = key2 - key1;
        else
            diff = key1 - key2;
        // key1: deg(s)+deg(t)
        // key2: diff(deg(s),deg(t))
        arr.emplace_back(v1, v2, key1 + key2, diff);
    }
    auto cmpSortVertexPair = [](const SortVertexPair &l, const SortVertexPair &r) -> bool
    {
        if (l.key1 != r.key1)
            return l.key1 > r.key1;
        return l.key2 < r.key2;
    };
    std::sort(arr.begin(), arr.end(), cmpSortVertexPair);

    BatchSizeType i = 0;
    for (auto &p : arr)
    {
        vps[i++] = p.v1;
        vps[i++] = p.v2;
    }

    return uniform_split<BatchSetType>(total_vp_num, batch_num);
}

template <typename BatchSetType>
BatchSizeType *regroup(VertexIdType *&vps, BatchSizeType &total_vp_num,
                       BatchSizeType &batch_num,
                       const std::vector<std::vector<std::vector<VertexIdType>>> &p1s,
                       std::vector<BatchSizeType> &new2old_posInVPs,
                       size_t should_have_path_num)
{
    std::vector<BatchSizeType> old_new2old_posInVPs = new2old_posInVPs;
    std::vector<VertexIdType> new_vps_vec;
    new2old_posInVPs.clear();
    BatchSizeType jed = total_vp_num;
    if (old_new2old_posInVPs.empty() == 0)
        jed = (BatchSizeType)old_new2old_posInVPs.size();
    for (BatchSizeType j = 0, i; j < jed; ++j)
    {
        if (old_new2old_posInVPs.empty())
            i = j;
        else
            i = old_new2old_posInVPs[j];
        if (p1s[i].size() == should_have_path_num && p1s[i].back().size() >= 2)
        {
            new2old_posInVPs.push_back(i);
            new_vps_vec.push_back(vps[2 * j]);
            new_vps_vec.push_back(vps[2 * j + 1]);
        }
    }
    total_vp_num = new2old_posInVPs.size();
    VertexIdType *new_vps = new VertexIdType[total_vp_num * 2];
    for (size_t i = 0, ied = new_vps_vec.size(); i < ied; ++i)
        new_vps[i] = new_vps_vec[i];
    delete[] vps;
    vps = new_vps;

    return uniform_split<BatchSetType>(total_vp_num, batch_num);
}

#endif // _PARTITION_H