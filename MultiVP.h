#ifndef _MULTIVP_H
#define _MULTIVP_H

#include "BatchSet.h"
#include "Baseline.h"
#include "Partition.h"

template <typename BatchSetType>
class OneBatchPath2
{
    std::vector<BatchSetType> in_which_p1_as_middle, in_which_p1_as_s, in_which_p1_as_t;
    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> preHop, nextHop;

public:
    VertexIdType *vps;
    BatchSizeType vp_num;
    Graph &g;
    BatchSizeType new_vp_pos_b;

    SignedVertexIdType st_p1v_e;
    const std::vector<std::vector<std::vector<VertexIdType>>> &p1s;
    const std::vector<BatchSizeType> &new2old_posInVPs;
    const std::vector<SignedVertexIdType> &is_p1v;

    OneBatchPath2(VertexIdType *vps, BatchSizeType vp_num, Graph &g, BatchSizeType new_vp_pos_b,
                  SignedVertexIdType p1v_e, SignedVertexIdType st_p1v_e,
                  const std::vector<std::vector<std::vector<VertexIdType>>> &p1s, const std::vector<BatchSizeType> &new2old_posInVPs,
                  const std::vector<SignedVertexIdType> &is_p1v) : in_which_p1_as_middle(p1v_e),
                                                                   in_which_p1_as_s(st_p1v_e),
                                                                   in_which_p1_as_t(st_p1v_e),
                                                                   vps(vps), vp_num(vp_num), g(g), new_vp_pos_b(new_vp_pos_b),
                                                                   st_p1v_e(st_p1v_e), p1s(p1s), new2old_posInVPs(new2old_posInVPs), is_p1v(is_p1v)
    {
    }

    long long record_p1_info_for_batch();

    long long path2(std::vector<std::deque<VertexIdType>> &p2s);
    long long path2_1step2(std::vector<std::deque<VertexIdType>> &p2s);
};

class MultiVP
{
    std::vector<std::vector<std::vector<VertexIdType>>> p1s;
    std::vector<SignedVertexIdType> is_p1v;
    std::vector<BatchSizeType> new2old_posInVPs;

public:
    VertexIdType *vps;
    BatchSizeType total_vp_num;
    Graph &g;
    const std::string &outpf;

    MultiVP(VertexIdType *vps, BatchSizeType total_vp_num, Graph &g, const std::string &outpf) : vps(vps), total_vp_num(total_vp_num), g(g), outpf(outpf)
    {
    }

#ifdef MERGE_PATH_RECORE
    template <typename BatchSetType>
    static void track_path_from_joint(BatchSizeType vp_num, std::vector<std::deque<VertexIdType>>::iterator p1s_it,
                                      typename std::unordered_map<VertexIdType, BatchSetType> &joints,
                                      std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> &prec,
                                      std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> &succ);
#endif

    template <typename BatchSetType>
    static long long path1(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);

    void map_p1(SignedVertexIdType &p1v_e, SignedVertexIdType &st_p1v_e);

    // group_method:
    // 0: uniform
    // 1: degree descending
    // regroup_method:
    // 0: no reorder, just drop vp with no valid p1, then do uniform split
    // path2_method
    // 0: new-paradim, new-graph-represent
    // 1: new-paradim, new-graph-represent, 1step2
    // BatchSetType1: for path1
    // BatchSetType2: for path2
    template <typename BatchSetType1, typename BatchSetType2>
    long long inner_run(int group_method, int regroup_method, int path2_method, int k)
    {
        long long duration = 0;
        auto t1 = std::chrono::steady_clock::now();

        p1s.resize(total_vp_num);
        BatchSizeType ori_total_vp_num = total_vp_num;
        std::cout << "group" << std::endl;
        BatchSizeType batch_num = 0;
        VertexIdType *ranges;
        if (group_method == 0)
            ranges = uniform<BatchSetType1>(g, vps, total_vp_num, batch_num);
        else if (group_method == 1)
            ranges = degreeDescending<BatchSetType1>(g, vps, total_vp_num, batch_num);
        else
        {
            std::cout << "unsupported group_method: " << group_method << std::endl;
            return 0;
        }
        auto t2 = std::chrono::steady_clock::now();
        duration += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
        for (BatchSizeType i = 0; i < batch_num; ++i)
        {
            BatchSizeType b = ranges[2 * i], vp_num = ranges[2 * i + 1] - b;
            std::cout << "p1 for " << b << "-" << b + vp_num << std::endl;
            duration +=
                path1<BatchSetType1>(vps + 2 * b, vp_num, g, p1s.begin() + b, 0);
        }
        t1 = std::chrono::steady_clock::now();
        delete[] ranges;

        for (int kk = 1; kk < k; ++kk)
        {
            std::cout << "path" << kk + 1 << std::endl;
            std::cout << "regroup" << std::endl;
            if (regroup_method == 0)
                ranges = regroup<BatchSetType2>(vps, total_vp_num, batch_num, p1s, new2old_posInVPs, kk);
            else
            {
                std::cout << "unsupported regroup_method: " << regroup_method << std::endl;
                return 0;
            }
            if (total_vp_num == 0)
                break;

            SignedVertexIdType p1v_e, st_p1v_e;
            if (is_p1v.size() < g.get_N())
                is_p1v.resize(g.get_N(), -1);
            else
                is_p1v.assign(g.get_N(), -1);
            map_p1(p1v_e, st_p1v_e);
            t2 = std::chrono::steady_clock::now();
            duration += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();

            for (BatchSizeType i = 0; i < batch_num; ++i)
            {
                BatchSizeType b = ranges[2 * i], vp_num = ranges[2 * i + 1] - b;
                std::cout << "p" << kk + 1 << " for " << b << "-" << b + vp_num << std::endl;
                VertexIdType *local_vps = vps + 2 * b;

                OneBatchPath2<BatchSetType2> ins(local_vps, vp_num, g, b, p1v_e, st_p1v_e, p1s, new2old_posInVPs, is_p1v);
                std::vector<std::deque<VertexIdType>> p2s(vp_num);
                duration +=
                    ins.record_p1_info_for_batch();
                if (path2_method == 0)
                    duration +=
                        ins.path2(p2s);
                else if (path2_method == 1)
                    duration +=
                        ins.path2_1step2(p2s);
                else
                {
                    std::cout << "unsupported path2_method: " << path2_method << std::endl;
                    return 0;
                }

                t1 = std::chrono::steady_clock::now();
                std::cout << "adjust for " << b << "-" << b + vp_num << std::endl;
                for (BatchSizeType j = 0, old; j < vp_num; ++j)
                {
                    old = new2old_posInVPs[b + j];
                    if (int(p1s[old].size()) == kk && p1s[old].back().size() >= 2 && p2s[j].size() >= 2)
                        OneVP::adjust_p1p2(g.get_N(), p1s[old], p2s[j]);
                }
                t2 = std::chrono::steady_clock::now();
                duration += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
            }
            delete[] ranges;
        }

        std::cout << "output" << std::endl;
        for (BatchSizeType j = 0; j < ori_total_vp_num; ++j)
        {
            bool fail = 0;
            if (int(p1s[j].size()) < k)
                fail = 1;
            OneVP::output_paths(fail, vps[j * 2], vps[2 * j + 1], p1s[j], g, outpf);
        }
        delete[] vps;

        return duration;
    }

    template <typename BatchSetType1, typename BatchSetType2>
    long long ms_uniform(int k) { return inner_run<BatchSetType1, BatchSetType2>(0, 0, 0, k); }
    template <typename BatchSetType1, typename BatchSetType2>
    long long ms(int k) { return inner_run<BatchSetType1, BatchSetType2>(1, 0, 0, k); }
    template <typename BatchSetType1, typename BatchSetType2>
    long long ms_concentrate_uniform(int k) { return inner_run<BatchSetType1, BatchSetType2>(0, 0, 1, k); }
    template <typename BatchSetType1, typename BatchSetType2>
    long long ms_concentrate(int k) { return inner_run<BatchSetType1, BatchSetType2>(1, 0, 1, k); }
};

#endif // _MULTIVP_H