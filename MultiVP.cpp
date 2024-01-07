#include "MultiVP.h"

template class OneBatchPath2<BuiltinUIntWrapper<uint8_t>>;
template class OneBatchPath2<BuiltinUIntWrapper<uint16_t>>;
template class OneBatchPath2<BuiltinUIntWrapper<uint32_t>>;
template class OneBatchPath2<BuiltinUIntWrapper<uint64_t>>;
template class OneBatchPath2<M128iWrapper>;

template long long MultiVP::path1<BuiltinUIntWrapper<uint8_t>>(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);
template long long MultiVP::path1<BuiltinUIntWrapper<uint16_t>>(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);
template long long MultiVP::path1<BuiltinUIntWrapper<uint32_t>>(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);
template long long MultiVP::path1<BuiltinUIntWrapper<uint64_t>>(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);
template long long MultiVP::path1<M128iWrapper>(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos);

template <typename BatchSetType>
void MultiVP::track_path_from_joint(BatchSizeType vp_num, std::vector<std::deque<VertexIdType>>::iterator p1s_it,
                                    typename std::unordered_map<VertexIdType, BatchSetType> &joints,
                                    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> &prec,
                                    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> &succ)
{
    struct VertexWithiset
    {
        BatchSetType iset;
        VertexIdType v;
        VertexWithiset() : v((VertexIdType)-1ull) {}
        VertexWithiset(VertexIdType v, const BatchSetType &iset) : iset(iset), v(v) {}
    };
    std::queue<VertexWithiset> q;

    // prec
    VertexIdType v, u;
    for (auto &pr : joints)
    {
        v = pr.first;
        q.emplace(v, pr.second);
        typename BatchSetType::OnesIdx ones(pr.second, vp_num);
        for (auto ins = ones.next(); ~ins; ins = ones.next())
            (p1s_it + ins)->push_front(v);
    }
    while (q.empty() == 0)
    {
        VertexWithiset curr = q.front();
        q.pop();
        v = curr.v;
        for (auto &pr : prec[v])
        {
            BatchSetType iset(curr.iset, pr.second); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                u = pr.first;
                typename BatchSetType::OnesIdx ones(iset, vp_num);
                for (auto ins = ones.next(); ~ins; ins = ones.next())
                    (p1s_it + ins)->push_front(u);
                q.emplace(u, iset);
            }
        }
    }

    // succ
    for (auto &pr : joints)
        q.emplace(pr.first, pr.second);
    while (q.empty() == 0)
    {
        VertexWithiset curr = q.front();
        q.pop();
        v = curr.v;
        for (auto &pr : succ[v])
        {
            BatchSetType iset(curr.iset, pr.second); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                u = pr.first;
                typename BatchSetType::OnesIdx ones(iset, vp_num);
                for (auto ins = ones.next(); ~ins; ins = ones.next())
                    (p1s_it + ins)->push_back(u);
                q.emplace(u, iset);
            }
        }
    }
}

template <typename BatchSetType>
long long MultiVP::path1(VertexIdType *vps, BatchSizeType vp_num, Graph &g, std::vector<std::vector<std::vector<VertexIdType>>>::iterator vp1s_it, int p1_pos)
{
    auto t1 = std::chrono::steady_clock::now();

    std::vector<std::deque<VertexIdType>> dp1s(vp_num);
    auto p1s_it = dp1s.begin();

    VertexIdType *nbrs = g.get_nbrs();
    EdgeNumType *offset = g.get_offset();
    VertexIdType N = g.get_N();

    std::vector<BatchSetType> s_q[2], t_q[2];
    int q_idx = 0, q_next_idx = 1;
    std::vector<BatchSetType> s_seen(N), t_seen(N);
    BatchSetType undone, q_not_empty;
    undone.set_all();
    bool all_done = 0;
    typename std::unordered_map<VertexIdType, BatchSetType> joint;
    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> prec, succ;
    for (int i = 0; i < 2; ++i)
    {
        s_q[i].resize(N);
        t_q[i].resize(N);
    }
    VertexIdType s, t;
    for (BatchSizeType i = 0; i < vp_num; i++)
    {
        s = vps[i * 2];
        t = vps[i * 2 + 1];
        s_q[0][s].set_idx(i);
        s_seen[s].set_idx(i);
        t_q[0][t].set_idx(i);
        t_seen[t].set_idx(i);
    }

    VertexIdType u;
    long long depth = 0, frontier_num = 0;
    while (all_done == 0)
    {
        {
            typename BatchSetType::OnesIdx ones(undone, vp_num);
            std::cout << "depth=" << depth++ << " undone_vp_num=" << ones.size() << " ";
            frontier_num = 0;
        }
        q_not_empty.zero_all();

        // s_q
        for (VertexIdType v = 0; v < N && all_done == 0; ++v)
        {
            BatchSetType iset(s_q[q_idx][v], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                ++frontier_num;
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, s_seen[u], 0); // a&(~b)
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        s_seen[u] |= new_iset;
                        // prec[u] = v;
                        prec[u][v] |= new_iset;
                        BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
                            // joint = u;
                            joint[u] |= cross;
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }
                        s_q[q_next_idx][u] |= new_iset;
                        q_not_empty |= new_iset;
                    }
                }
                s_q[q_idx][v].zero_all();
            }
        }

        std::cout << "s_frontier_num=" << frontier_num << " ";
        frontier_num = 0;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }
        q_not_empty.zero_all();

        // t_q
        for (VertexIdType v = 0; v < N && all_done == 0; ++v)
        {
            BatchSetType iset(t_q[q_idx][v], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                ++frontier_num;
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, t_seen[u], 0); // a&(~b)
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        t_seen[u] |= new_iset;
                        // succ[u] = v;
                        succ[u][v] |= new_iset;
                        BatchSetType cross(t_seen[u], s_seen[u]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
                            // joint = u;
                            joint[u] |= cross;
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }
                        t_q[q_next_idx][u] |= new_iset;
                        q_not_empty |= new_iset;
                    }
                }
                t_q[q_idx][v].zero_all();
            }
        }

        std::cout << "t_frontier_num=" << frontier_num << std::endl;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }

        q_idx = q_next_idx;
        q_next_idx = 1 - q_next_idx;
    }

    std::cout << std::endl;

    MultiVP::track_path_from_joint(vp_num, p1s_it, joint, prec, succ);
    for (BatchSizeType i = 0; i < vp_num; i++)
    {
        (vp1s_it + i)->resize(p1_pos + 1);
        std::copy(dp1s[i].begin(), dp1s[i].end(), std::back_inserter((*(vp1s_it + i))[p1_pos]));
    }

    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

void MultiVP::map_p1(SignedVertexIdType &p1v_e, SignedVertexIdType &st_p1v_e)
{
    SignedVertexIdType p1v = 0;
    // map s&t
    VertexIdType v;
    for (BatchSizeType i = 0; i < total_vp_num; ++i)
    {
        v = vps[2 * i];
        if (is_p1v[v] < 0)
            is_p1v[v] = p1v++;
        v = vps[2 * i + 1];
        if (is_p1v[v] < 0)
            is_p1v[v] = p1v++;
    }
    st_p1v_e = p1v;
    // map other p1v
    for (BatchSizeType i = 0, old; i < total_vp_num; ++i)
    {
        old = new2old_posInVPs[i];
        for (size_t k = 0, ed = p1s[old].size(); k < ed; ++k)
        {
            size_t sz = p1s[old][k].size();
            if (sz >= 2)
            {
                --sz;
                for (size_t j = 1; j < sz; ++j)
                {
                    v = p1s[old][k][j];
                    if (is_p1v[v] < 0)
                        is_p1v[v] = p1v++;
                }
            }
        }
    }
    p1v_e = p1v;
}

template <typename BatchSetType>
long long OneBatchPath2<BatchSetType>::record_p1_info_for_batch()
{
    auto t1 = std::chrono::steady_clock::now();

    for (BatchSizeType new_pos = new_vp_pos_b, i = 0; i < vp_num; ++new_pos, ++i)
    {
        BatchSetType iset;
        iset.set_idx(i);
        BatchSizeType old_pos = new2old_posInVPs[new_pos];
        VertexIdType s = vps[2 * i], t = vps[2 * i + 1];
        in_which_p1_as_s[is_p1v[s]] |= iset;
        in_which_p1_as_t[is_p1v[t]] |= iset;
        for (size_t k = 0, ked = p1s[old_pos].size(); k < ked; ++k)
        {
            VertexIdType pre = s, curr;
            for (size_t j = 1, sz = p1s[old_pos][k].size(); j < sz; ++j)
            {
                if (j > 1)
                    in_which_p1_as_middle[is_p1v[pre]] |= iset;
                curr = p1s[old_pos][k][j];
                nextHop[pre][curr] |= iset;
                preHop[curr][pre] |= iset;
                pre = curr;
            }
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

template <typename BatchSetType>
long long OneBatchPath2<BatchSetType>::path2(std::vector<std::deque<VertexIdType>> &p2s)
{
    auto t1 = std::chrono::steady_clock::now();
    VertexIdType *nbrs = g.get_nbrs();
    EdgeNumType *offset = g.get_offset();
    VertexIdType N = g.get_N(), new_N = N * 2;

    std::vector<BatchSetType> s_q[2], t_q[2];
    int q_idx = 0, q_next_idx = 1;
    std::vector<BatchSetType> s_seen(new_N), t_seen(new_N);
    BatchSetType undone, q_not_empty;
    undone.set_all();
    bool all_done = 0;
    typename std::unordered_map<VertexIdType, BatchSetType> joint;
    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> prec, succ;
    for (int i = 0; i < 2; ++i)
    {
        s_q[i].resize(new_N);
        t_q[i].resize(new_N);
    }
    VertexIdType s, t;
    for (BatchSizeType i = 0; i < vp_num; i++)
    {
        s = vps[i * 2];
        t = vps[i * 2 + 1];
        s_q[0][s].set_idx(i);
        s_seen[s].set_idx(i);
        t_q[0][t].set_idx(i);
        t_seen[t].set_idx(i);
    }

    VertexIdType u, n, v;
    SignedVertexIdType p1v, p1u;
    typename std::unordered_map<VertexIdType, BatchSetType>::iterator it;
    long long depth = 0, frontier_num = 0;
    while (all_done == 0)
    {
        {
            typename BatchSetType::OnesIdx ones(undone, vp_num);
            std::cout << "depth=" << depth++ << " undone_vp_num=" << ones.size() << " ";
            frontier_num = 0;
        }
        q_not_empty.zero_all();

        // s_q: out_nbr
        for (VertexIdType v0 = 0; v0 < new_N && all_done == 0; ++v0)
        {
            BatchSetType iset(s_q[q_idx][v0], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                frontier_num++;
                if (v0 >= N)
                    v = v0 - N;
                else
                    v = v0;

                // 1. v is not in p1, or v is t // out_nbr: nbr_in_g
                // 2. v is s // out_nbr: nbr_in_g-next_hop
                // 3. v is p1-middle-vertex
                // 3.1. if v is v.OUT (v0!=v) // out_nbr: nbr_in_g-next_hop, v.IN
                // 3.2. if v is v.IN (v0==v) // out_nbr: pre_hop.OUT

                p1v = is_p1v[v];
                // 3.2. v is v.IN
                if (p1v >= 0 && v == v0)
                {
                    BatchSetType is_vin(iset, in_which_p1_as_middle[p1v]); // a&b
                    if (is_vin.is_all_zero(vp_num) == 0)
                    {
                        for (auto &pr : preHop.at(v))
                        {
                            n = pr.first;
                            BatchSetType vhasn(is_vin, s_seen[n + N], pr.second); // (a&(~b))&c
                            if (vhasn.is_all_zero(vp_num) == 0)
                            {
                                s_seen[n + N] |= vhasn;
                                // prec[n + N] = v;
                                prec[n + N][v] |= vhasn;

                                BatchSetType cross(s_seen[n + N], t_seen[n + N]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[n + N] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }
                                s_q[q_next_idx][n + N] |= vhasn;
                                q_not_empty |= vhasn;
                            }
                        }

                        iset.andnot(is_vin);
                    }
                }

                bool some_skip = 0;
                typename std::unordered_map<VertexIdType, BatchSetType> *next_hops = 0;
                typename std::unordered_map<VertexIdType, BatchSetType>::iterator ed;
                BatchSetType skip_next_hop;
                if (p1v >= 0)
                {
                    if (p1v < st_p1v_e)
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v] | in_which_p1_as_s[p1v]);
                    else
                        skip_next_hop = iset & in_which_p1_as_middle[p1v];
                    if (skip_next_hop.is_all_zero(vp_num) == 0)
                    {
                        next_hops = &(nextHop.at(v));
                        ed = next_hops->end();
                        some_skip = 1;
                    }
                }
                // out_nbr: nbr_in_g(-next_hop)
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, s_seen[u], 0); // a&(~b)
                    p1u = is_p1v[u];
                    // for these vps, they need to skip u: u is v's nexthop in some p1
                    if (some_skip && p1u >= 0 && ((it = next_hops->find(u)) != ed))
                        new_iset.andnot(skip_next_hop & it->second);
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        if (p1u >= 0)
                        {
                            // u is not in p1 || u is p1-terminal (i.e. u is not p1-middle-vertex): normal enqueue
                            BatchSetType normal(new_iset, in_which_p1_as_middle[p1u], 0); // a&(~b)
                            if (normal.is_all_zero(vp_num) == 0)
                            {
                                s_seen[u] |= normal;

                                if (p1v >= 0)
                                {
                                    // if (v_is_middle_p1)
                                    //     prec[u] = v + N;
                                    // else
                                    //     prec[u] = v;
                                    BatchSetType v_is_middle_p1(normal, in_which_p1_as_middle[p1v]);        // a&b
                                    BatchSetType v_is_not_middle_p1(normal, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            prec[ins][u] = v + N;
                                    }
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                        for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                            prec[ins][u] = v;
                                    }
#else
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v + N] |= v_is_middle_p1;
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v] |= v_is_not_middle_p1;
#endif
                                }
                                else
                                {
                                    // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(normal, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v;
#else
                                    prec[u][v] |= normal;
#endif
                                }

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                s_q[q_next_idx][u] |= normal;
                                q_not_empty |= normal;
                            }
                            // u is p1-middle-vertex: then u is u.IN
                            BatchSetType umiddle(new_iset, in_which_p1_as_middle[p1u]); // a&b
                            if (umiddle.is_all_zero(vp_num) == 0)
                            {
                                s_seen[u] |= umiddle;

                                if (p1v >= 0)
                                {
                                    // if (v_is_middle_p1)
                                    //     prec[u] = v + N;
                                    // else
                                    //     prec[u] = v;
                                    BatchSetType v_is_middle_p1(umiddle, in_which_p1_as_middle[p1v]);        // a&b
                                    BatchSetType v_is_not_middle_p1(umiddle, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            prec[ins][u] = v + N;
                                    }
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                        for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                            prec[ins][u] = v;
                                    }
#else
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v + N] |= v_is_middle_p1;
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v] |= v_is_not_middle_p1;
#endif
                                }
                                else
                                {
                                    // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(umiddle, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v;
#else
                                    prec[u][v] |= umiddle;
#endif
                                }

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }
                                s_q[q_next_idx][u] |= umiddle;
                                q_not_empty |= umiddle;
                            }
                        }
                        // u is not in p1: normal enqueue
                        else
                        {
                            s_seen[u] |= new_iset;

                            if (p1v >= 0)
                            {
                                // if (v_is_middle_p1)
                                //     prec[u] = v + N;
                                // else
                                //     prec[u] = v;
                                BatchSetType v_is_middle_p1(new_iset, in_which_p1_as_middle[p1v]);        // a&b
                                BatchSetType v_is_not_middle_p1(new_iset, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                {
                                    typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v + N;
                                }
                                if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                {
                                    typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                    for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                        prec[ins][u] = v;
                                }
#else
                                if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    prec[u][v + N] |= v_is_middle_p1;
                                if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    prec[u][v] |= v_is_not_middle_p1;
#endif
                            }
                            else
                            {
                                // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(new_iset, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    prec[ins][u] = v;
#else
                                prec[u][v] |= new_iset;
#endif
                            }

                            BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                            cross &= undone;
                            if (cross.is_all_zero(vp_num) == 0)
                            {
                                undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(cross, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    joint[ins] = u;
#else
                                joint[u] |= cross;
#endif
                                if (undone.is_all_zero(vp_num))
                                {
                                    all_done = 1;
                                    break;
                                }
                            }

                            s_q[q_next_idx][u] |= new_iset;
                            q_not_empty |= new_iset;
                        }
                    }
                }

                // out_nbr: v.IN
                if (p1v >= 0)
                {
                    BatchSetType has_vin(iset, in_which_p1_as_middle[p1v]); // a&b
                    has_vin &= undone;
                    has_vin.andnot(s_seen[v]);
                    if (has_vin.is_all_zero(vp_num) == 0)
                    {
                        s_seen[v] |= has_vin;
// prec[v] = v + N;
#ifndef MERGE_PATH_RECORE
                        typename BatchSetType::OnesIdx ones(has_vin, vp_num);
                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                            prec[ins][v] = v + N;
#else
                        prec[v][v + N] |= has_vin;
#endif

                        BatchSetType cross(s_seen[v], t_seen[v]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                joint[ins] = v;
#else
                            joint[v] |= cross;
#endif
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }
                        s_q[q_next_idx][v] |= has_vin;
                        q_not_empty |= has_vin;
                    }
                }
                s_q[q_idx][v0].zero_all();
            }
        }

        std::cout << "s_frontier_num=" << frontier_num << " ";
        frontier_num = 0;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }
        q_not_empty.zero_all();

        // t_q: in_nbr
        for (VertexIdType v0 = 0; v0 < new_N && all_done == 0; ++v0)
        {
            BatchSetType iset(t_q[q_idx][v0], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                ++frontier_num;

                if (v0 >= N)
                    v = v0 - N;
                else
                    v = v0;

                // 1. v is not in p1, or v is s // in_nbr: nbr_in_g
                // 2. v is t // in_nbr: nbr_in_g-pre_hop
                // 3. v is p1-middle-vertex
                // 3.1. if v is v.IN (v0==v) // in_nbr: nbr_in_g-pre_hop, v.OUT
                // 3.2. if v is v.OUT (v0!=v) // in_nbr: next_hop.IN

                p1v = is_p1v[v];
                // 3.2. v is v.OUT
                if (p1v >= 0 && v0 != v)
                {
                    BatchSetType is_vout(iset, in_which_p1_as_middle[p1v]); // a&b
                    if (is_vout.is_all_zero(vp_num) == 0)
                    {
                        for (auto &pr : nextHop.at(v))
                        {
                            n = pr.first;
                            BatchSetType vhasn(is_vout, t_seen[n], pr.second); // (a&(~b))&c
                            if (vhasn.is_all_zero(vp_num) == 0)
                            {
                                t_seen[n] |= vhasn;
// succ[n] = v + N;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(vhasn, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    succ[ins][n] = v + N;
#else
                                succ[n][v + N] |= vhasn;
#endif

                                BatchSetType cross(s_seen[n], t_seen[n]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
// joint = n;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(cross, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        joint[ins] = n;
#else
                                    joint[n] |= cross;
#endif
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }
                                t_q[q_next_idx][n] |= vhasn;
                                q_not_empty |= vhasn;
                            }
                        }

                        iset.andnot(is_vout);
                    }
                }

                bool some_skip = 0;
                typename std::unordered_map<VertexIdType, BatchSetType> *next_hops = 0;
                typename std::unordered_map<VertexIdType, BatchSetType>::iterator ed;
                BatchSetType skip_next_hop;
                if (p1v >= 0)
                {
                    if (p1v < st_p1v_e)
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v] | in_which_p1_as_t[p1v]);
                    else
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v]);
                    if (skip_next_hop.is_all_zero(vp_num) == 0)
                    {
                        next_hops = &(preHop.at(v));
                        ed = next_hops->end();
                        some_skip = 1;
                    }
                }
                // in_nbr: nbr_in_g(-pre_hop)
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, t_seen[u], 0); // a&(~b)
                    p1u = is_p1v[u];
                    // for these vps, they need to skip u: u is v's prehop in some p1
                    if (some_skip && p1u >= 0 && ((it = next_hops->find(u)) != ed))
                        new_iset.andnot(skip_next_hop & it->second);
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        // u is in some's p1
                        if (p1u >= 0)
                        {
                            // u is not in p1 || u is p1-terminal (i.e. u is not p1-middle-vertex): normal enqueue
                            BatchSetType normal(new_iset, in_which_p1_as_middle[p1u], 0); // a&(~b)
                            if (normal.is_all_zero(vp_num) == 0)
                            {
                                t_seen[u] |= normal;
// succ[u] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones2(normal, vp_num);
                                for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                    succ[ins][u] = v;
#else
                                succ[u][v] |= normal;
#endif

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                t_q[q_next_idx][u] |= normal;
                                q_not_empty |= normal;
                            }
                            // u is p1-middle-vertex: then u is u.OUT
                            BatchSetType umiddle(new_iset, t_seen[u + N], in_which_p1_as_middle[p1u]); // (a&(~b))&c
                            if (umiddle.is_all_zero(vp_num) == 0)
                            {
                                t_seen[u + N] |= umiddle;

// succ[u + N] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones2(umiddle, vp_num);
                                for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                    succ[ins][u + N] = v;
#else
                                succ[u + N][v] |= umiddle;
#endif

                                BatchSetType cross(s_seen[u + N], t_seen[u + N]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(cross, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        joint[ins] = u + N;
#else
                                    joint[u + N] |= cross;
#endif
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                t_q[q_next_idx][u + N] |= umiddle;
                                q_not_empty |= umiddle;
                            }
                        }
                        // u is not in p1: normal enqueue
                        else
                        {
                            t_seen[u] |= new_iset;
// succ[u] = v;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones2(new_iset, vp_num);
                            for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                succ[ins][u] = v;
#else
                            succ[u][v] |= new_iset;
#endif

                            BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                            cross &= undone;
                            if (cross.is_all_zero(vp_num) == 0)
                            {
                                undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(cross, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    joint[ins] = u;
#else
                                joint[u] |= cross;
#endif
                                if (undone.is_all_zero(vp_num))
                                {
                                    all_done = 1;
                                    break;
                                }
                            }

                            t_q[q_next_idx][u] |= new_iset;
                            q_not_empty |= new_iset;
                        }
                    }
                }

                // in_nbr: v.OUT
                if (p1v >= 0)
                {
                    BatchSetType has_vin(iset, in_which_p1_as_middle[p1v]); // a&b
                    has_vin &= undone;
                    has_vin.andnot(t_seen[v + N]);
                    if (has_vin.is_all_zero(vp_num) == 0)
                    {
                        t_seen[v + N] |= has_vin;
// succ[v + N] = v;
#ifndef MERGE_PATH_RECORE
                        typename BatchSetType::OnesIdx ones(has_vin, vp_num);
                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                            succ[ins][v + N] = v;
#else
                        succ[v + N][v] |= has_vin;
#endif

                        BatchSetType cross(s_seen[v + N], t_seen[v + N]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                joint[ins] = v + N;
#else
                            joint[v + N] |= cross;
#endif
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }

                        t_q[q_next_idx][v + N] |= has_vin;
                        q_not_empty |= has_vin;
                    }
                }
                t_q[q_idx][v0].zero_all();
            }
        }

        std::cout << "t_frontier_num=" << frontier_num << std::endl;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }

        q_idx = q_next_idx;
        q_next_idx = 1 - q_next_idx;
    }

    std::cout << std::endl;

#ifndef MERGE_PATH_RECORE
    for (BatchSizeType i = 0; i < vp_num; ++i)
    {
        if (~joint[i])
        {
            OneVP::track_path_from_joint(p2s[i], joint[i], prec[i], succ[i]);
        }
    }
#else
    MultiVP::track_path_from_joint(vp_num, p2s.begin(), joint, prec, succ);
#endif
    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

template <typename BatchSetType>
long long OneBatchPath2<BatchSetType>::path2_1step2(std::vector<std::deque<VertexIdType>> &p2s)
{
    auto t1 = std::chrono::steady_clock::now();
    VertexIdType *nbrs = g.get_nbrs();
    EdgeNumType *offset = g.get_offset();
    VertexIdType N = g.get_N();

    std::vector<BatchSetType> s_q[2], t_q[2];
    int q_idx = 0, q_next_idx = 1;
    std::vector<BatchSetType> s_seen(N * 2), t_seen(N * 2);
    BatchSetType undone, q_not_empty;
    undone.set_all();
    bool all_done = 0;
#ifndef MERGE_PATH_RECORE
    std::unordered_map<VertexIdType, VertexIdType> prec[vp_num], succ[vp_num];
    VertexIdType joint[vp_num];
    memset(joint, -1, vp_num * sizeof(VertexIdType));
#else
    typename std::unordered_map<VertexIdType, BatchSetType> joint;
    std::unordered_map<VertexIdType, typename std::unordered_map<VertexIdType, BatchSetType>> prec, succ;
#endif
    for (int i = 0; i < 2; ++i)
    {
        s_q[i].resize(N);
        t_q[i].resize(N);
    }
    VertexIdType s, t;
    for (BatchSizeType i = 0; i < vp_num; i++)
    {
        s = vps[i * 2];
        t = vps[i * 2 + 1];
        s_q[0][s].set_idx(i);
        s_seen[s].set_idx(i);
        t_q[0][t].set_idx(i);
        t_seen[t].set_idx(i);
    }

    VertexIdType u, n;
    SignedVertexIdType p1v, p1u;
    typename std::unordered_map<VertexIdType, BatchSetType>::iterator it;
    long long depth = 0, frontier_num = 0;
    while (all_done == 0)
    {
        {
            typename BatchSetType::OnesIdx ones(undone, vp_num);
            std::cout << "depth=" << depth++ << " undone_vp_num=" << ones.size() << " ";
            frontier_num = 0;
        }
        q_not_empty.zero_all();

        // s_q: out_nbr
        for (VertexIdType v = 0; v < N && all_done == 0; ++v)
        {
            BatchSetType iset(s_q[q_idx][v], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                frontier_num++;
                // 1. v is not in p1, or v is t // out_nbr: nbr_in_g
                // 2. v is s // out_nbr: nbr_in_g-next_hop
                // 3. v is p1-middle-vertex
                // v is v.OUT // out_nbr: nbr_in_g-next_hop, v.IN

                p1v = is_p1v[v];
                bool some_skip = 0;
                typename std::unordered_map<VertexIdType, BatchSetType> *next_hops = 0;
                typename std::unordered_map<VertexIdType, BatchSetType>::iterator ed;
                BatchSetType skip_next_hop;
                if (p1v >= 0)
                {
                    if (p1v < st_p1v_e)
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v] | in_which_p1_as_s[p1v]);
                    else
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v]);
                    if (skip_next_hop.is_all_zero(vp_num) == 0)
                    {
                        next_hops = &(nextHop.at(v));
                        ed = next_hops->end();
                        some_skip = 1;
                    }
                }
                // out_nbr: nbr_in_g(-next_hop)
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, s_seen[u], 0); // a&(~b)
                    p1u = is_p1v[u];
                    // for these vps, they need to skip u: u is v's nexthop in some p1
                    if (some_skip && p1u >= 0 && ((it = next_hops->find(u)) != ed))
                        new_iset.andnot(skip_next_hop & it->second);
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        if (p1u >= 0)
                        {
                            // u is not in p1 || u is p1-terminal (i.e. u is not p1-middle-vertex): normal enqueue
                            BatchSetType normal(new_iset, in_which_p1_as_middle[p1u], 0); // a&(~b)
                            if (normal.is_all_zero(vp_num) == 0)
                            {
                                s_seen[u] |= normal;

                                if (p1v >= 0)
                                {
                                    // if (v_is_middle_p1)
                                    //     prec[u] = v + N;
                                    // else
                                    //     prec[u] = v;
                                    BatchSetType v_is_middle_p1(normal, in_which_p1_as_middle[p1v]);        // a&b
                                    BatchSetType v_is_not_middle_p1(normal, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            prec[ins][u] = v + N;
                                    }
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                        for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                            prec[ins][u] = v;
                                    }
#else
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v + N] |= v_is_middle_p1;
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v] |= v_is_not_middle_p1;
#endif
                                }
                                else
                                {
                                    // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(normal, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v;
#else
                                    prec[u][v] |= normal;
#endif
                                }

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                s_q[q_next_idx][u] |= normal;
                                q_not_empty |= normal;
                            }
                            // u is p1-middle-vertex: then u is u.IN
                            BatchSetType umiddle(new_iset, in_which_p1_as_middle[p1u]); // a&b
                            if (umiddle.is_all_zero(vp_num) == 0)
                            {
                                s_seen[u] |= umiddle;

                                if (p1v >= 0)
                                {
                                    // if (v_is_middle_p1)
                                    //     prec[u] = v + N;
                                    // else
                                    //     prec[u] = v;
                                    BatchSetType v_is_middle_p1(umiddle, in_which_p1_as_middle[p1v]);        // a&b
                                    BatchSetType v_is_not_middle_p1(umiddle, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            prec[ins][u] = v + N;
                                    }
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    {
                                        typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                        for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                            prec[ins][u] = v;
                                    }
#else
                                    if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v + N] |= v_is_middle_p1;
                                    if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                        prec[u][v] |= v_is_not_middle_p1;
#endif
                                }
                                else
                                {
                                    // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(umiddle, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v;
#else
                                    prec[u][v] |= umiddle;
#endif
                                }

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                umiddle &= undone;
                                for (auto &pr : preHop.at(u))
                                {
                                    n = pr.first;
                                    BatchSetType uhasn(umiddle, s_seen[n + N], pr.second); // (a&(~b))&c
                                    if (uhasn.is_all_zero(vp_num) == 0)
                                    {
                                        s_seen[n + N] |= uhasn;
// prec[n + N] = u;
#ifndef MERGE_PATH_RECORE
                                        typename BatchSetType::OnesIdx ones(uhasn, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            prec[ins][n + N] = u;
#else
                                        prec[n + N][u] |= uhasn;
#endif

                                        BatchSetType cross(s_seen[n + N], t_seen[n + N]); // a&b
                                        cross &= undone;
                                        if (cross.is_all_zero(vp_num) == 0)
                                        {
                                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                                joint[ins] = n + N;
#else
                                            joint[n + N] |= cross;
#endif
                                            if (undone.is_all_zero(vp_num))
                                            {
                                                all_done = 1;
                                                break;
                                            }
                                        }
                                    }

                                    s_q[q_next_idx][n] |= uhasn;
                                    q_not_empty |= uhasn;
                                }
                            }
                        }
                        // u is not in p1: normal enqueue
                        else
                        {
                            s_seen[u] |= new_iset;

                            if (p1v >= 0)
                            {
                                // if (v_is_middle_p1)
                                //     prec[u] = v + N;
                                // else
                                //     prec[u] = v;
                                BatchSetType v_is_middle_p1(new_iset, in_which_p1_as_middle[p1v]);        // a&b
                                BatchSetType v_is_not_middle_p1(new_iset, in_which_p1_as_middle[p1v], 0); // a&(~b)
#ifndef MERGE_PATH_RECORE
                                if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                {
                                    typename BatchSetType::OnesIdx ones(v_is_middle_p1, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        prec[ins][u] = v + N;
                                }
                                if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                {
                                    typename BatchSetType::OnesIdx ones2(v_is_not_middle_p1, vp_num);
                                    for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                        prec[ins][u] = v;
                                }
#else
                                if (v_is_middle_p1.is_all_zero(vp_num) == 0)
                                    prec[u][v + N] |= v_is_middle_p1;
                                if (v_is_not_middle_p1.is_all_zero(vp_num) == 0)
                                    prec[u][v] |= v_is_not_middle_p1;
#endif
                            }
                            else
                            {
                                // prec[u] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(new_iset, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    prec[ins][u] = v;
#else
                                prec[u][v] |= new_iset;
#endif
                            }

                            BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                            cross &= undone;
                            if (cross.is_all_zero(vp_num) == 0)
                            {
                                undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(cross, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    joint[ins] = u;
#else
                                joint[u] |= cross;
#endif
                                if (undone.is_all_zero(vp_num))
                                {
                                    all_done = 1;
                                    break;
                                }
                            }

                            s_q[q_next_idx][u] |= new_iset;
                            q_not_empty |= new_iset;
                        }
                    }
                }

                // out_nbr: v.IN
                if (p1v >= 0)
                {
                    BatchSetType has_vin(iset, in_which_p1_as_middle[p1v]); // a&b
                    has_vin &= undone;
                    has_vin.andnot(s_seen[v]);
                    if (has_vin.is_all_zero(vp_num) == 0)
                    {
                        s_seen[v] |= has_vin;
// prec[v] = v + N;
#ifndef MERGE_PATH_RECORE
                        typename BatchSetType::OnesIdx ones(has_vin, vp_num);
                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                            prec[ins][v] = v + N;
#else
                        prec[v][v + N] |= has_vin;
#endif

                        BatchSetType cross(s_seen[v], t_seen[v]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                joint[ins] = v;
#else
                            joint[v] |= cross;
#endif
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }

                        has_vin &= undone;
                        for (auto &pr : preHop.at(v))
                        {
                            n = pr.first;
                            BatchSetType vhasn(has_vin, s_seen[n + N], pr.second); // (a&(~b))&c
                            if (vhasn.is_all_zero(vp_num) == 0)
                            {
                                s_seen[n + N] |= vhasn;
// prec[n + N] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(vhasn, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    prec[ins][n + N] = v;
#else
                                prec[n + N][v] |= vhasn;
#endif

                                BatchSetType cross(s_seen[n + N], t_seen[n + N]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(cross, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        joint[ins] = n + N;
#else
                                    joint[n + N] |= cross;
#endif
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }
                                s_q[q_next_idx][n] |= vhasn;
                                q_not_empty |= vhasn;
                            }
                        }
                    }
                }
                s_q[q_idx][v].zero_all();
            }
        }

        std::cout << "s_frontier_num=" << frontier_num << " ";
        frontier_num = 0;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }
        q_not_empty.zero_all();

        // t_q: in_nbr
        for (VertexIdType v = 0; v < N && all_done == 0; ++v)
        {
            BatchSetType iset(t_q[q_idx][v], undone); // a&b
            if (iset.is_all_zero(vp_num) == 0)
            {
                ++frontier_num;
                // 1. v is not in p1, or v is s // in_nbr: nbr_in_g
                // 2. v is t // in_nbr: nbr_in_g-pre_hop
                // 3. v is p1-middle-vertex
                // v is v.IN // in_nbr: nbr_in_g-pre_hop, v.OUT

                p1v = is_p1v[v];
                bool some_skip = 0;
                typename std::unordered_map<VertexIdType, BatchSetType> *next_hops = 0;
                typename std::unordered_map<VertexIdType, BatchSetType>::iterator ed;
                BatchSetType skip_next_hop;
                if (p1v >= 0)
                {
                    if (p1v < st_p1v_e)
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v] | in_which_p1_as_t[p1v]);
                    else
                        skip_next_hop = iset & (in_which_p1_as_middle[p1v]);
                    if (skip_next_hop.is_all_zero(vp_num) == 0)
                    {
                        next_hops = &(preHop.at(v));
                        ed = next_hops->end();
                        some_skip = 1;
                    }
                }
                // in_nbr: nbr_in_g(-pre_hop)
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    BatchSetType new_iset(iset, t_seen[u], 0); // a&(~b)
                    p1u = is_p1v[u];
                    // for these vps, they need to skip u: u is v's prehop in some p1
                    if (some_skip && p1u >= 0 && ((it = next_hops->find(u)) != ed))
                        new_iset.andnot(skip_next_hop & it->second);
                    if (new_iset.is_all_zero(vp_num) == 0)
                    {
                        if (p1u >= 0)
                        {
                            // u is not in p1 || u is p1-terminal (i.e. u is not p1-middle-vertex): normal enqueue
                            BatchSetType normal(new_iset, in_which_p1_as_middle[p1u], 0); // a&(~b)
                            if (normal.is_all_zero(vp_num) == 0)
                            {
                                t_seen[u] |= normal;
// succ[u] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones2(normal, vp_num);
                                for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                    succ[ins][u] = v;
#else
                                succ[u][v] |= normal;
#endif

                                BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
                                    // joint = u;
                                    joint[u] |= cross;
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                t_q[q_next_idx][u] |= normal;
                                q_not_empty |= normal;
                            }
                            // u is p1-middle-vertex: then u is u.OUT
                            BatchSetType umiddle(new_iset, t_seen[u + N], in_which_p1_as_middle[p1u]); // (a&(~b))&c
                            if (umiddle.is_all_zero(vp_num) == 0)
                            {
                                t_seen[u + N] |= umiddle;

// succ[u + N] = v;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones2(umiddle, vp_num);
                                for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                    succ[ins][u + N] = v;
#else
                                succ[u + N][v] |= umiddle;
#endif

                                BatchSetType cross(s_seen[u + N], t_seen[u + N]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(cross, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        joint[ins] = u + N;
#else
                                    joint[u + N] |= cross;
#endif
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }

                                umiddle &= undone;
                                for (auto &pr : nextHop.at(u))
                                {
                                    n = pr.first;
                                    BatchSetType uhasn(umiddle, t_seen[n], pr.second); // (a&(~b))&c
                                    if (uhasn.is_all_zero(vp_num) == 0)
                                    {
                                        t_seen[n] |= uhasn;
// succ[n] = u + N;
#ifndef MERGE_PATH_RECORE
                                        typename BatchSetType::OnesIdx ones(uhasn, vp_num);
                                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                                            succ[ins][n] = u + N;
#else
                                        succ[n][u + N] |= uhasn;
#endif

                                        BatchSetType cross(s_seen[n], t_seen[n]); // a&b
                                        cross &= undone;
                                        if (cross.is_all_zero(vp_num) == 0)
                                        {
                                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                                joint[ins] = n;
#else
                                            joint[n] |= cross;
#endif
                                            if (undone.is_all_zero(vp_num))
                                            {
                                                all_done = 1;
                                                break;
                                            }
                                        }
                                    }

                                    t_q[q_next_idx][n] |= uhasn;
                                    q_not_empty |= uhasn;
                                }
                            }
                        }
                        // u is not in p1: normal enqueue
                        else
                        {
                            t_seen[u] |= new_iset;
// succ[u] = v;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones2(new_iset, vp_num);
                            for (auto ins = ones2.next(); ~ins; ins = ones2.next())
                                succ[ins][u] = v;
#else
                            succ[u][v] |= new_iset;
#endif

                            BatchSetType cross(s_seen[u], t_seen[u]); // a&b
                            cross &= undone;
                            if (cross.is_all_zero(vp_num) == 0)
                            {
                                undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(cross, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    joint[ins] = u;
#else
                                joint[u] |= cross;
#endif
                                if (undone.is_all_zero(vp_num))
                                {
                                    all_done = 1;
                                    break;
                                }
                            }

                            t_q[q_next_idx][u] |= new_iset;
                            q_not_empty |= new_iset;
                        }
                    }
                }

                // in_nbr: v.OUT
                if (p1v >= 0)
                {
                    BatchSetType has_vin(iset, in_which_p1_as_middle[p1v]); // a&b
                    has_vin &= undone;
                    has_vin.andnot(t_seen[v + N]);
                    if (has_vin.is_all_zero(vp_num) == 0)
                    {
                        t_seen[v + N] |= has_vin;
// succ[v + N] = v;
#ifndef MERGE_PATH_RECORE
                        typename BatchSetType::OnesIdx ones(has_vin, vp_num);
                        for (auto ins = ones.next(); ~ins; ins = ones.next())
                            succ[ins][v + N] = v;
#else
                        succ[v + N][v] |= has_vin;
#endif

                        BatchSetType cross(s_seen[v + N], t_seen[v + N]); // a&b
                        cross &= undone;
                        if (cross.is_all_zero(vp_num) == 0)
                        {
                            undone.andnot(cross);
// joint = u;
#ifndef MERGE_PATH_RECORE
                            typename BatchSetType::OnesIdx ones(cross, vp_num);
                            for (auto ins = ones.next(); ~ins; ins = ones.next())
                                joint[ins] = v + N;
#else
                            joint[v + N] |= cross;
#endif
                            if (undone.is_all_zero(vp_num))
                            {
                                all_done = 1;
                                break;
                            }
                        }

                        has_vin &= undone;
                        for (auto &pr : nextHop.at(v))
                        {
                            n = pr.first;
                            BatchSetType vhasn(has_vin, t_seen[n], pr.second); // (a&(~b))&c
                            if (vhasn.is_all_zero(vp_num) == 0)
                            {
                                t_seen[n] |= vhasn;
// succ[n] = v + N;
#ifndef MERGE_PATH_RECORE
                                typename BatchSetType::OnesIdx ones(vhasn, vp_num);
                                for (auto ins = ones.next(); ~ins; ins = ones.next())
                                    succ[ins][n] = v + N;
#else
                                succ[n][v + N] |= vhasn;
#endif

                                BatchSetType cross(s_seen[n], t_seen[n]); // a&b
                                cross &= undone;
                                if (cross.is_all_zero(vp_num) == 0)
                                {
                                    undone.andnot(cross);
// joint = n;
#ifndef MERGE_PATH_RECORE
                                    typename BatchSetType::OnesIdx ones(cross, vp_num);
                                    for (auto ins = ones.next(); ~ins; ins = ones.next())
                                        joint[ins] = n;
#else
                                    joint[n] |= cross;
#endif
                                    if (undone.is_all_zero(vp_num))
                                    {
                                        all_done = 1;
                                        break;
                                    }
                                }
                                t_q[q_next_idx][n] |= vhasn;
                                q_not_empty |= vhasn;
                            }
                        }
                    }
                }
                t_q[q_idx][v].zero_all();
            }
        }

        std::cout << "t_frontier_num=" << frontier_num << std::endl;

        if (all_done)
            break;
        undone &= q_not_empty; // only undone vp whose next_queue is not empty needs to explore
        if (undone.is_all_zero(vp_num))
        {
            all_done = 1;
            break;
        }

        q_idx = q_next_idx;
        q_next_idx = 1 - q_next_idx;
    }

    std::cout << std::endl;

#ifndef MERGE_PATH_RECORE
    for (BatchSizeType i = 0; i < vp_num; ++i)
    {
        if (~joint[i])
        {
            OneVP::track_path_from_joint(p2s[i], joint[i], prec[i], succ[i]);
        }
    }
#else
    MultiVP::track_path_from_joint(vp_num, p2s.begin(), joint, prec, succ);
#endif

    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}
