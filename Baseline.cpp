#include "Baseline.h"

void OneVP::track_path_from_joint(std::deque<VertexIdType> &p1, VertexIdType joint, std::unordered_map<VertexIdType, VertexIdType> &prec, std::unordered_map<VertexIdType, VertexIdType> &succ)
{
    VertexIdType v = joint, n;
    std::unordered_map<VertexIdType, VertexIdType>::iterator ed = prec.end(), it;
    p1.push_front(v);
    while ((it = prec.find(v)) != ed)
    {
        n = it->second;
        p1.push_front(n);
        v = n;
    }
    ed = succ.end();
    v = joint;
    while ((it = succ.find(v)) != ed)
    {
        n = it->second;
        p1.push_back(n);
        v = n;
    }
}

long long OneVP::path1(std::vector<VertexIdType> &p1, const std::string &outf)
{
    auto t1 = std::chrono::steady_clock::now();

    VertexIdType *nbrs = g.get_nbrs();
    EdgeNumType *offset = g.get_offset();

    std::queue<VertexIdType> s_q, t_q;
    std::unordered_set<VertexIdType> s_seen, t_seen;
    std::unordered_map<VertexIdType, VertexIdType> prec, succ;
    s_q.push(s);
    t_q.push(t);
    s_seen.insert(s);
    t_seen.insert(t);

    VertexIdType v, u;
    VertexIdType joint = (VertexIdType)-1ull;
    size_t qsz;
    while (joint == (VertexIdType)-1ull && s_q.empty() == 0 && t_q.empty() == 0)
    {
        // s_q
        qsz = s_q.size();
        while (joint == (VertexIdType)-1ull && qsz)
        {
            --qsz;
            v = s_q.front();
            s_q.pop();

            for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
            {
                u = nbrs[i];
                if (s_seen.find(u) == s_seen.end())
                {
                    s_seen.insert(u);
                    prec[u] = v;
                    if (t_seen.find(u) != t_seen.end())
                    {
                        joint = u;
                        break;
                    }
                    s_q.push(u);
                }
            }
        }

        if (~joint)
            break;

        // t_q
        qsz = t_q.size();
        while (joint == (VertexIdType)-1ull && qsz)
        {
            --qsz;
            v = t_q.front();
            t_q.pop();

            for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
            {
                u = nbrs[i];
                if (t_seen.find(u) == t_seen.end())
                {
                    t_seen.insert(u);
                    succ[u] = v;
                    if (s_seen.find(u) != s_seen.end())
                    {
                        joint = u;
                        break;
                    }
                    t_q.push(u);
                }
            }
        }
    }

    if (~joint)
    {
        std::deque<VertexIdType> p1d;
        track_path_from_joint(p1d, joint, prec, succ);
        std::copy(p1d.begin(), p1d.end(), std::back_inserter(p1));

        if (outf.empty() == 0)
        {
            std::ofstream fout(outf, std::ios::app);
            if (fout.is_open() == 0)
            {
                std::cout << "file opened failed: " << outf << std::endl;
                return 0;
            }
            fout << "\t" << s << " " << t << std::endl;
            for (VertexIdType v : p1)
                fout << v << " ";
            fout << std::endl;
        }
    }

    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

long long OneVP::record_p1_for_brutal_path2(std::vector<std::vector<VertexIdType>>::iterator p1_itb, std::vector<std::vector<VertexIdType>>::iterator p1_ite, BrutalNbrUpdate &g_comma)
{
    auto t1 = std::chrono::steady_clock::now();

    VertexIdType N = g.get_N();
    EdgeNumType *offset = g.get_offset();
    VertexIdType *nbrs = g.get_nbrs();

    VertexIdType inner_id = 0;
    std::vector<std::vector<VertexIdType>> in_nbrs, out_nbrs;
    auto &vid2conid = g_comma.vid2conid;
    auto map_vertex = [&vid2conid, N, &inner_id, &in_nbrs, &out_nbrs](VertexIdType v, bool map_double) -> VertexIdType
    {
        VertexIdType inner_v;
        std::unordered_map<VertexIdType, VertexIdType>::iterator it;
        if ((it = vid2conid.find(v)) != vid2conid.end())
            inner_v = it->second;
        else
        {
            inner_v = inner_id++;
            vid2conid[v] = inner_v;
            if (map_double)
                vid2conid[v + N] = inner_id++;
            in_nbrs.resize(inner_id);
            out_nbrs.resize(inner_id);
        }
        return inner_v;
    };

    p1v.clear();
    for (auto it = p1_itb; it < p1_ite; ++it)
        for (VertexIdType v : *it)
            p1v.insert(v);

    VertexIdType pre = s, inner_pre, n;
    // s
    {
        pre = s;
        inner_pre = map_vertex(pre, 0);
        std::unordered_set<VertexIdType> s_next_hops;
        for (auto it = p1_itb; it < p1_ite; ++it)
            s_next_hops.insert((*it)[1]);
        for (EdgeNumType k = offset[pre], ked = offset[pre + 1]; k < ked; ++k)
        {
            n = nbrs[k];
            if (n != s && n != t && p1v.find(n) != p1v.end())
                in_nbrs[inner_pre].push_back(n + N);
            else
                in_nbrs[inner_pre].push_back(n);
            if (s_next_hops.find(n) == s_next_hops.end())
                out_nbrs[inner_pre].push_back(n);
        }
    }
    // t
    {
        pre = t;
        inner_pre = map_vertex(pre, 0);
        std::unordered_set<VertexIdType> t_pre_hops;
        for (auto it = p1_itb; it < p1_ite; ++it)
            t_pre_hops.insert(*((it->end()) - 2));
        for (EdgeNumType k = offset[pre], ked = offset[pre + 1]; k < ked; ++k)
        {
            n = nbrs[k];
            if (t_pre_hops.find(n) == t_pre_hops.end())
            {
                if (n != s && n != t && p1v.find(n) != p1v.end())
                    in_nbrs[inner_pre].push_back(n + N);
                else
                    in_nbrs[inner_pre].push_back(n);
                out_nbrs[inner_pre].push_back(n);
            }
            else
            {
                if (n != s)
                    out_nbrs[inner_pre].push_back(n + N);
                else
                    out_nbrs[inner_pre].push_back(n);
            }
        }
    }
    // middle vertex in p1
    VertexIdType curr, prepre = (VertexIdType)-1ull;
    for (auto it = p1_itb; it < p1_ite; ++it)
    {
        pre = s;
        prepre = (VertexIdType)-1ull;
        auto &p1 = *it;
        for (size_t j = 1, sz = p1.size(); j < sz; ++j)
        {
            curr = p1[j];

            if (j > 1)
            {
                inner_pre = map_vertex(pre, 1);
                for (EdgeNumType k = offset[pre], ked = offset[pre + 1]; k < ked; ++k)
                {
                    n = nbrs[k];
                    if (n != prepre)
                    {
                        if (n != s && n != t && p1v.find(n) != p1v.end())
                            in_nbrs[inner_pre].push_back(n + N);
                        else
                            in_nbrs[inner_pre].push_back(n);
                    }
                    else
                        in_nbrs[inner_pre].push_back(pre + N);
                }
                if (prepre != s)
                    out_nbrs[inner_pre].push_back(prepre + N);
                else
                    out_nbrs[inner_pre].push_back(prepre);
                ++inner_pre;
                for (EdgeNumType k = offset[pre], ked = offset[pre + 1]; k < ked; ++k)
                {
                    n = nbrs[k];
                    if (n != curr)
                        out_nbrs[inner_pre].push_back(n);
                    else
                        out_nbrs[inner_pre].push_back(pre);
                }
                in_nbrs[inner_pre].push_back(curr);
            }

            prepre = pre;
            pre = curr;
        }
    }

    // -> CSR
    EdgeNumType e = 0;
    for (VertexIdType j = 0; j < inner_id; ++j)
    {
        g_comma.in_offset.push_back(e);
        e += in_nbrs[j].size();
        for (VertexIdType v : in_nbrs[j])
            g_comma.in_nbrs.push_back(v);
    }
    g_comma.in_offset.push_back(e);
    e = 0;
    for (VertexIdType j = 0; j < inner_id; ++j)
    {
        g_comma.out_offset.push_back(e);
        e += out_nbrs[j].size();
        for (VertexIdType v : out_nbrs[j])
            g_comma.out_nbrs.push_back(v);
    }
    g_comma.out_offset.push_back(e);

    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

long long OneVP::brutal_path2(BrutalNbrUpdate &g_comma)
{
    auto t1 = std::chrono::steady_clock::now();

    p2_before.clear();
    VertexIdType *nbrs = g.get_nbrs();
    EdgeNumType *offset = g.get_offset();
    VertexIdType N = g.get_N();

    std::queue<VertexIdType> s_q, t_q;
    std::unordered_set<VertexIdType> s_seen, t_seen;
    std::unordered_map<VertexIdType, VertexIdType> prec, succ;
    s_q.push(s);
    t_q.push(t);
    s_seen.insert(s);
    t_seen.insert(t);

    VertexIdType v, u, v0;
    VertexIdType joint = (VertexIdType)-1ull;
    size_t qsz;
    while (joint == (VertexIdType)-1ull && s_q.empty() == 0 && t_q.empty() == 0)
    {
        // s_q
        qsz = s_q.size();
        while (joint == (VertexIdType)-1ull && qsz)
        {
            --qsz;
            v = s_q.front();
            s_q.pop();

            if (v >= N)
                v0 = v - N;
            else
                v0 = v;
            if (p1v.find(v0) != p1v.end())
            {
                auto &vid2conid = g_comma.vid2conid;
                VertexIdType inner_v = vid2conid.at(v);
                for (EdgeNumType i = g_comma.out_offset[inner_v], ied = g_comma.out_offset[inner_v + 1]; i < ied; ++i)
                {
                    u = g_comma.out_nbrs[i];
                    if (s_seen.find(u) == s_seen.end())
                    {
                        s_seen.insert(u);
                        prec[u] = v;
                        if (t_seen.find(u) != t_seen.end())
                        {
                            joint = u;
                            break;
                        }
                        s_q.push(u);
                    }
                }
            }
            else
            {
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    if (s_seen.find(u) == s_seen.end())
                    {
                        s_seen.insert(u);
                        prec[u] = v;
                        if (t_seen.find(u) != t_seen.end())
                        {
                            joint = u;
                            break;
                        }
                        s_q.push(u);
                    }
                }
            }
        }

        if (~joint)
            break;

        // t_q
        qsz = t_q.size();
        while (joint == (VertexIdType)-1ull && qsz)
        {
            --qsz;
            v = t_q.front();
            t_q.pop();

            if (v >= N)
                v0 = v - N;
            else
                v0 = v;
            if (p1v.find(v0) != p1v.end())
            {
                auto &vid2conid = g_comma.vid2conid;
                VertexIdType inner_v = vid2conid.at(v);
                for (EdgeNumType i = g_comma.in_offset[inner_v], ied = g_comma.in_offset[inner_v + 1]; i < ied; ++i)
                {
                    u = g_comma.in_nbrs[i];
                    if (t_seen.find(u) == t_seen.end())
                    {
                        t_seen.insert(u);
                        succ[u] = v;
                        if (s_seen.find(u) != s_seen.end())
                        {
                            joint = u;
                            break;
                        }
                        t_q.push(u);
                    }
                }
            }
            else
            {
                for (EdgeNumType i = offset[v], ied = offset[v + 1]; i < ied; ++i)
                {
                    u = nbrs[i];
                    if (u != s && u != t && p1v.find(u) != p1v.end())
                        u += N;
                    if (t_seen.find(u) == t_seen.end())
                    {
                        t_seen.insert(u);
                        succ[u] = v;
                        if (s_seen.find(u) != s_seen.end())
                        {
                            joint = u;
                            break;
                        }
                        t_q.push(u);
                    }
                }
            }
        }
    }

    if (~joint)
    {
        track_path_from_joint(p2_before, joint, prec, succ);
    }

    auto t2 = std::chrono::steady_clock::now();
    return (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
}

void OneVP::adjust_p1p2()
{
    adjust_p1p2(g.get_N(), ps, p2_before);
}

void OneVP::adjust_p1p2(VertexIdType N, std::vector<std::vector<VertexIdType>> &ps, std::deque<VertexIdType> &p2_before)
{
    VertexIdType s = p2_before.front();
    std::unordered_map<VertexIdType, VertexIdType> next_hop;
    std::vector<VertexIdType> s_next_hop;
    s_next_hop.reserve(ps.size() + 1);
    VertexIdType v, n;
    for (auto &p1_before : ps)
    {
        for (size_t i = 1, sz = p1_before.size(); i < sz; ++i)
        {
            n = p1_before[i];
            if (i > 1)
            {
                next_hop[v] = v + N;
                next_hop[v + N] = n;
            }
            else
                s_next_hop.push_back(n);
            v = n;
        }
    }
    std::unordered_map<VertexIdType, VertexIdType>::iterator it;
    s_next_hop.push_back(p2_before[1]);
    n = p2_before[1];
    for (size_t i = 2, sz = p2_before.size(); i < sz; ++i)
    {
        v = p2_before[i];
        if (((it = next_hop.find(v)) != next_hop.end()) && (it->second == n))
        {
            next_hop.erase(it);
        }
        else
            next_hop[n] = v;
        n = v;
    }

    int k = int(ps.size() + 1);
    ps.clear();
    ps.resize(k);
    auto ed = next_hop.end();
    for (int j = 0; j < k; ++j)
    {
        ps[j].push_back(s);
        v = s_next_hop[j];
        ps[j].push_back(v);
        while ((it = next_hop.find(v)) != ed)
        {
            v = it->second;
            if (v < N)
                ps[j].push_back(v);
        }
    }
}

void OneVP::output_paths(bool fail)
{
    output_paths(fail, s, t, ps, g, outpf);
}
void OneVP::output_paths(bool fail, VertexIdType s, VertexIdType t, std::vector<std::vector<VertexIdType>> &ps, Graph &g, const std::string &outpf)
{
    std::ofstream fout(outpf, std::ios::app);
    if (fout.is_open())
        output_paths_inner(fail, fout, s, t, ps, g);
    else
        output_paths_inner(fail, std::cout, s, t, ps, g);
}

void OneVP::output_paths_inner(bool fail, std::ostream &fout, VertexIdType s, VertexIdType t, std::vector<std::vector<VertexIdType>> &ps, Graph &g)
{
    fout << "\t" << s << " " << t << std::endl;
    if (fail)
    {
        fout << "no solution" << std::endl;
    }
    else
    {
        int k = int(ps.size());
        for (int j = 0; j < k; ++j)
        {
            for (auto v : ps[j])
                fout << v << " ";
            fout << std::endl;
        }
    }
}
