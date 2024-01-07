#ifndef _BASELINE_H
#define _BASELINE_H

#include "Graph.h"

struct BrutalNbrUpdate
{
    std::unordered_map<VertexIdType, VertexIdType> vid2conid;
    // CSR
    std::vector<EdgeNumType> in_offset, out_offset;
    std::vector<VertexIdType> in_nbrs, out_nbrs;
};

class OneVP
{
    std::unordered_set<VertexIdType> p1v;
    std::deque<VertexIdType> p2_before;
    std::vector<std::vector<VertexIdType>> ps;

public:
    VertexIdType s, t;
    Graph &g;
    std::string &outpf;

    OneVP(VertexIdType s, VertexIdType t, Graph &g, std::string &outpf) : s(s), t(t), g(g), outpf(outpf) {}

    static void track_path_from_joint(std::deque<VertexIdType> &p1, VertexIdType joint, std::unordered_map<VertexIdType, VertexIdType> &prec, std::unordered_map<VertexIdType, VertexIdType> &succ);

    long long path1(std::vector<VertexIdType> &p1, const std::string &outf = "");

    long long record_p1_for_brutal_path2(std::vector<std::vector<VertexIdType>>::iterator p1_itb, std::vector<std::vector<VertexIdType>>::iterator p1_ite, BrutalNbrUpdate &g_comma);

    long long brutal_path2(BrutalNbrUpdate &g_comma);

    void adjust_p1p2();
    static void adjust_p1p2(VertexIdType N, std::vector<std::vector<VertexIdType>> &ps, std::deque<VertexIdType> &p2_before);

    void output_paths(bool fail);
    static void output_paths(bool fail, VertexIdType s, VertexIdType t, std::vector<std::vector<VertexIdType>> &ps, Graph &g, const std::string &outpf);
    static void output_paths_inner(bool fail, std::ostream &fout, VertexIdType s, VertexIdType t, std::vector<std::vector<VertexIdType>> &ps, Graph &g);

    // method:
    // 0: maxflow
    long long inner_run(int method, int k)
    {
        if (method != 0)
        {
            std::cout << "unsupported method: " << method << std::endl;
            return -1;
        }

        long long duration = 0;
        bool fail = 0;

        ps.emplace_back();
        duration +=
            path1(ps[0]);

        if (ps[0].empty())
            fail = 1;
        else
        {
            for (int kk = 1; kk < k && fail == 0; ++kk)
            {
                if (method == 0)
                {
                    BrutalNbrUpdate g_comma;
                    duration +=
                        record_p1_for_brutal_path2(ps.begin(), ps.end(), g_comma);
                    duration +=
                        brutal_path2(g_comma);

                    auto t1 = std::chrono::steady_clock::now();
                    if (p2_before.empty() == 0)
                        adjust_p1p2();
                    else
                        fail = 1;
                    auto t2 = std::chrono::steady_clock::now();
                    duration += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
                }
            }
        }
        auto t1 = std::chrono::steady_clock::now();
        output_paths(fail);

        auto t2 = std::chrono::steady_clock::now();
        duration += (std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)).count();
        return duration;
    }

    long long maxflow(int k) { return inner_run(0, k); }
};

#endif // _BASELINE_H