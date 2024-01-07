#include "Graph.h"

void Graph::load(const std::string &gf)
{
    MmapedFile mp(gf.c_str(), O_RDONLY);
    if (mp.is_mapped() == 0)
        return;
    uint8_t *addr = (uint8_t *)mp.get_mapping();
    long long fsz = mp.get_fsz();

    std::vector<std::vector<VertexIdType>> vecnbrs;
    auto insert_nbr_inner = [&vecnbrs](VertexIdType v, VertexIdType nbr)
    {
        if (vecnbrs.size() <= v)
            vecnbrs.resize(v + 1);
        vecnbrs[v].push_back(nbr);
    };
    auto insert_nbr = [&insert_nbr_inner](VertexIdType v1, VertexIdType v2)
    {
        insert_nbr_inner(v1, v2);
        insert_nbr_inner(v2, v1);
    };

    EdgeNumType M0 = 0;
    long long i = 0;
    VertexIdType v1 = (VertexIdType)-1ull;
    while (i < fsz && addr[i] == '#')
    {
        while (i < fsz && addr[i] != '\n' && addr[i] != '\r')
            ++i;
        while (i < fsz && (addr[i] == '\n' || addr[i] == '\r'))
            ++i;
    }

    while (i < fsz && (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r'))
        ++i;
    bool already_read_first_v = 0;
    while (i < fsz)
    {
        VertexIdType id = 0;
        while (i < fsz && addr[i] >= '0' && addr[i] <= '9')
        {
            id = id * 10 + (addr[i] - '0');
            ++i;
        }
        if (i >= fsz)
        {
            if (v1 != id)
            {
                insert_nbr(v1, id);
                ++M0;
            }
            break;
        }
        if (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r')
        {
            if (already_read_first_v == 0)
            {
                v1 = id;
                already_read_first_v = 1;
            }
            else
            {
                if (v1 != id)
                {
                    insert_nbr(v1, id);
                    ++M0;
                }
                already_read_first_v = 0;
            }
        }
        else
        {
            std::cout << "wrong file format. unexpected char: " << addr[i] << std::endl;
            nbrs = 0;
            return;
        }
        ++i;
        while (i < fsz && (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r'))
            ++i;
    }

    N = vecnbrs.size();
    M0 *= 2;
    nbrs = new VertexIdType[M0];
    this->M = M0;
    if (unlikely(nbrs == NULL))
    {
        std::cout << "nbrs allocated failed. edges num(before undirecd dup):" << M0 << std::endl;
        nbrs = 0;
        return;
    }
    offset = new EdgeNumType[N + 1];
    if (unlikely(offset == NULL))
    {
        std::cout << "offset allocated failed. edges num(before undirecd dup):" << M0 << " vertex num:" << N << std::endl;
        nbrs = 0;
        return;
    }

    M0 = 0;
    for (VertexIdType v = 0; v < N; ++v)
    {
        offset[v] = M0;
        for (VertexIdType n : vecnbrs[v])
            nbrs[M0++] = n;
    }
    offset[N] = M0;
}

VertexIdType *Graph::load_vps(const std::string &vfname, BatchSizeType &problem_num)
{
    std::ifstream fin(vfname);
    if (fin.is_open() == 0)
    {
        std::cout << vfname << " opened failed!" << std::endl;
        return 0;
    }

    BatchSizeType totol_num;
    fin >> totol_num;
    {
        std::string nouse;
        getline(fin, nouse);
    }
    if (problem_num)
        problem_num = std::min(problem_num, totol_num);
    else
        problem_num = totol_num;
    VertexIdType *vps = new VertexIdType[problem_num * 2];

    VertexIdType s, t;
    for (BatchSizeType i = 0; (i < problem_num) && fin >> s >> t; ++i)
    {
        vps[2 * i] = s;
        vps[2 * i + 1] = t;
    }
    return vps;
}

VertexIdType Graph::max_vid(const std::string &f)
{
    MmapedFile mp(f.c_str(), O_RDONLY);
    if (mp.is_mapped() == 0)
    {
        return VertexIdType(-1ull);
    }
    uint8_t *addr = (uint8_t *)mp.get_mapping();
    long long fsz = mp.get_fsz();

    auto return_max_v = [](VertexIdType v, VertexIdType max_v)
    {
        if (v > max_v || max_v == (VertexIdType)-1ull)
            return v;
        return max_v;
    };

    long long i = 0;
    VertexIdType v1 = (VertexIdType)-1ull, max_v = v1;
    while (i < fsz && addr[i] == '#')
    {
        while (i < fsz && addr[i] != '\n' && addr[i] != '\r')
            ++i;
        while (i < fsz && (addr[i] == '\n' || addr[i] == '\r'))
            ++i;
    }

    while (i < fsz && (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r'))
        ++i;
    bool already_read_first_v = 0;
    while (i < fsz)
    {
        VertexIdType id = 0;
        while (i < fsz && addr[i] >= '0' && addr[i] <= '9')
        {
            id = id * 10 + (addr[i] - '0');
            ++i;
        }
        if (i >= fsz)
        {
            max_v = return_max_v(v1, max_v);
            max_v = return_max_v(id, max_v);
            break;
        }
        if (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r')
        {
            if (already_read_first_v == 0)
            {
                v1 = id;
                already_read_first_v = 1;
            }
            else
            {
                max_v = return_max_v(v1, max_v);
                max_v = return_max_v(id, max_v);
                already_read_first_v = 0;
            }
        }
        else
        {
            std::cout << "wrong file format. unexpected char: " << addr[i] << std::endl;
            return VertexIdType(-1ull);
        }
        ++i;
        while (i < fsz && (addr[i] == ' ' || addr[i] == '\t' || addr[i] == '\n' || addr[i] == '\r'))
            ++i;
    }

    return max_v;
}

void Graph::generate_rand_vertex_pairs(const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num)
{
    VertexIdType N = max_vid(f);
    if (~N)
    {
        auto identical_map = [](VertexIdType v)
        { return v; };
        generate_rand_vertex_pairs(N, identical_map, f, vp_output_dir, vp_num);
    }
    else
        std::cout << "read N failed! stop gen" << std::endl;
}
void Graph::generate_rand_vertex_pairs_above_degree_bound(EdgeNumType degree_lower_bound, const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num)
{
    std::cout << "load graph..." << std::endl;
    Graph g(f);
    std::cout << "load graph done" << std::endl;
    if (g.is_loaded())
    {
        std::vector<VertexIdType> newid2vid;
        EdgeNumType *offset = g.offset;
        for (VertexIdType v = 0, N = g.N; v < N; ++v)
            if (offset[v + 1] - offset[v] >= degree_lower_bound)
                newid2vid.push_back(v);
        VertexIdType rest_v_num = VertexIdType(newid2vid.size());
        if (rest_v_num * (rest_v_num - 1) / 2 < vp_num)
        {
            std::cout << "vertices with degree >= " << degree_lower_bound << " is too few! only have " << rest_v_num << " of them. cannot make up " << vp_num << " different vps! (at most make up " << rest_v_num * (rest_v_num - 1) / 2 << ") stop gen" << std::endl;
            return;
        }
        auto map2orivid = [&newid2vid, rest_v_num](VertexIdType v) -> VertexIdType
        {
            if (v < rest_v_num)
                return newid2vid[v];
            else
            {
                std::cout << "invalid v: " << v << ". >=" << newid2vid.size() << std::endl;
                return (VertexIdType)-1ull;
            }
        };
        generate_rand_vertex_pairs(rest_v_num, map2orivid, f, vp_output_dir, vp_num);
    }
    else
        std::cout << "graph loaded failed! stop gen" << std::endl;
}
void Graph::generate_rand_vertex_pairs(VertexIdType N, std::function<VertexIdType(VertexIdType)> map_vid, const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num)
{
    size_t i1 = f.find_last_of('.'), i2 = f.find_last_of('/');
    std::string g = f.substr(i2 + 1, i1 - i2 - 1);
    std::string vp_f = vp_output_dir + '/' + g + ".txt";

    std::ofstream fout(vp_f);
    fout << vp_num << std::endl;
    std::unordered_set<VertexIdType> already;
    VertexIdType s, t;
    auto gen_vp = [N, &s, &t]()
    {
        s = Rand(N);
        t = Rand(N);
        while (t == s)
            t = Rand(N);
    };
    for (BatchSizeType i = 0; i < vp_num; ++i)
    {
        gen_vp();
        while (already.count(s * N + t))
            gen_vp();
        already.insert(s * N + t);
        fout << map_vid(s) << ' ' << map_vid(t) << '\n';
        if (i % 100 == 0)
            std::cout << i << " ";
    }
    std::cout << "\nwrite vp to file " << vp_f << std::endl;
}
