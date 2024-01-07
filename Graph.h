#ifndef _GRAPH_H
#define _GRAPH_H

#include "tools.h"

class MmapedFile
{
    int fd;
    long long fsz;
    void *mapping;

public:
    MmapedFile(const char *f, int flags)
    {
        fd = open(f, flags);
        if (fd < 0)
        {
            std::cout << "file open failed: " << f << std::endl;
            return;
        }
        fsz = lseek(fd, 0, SEEK_END);
        if (unlikely(fsz == -1))
        {
            close(fd);
            fd = -1;
            std::cout << "file size get failed: " << f << std::endl;
            return;
        }

        mapping = mmap(0, fsz, PROT_READ, MAP_PRIVATE, fd, 0);
        if (mapping == MAP_FAILED)
        {
            close(fd);
            fd = -1;
            std::cout << "mmap failed: " << f << std::endl;
            return;
        }

        close(fd);
    }

    MmapedFile(MmapedFile &&other) : fd(other.fd), fsz(other.fsz), mapping(other.mapping)
    {
        other.fd = -1;
        other.fsz = 0;
        other.mapping = 0;
    }
    ~MmapedFile()
    {
        if (fd >= 0)
        {
            munmap(mapping, fsz);
        }
    }

    bool is_mapped() { return fd >= 0; }
    long long get_fsz() { return fsz; }
    void *get_mapping() { return mapping; }
};

class Graph
{
    VertexIdType N;
    EdgeNumType M;
    // CSR
    VertexIdType *nbrs;
    EdgeNumType *offset;

    static void generate_rand_vertex_pairs(VertexIdType N, std::function<VertexIdType(VertexIdType)> map_vid, const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num);

public:
    Graph() : N(0), M(0), nbrs(0), offset(0)
    {
    }
    Graph(const std::string &gf) : N(0), M(0), nbrs(0), offset(0)
    {
        load(gf);
    }
    ~Graph()
    {
        if (nbrs)
        {
            delete[] nbrs;
            delete[] offset;
        }
    }

    void load(const std::string &gf);
    bool is_loaded() const { return nbrs != 0; }

    void print_nbrs()
    {
        for (VertexIdType i = 0; i < N; ++i)
        {
            std::cout << i << ": ";
            for (EdgeNumType idx = offset[i]; idx < offset[i + 1]; ++idx)
                std::cout << nbrs[idx] << " ";
            std::cout << std::endl;
        }
    }

    VertexIdType get_N() { return N; }
    EdgeNumType get_M() { return M; }
    VertexIdType *get_nbrs() { return nbrs; }
    EdgeNumType *get_offset() { return offset; }

    static VertexIdType max_vid(const std::string &f);
    static VertexIdType *load_vps(const std::string &vfname, BatchSizeType &problem_num);

    static void generate_rand_vertex_pairs(const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num);
    static void generate_rand_vertex_pairs_above_degree_bound(EdgeNumType degree_lower_bound, const std::string &f, const std::string &vp_output_dir, BatchSizeType vp_num);
};

#endif //_GRAPH_H