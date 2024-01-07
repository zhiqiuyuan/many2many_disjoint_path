#include "Baseline.h"
#include "MultiVP.h"

std::vector<std::string> methodcode2str = {
    "maxflow",                // 0
    "ms_uniform",             // 1
    "ms_concentrate_uniform", // 2
    "ms",                     // 3
    "ms_concentrate",         // 4
};

#define CALL(COND, batch_sz, reg_width, ...)                                   \
    COND((batch_sz) == (reg_width))                                            \
    {                                                                          \
        __VA_ARGS__::init_masks();                                             \
        MultiVP ins(vps, problem_num, g, outpf);                               \
        if (method == 1)                                                       \
            total_t = ins.ms_uniform<__VA_ARGS__, __VA_ARGS__>(k);             \
        else if (method == 2)                                                  \
            total_t = ins.ms_concentrate_uniform<__VA_ARGS__, __VA_ARGS__>(k); \
        else if (method == 3)                                                  \
            total_t = ins.ms<__VA_ARGS__, __VA_ARGS__>(k);                     \
        else if (method == 4)                                                  \
            total_t = ins.ms_concentrate<__VA_ARGS__, __VA_ARGS__>(k);         \
    }

int main(int argc, char **argv)
{
    if (argc < 11)
    {
        std::cout << "Wrong usage. Read README.md." << std::endl;
        return -1;
    }

    std::string gfname, vfname, outdir;
    BatchSizeType problem_num = 0, batch_sz = 128;
    int method, k;

    {
        std::unordered_map<std::string, std::string> cmd;
        for (int16_t i = 1; i < argc; i += 2)
            cmd[std::string(argv[i])] = std::string(argv[i + 1]);
        gfname = cmd["-g"];
        vfname = cmd["-v"];
        outdir = cmd["-o"];
        method = std::stoi(cmd["-m"]);
        k = std::stoi(cmd["-k"]);
        if (method >= int(methodcode2str.size()))
        {
            std::cout << "unsupported method: " << method << std::endl;
            return -1;
        }

        if (outdir.empty() == 0 && outdir[outdir.size() - 1] != '/')
            outdir += '/';

        if (cmd.count("-c"))
            problem_num = std::stoul(cmd["-c"]);

        if (cmd.count("-b"))
        {
            batch_sz = std::stoull(cmd["-b"]);
            if (batch_sz < 8)
                batch_sz = 8;
            else
            {
                // round down to 2^(3..7)
                BatchSizeType supported_sz[] = {128, 64, 32, 16, 8};
                int i;
                for (i = 0; supported_sz[i] > batch_sz; ++i)
                    ;
                batch_sz = supported_sz[i];
            }
        }
    }

    std::cout << "[load vertex pairs...]" << std::endl;
    VertexIdType *vps = Graph::load_vps(vfname, problem_num);
    if (vps == 0)
        return -1;

    size_t idx0 = gfname.find_last_of('/') + 1, idx1 = gfname.find_last_of('.');
    std::string dgname = gfname.substr(idx0, idx1 - idx0);
    std::string outpf = outdir + std::to_string(k) + "_paths.txt";
    std::string outtf = outdir + std::to_string(k) + "_time.txt";

    std::ofstream tfout(outtf);
    if (tfout.is_open() == 0)
    {
        std::cout << outtf << " opened failed, write time to this instead: ";
        outtf = dgname + "_time.txt";
        std::cout << outtf << std::endl;
        tfout.open(outtf);
    }
    tfout << "[format]: unit: us" << std::endl;
    tfout << "load graph time" << std::endl;
    tfout << "total compute time" << std::endl;
    tfout << "vertex pair num" << std::endl;
    tfout << "compute time average over vertex pair" << std::endl;

    std::cout << "[load graph...]" << std::endl;
    auto start_time = std::chrono::steady_clock::now();
    Graph g(gfname);
    auto end_time = std::chrono::steady_clock::now();
    tfout << (std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time)).count() << std::endl;
    if (g.is_loaded() == 0)
    {
        std::cout << "graph loaded failed: " << gfname << std::endl;
        return -1;
    }
    std::cout << "\t|V|, |M|= " << g.get_N() << " " << g.get_M() << std::endl;

    std::cout << "[solve...]" << std::endl;
    std::cout << methodcode2str[method] << std::endl;
    long long total_t = 0;
    if (method == 0)
    {
        VertexIdType s, t;
        for (BatchSizeType i = 0; i < problem_num; ++i)
        {
            s = vps[2 * i];
            t = vps[2 * i + 1];
            OneVP ins(s, t, g, outpf);
            if (method == 0)
                total_t +=
                    ins.maxflow(k);

            if (i % 100 == 0)
                std::cout << i << std::endl;
        }
        delete[] vps;
    }
    else
    {
        CALL(if, batch_sz, 128, M128iWrapper)
        CALL(else if, batch_sz, 64, BuiltinUIntWrapper<uint64_t>)
        CALL(else if, batch_sz, 32, BuiltinUIntWrapper<uint32_t>)
        CALL(else if, batch_sz, 16, BuiltinUIntWrapper<uint16_t>)
        CALL(else if, batch_sz, 8, BuiltinUIntWrapper<uint8_t>)
        else
        {
            std::cout << "unsupported batch_sz: " << batch_sz << std::endl;
            return -1;
        }
    }

    tfout << total_t << std::endl;
    tfout << problem_num << std::endl;
    tfout << total_t / (long long)problem_num << std::endl;

    return 0;
}
