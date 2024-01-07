#include "Graph.h"

void generate_rand_vertex_pairs(int argc, char **argv)
{
    std::string f = argv[1];
    std::string vp_output_dir = argv[2];
    BatchSizeType vp_num = (BatchSizeType)std::stoull(argv[3]);
    EdgeNumType low_deg_bound = (EdgeNumType)std::stoull(argv[4]);
    Graph::generate_rand_vertex_pairs_above_degree_bound(low_deg_bound, f, vp_output_dir, vp_num);
}

// bin choice [args]
int main(int argc, char **argv)
{
    int choice = std::stoi(argv[1]);
    argc--;
    argv++;
    if (choice == 2)
        // bin 2 f vp_output_dir vp_num low_deg_bound
        generate_rand_vertex_pairs(argc, argv);
    else
        std::cout << "unsupported choice: " << choice << std::endl;
    return 0;
}