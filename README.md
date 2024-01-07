## build
```bash
mkdir build
cd build
cmake ..
make vdp -j
```
## run
Parameter:
- `-g`: path to the graph file
- `-v`: path to the vertex pair file
- `-o`: output folder
- `-m`: method number
- `-k`: parameter k
- `[-c]`: number of vertex pairs (optional, default to all vertex pairs in -v)
- `[-b]`: batch_size (optional, default to 128)

Method number:
- 0: maxflow
- 1: ms_uniform
- 2: ms_concentrate_uniform
- 3: ms
- 4: ms_concentrate

For example:
```bash
build/vdp -g test_data/graph.txt -v test_data/vertex_pairs.txt -o test_data -m 4 -k 2
```