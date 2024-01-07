typedef uint32_t VertexIdType;
typedef int32_t SignedVertexIdType;
typedef uint64_t EdgeNumType;
typedef uint32_t BatchSizeType;

#define MERGE_PATH_RECORE
// if set, merge path record for vp in same batch
// else, record path independently for each vp (just as maxflow)
