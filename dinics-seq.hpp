#ifndef DINICS_SEQ_HPP
#define DINICS_SEQ_HPP

#include <limits>
#include <ostream>
#include <unordered_map>
#include <vector>
#define UNSET std::numeric_limits<int>::max()
#define SOURCE 0
#define SINK 1
// #define MY_PRINT_ENABLED
#ifdef MY_PRINT_ENABLED
// If it's defined, define PRINTF as printf
#define PRINTF(format, ...) printf(format, ##__VA_ARGS__)
#else
// If it's not defined, define PRINTF as an empty macro
#define PRINTF(format, ...) ((void)0)
#endif

#define MAX_THREADS 8

struct Edge {
  int cap;
  int initial_cap;
};
struct Vertex {
  int index;
  int current_edge = 0;
  int current_edge_array = 0;
  int layer = UNSET;
  int parent = -1;
  std::vector<int> layered_dst;
  std::vector<int> layered_dst_array[MAX_THREADS];

  void reset() {
    this->current_edge = 0;
    this->current_edge_array = 0;
    this->layer = UNSET;
    this->parent = -1;
    this->layered_dst.clear();
    for (int i = 0; i < MAX_THREADS; i++) {
      this->layered_dst_array[i] = std::vector<int>();
    }
  }
};
namespace std {
template <> struct hash<Vertex> {
  std::size_t operator()(const Vertex &obj) const {
    // index is the unique identifier anyway
    return std::hash<int>()(obj.index);
  }
};
} // namespace std

struct Graph {
  double bfs_time = 0;
  double dfs_time = 0;

  std::vector<Vertex> vertices;
  std::vector<std::unordered_map<int, Edge>> neighbors;

  Graph(int size);
  void printEdges();
  void printEdgesVisualized();
  void validate();
  bool bfsParallel();
  bool isLayerReachable(const Vertex &srcVert, const Vertex &dstVert);
  bool visitVertexParallel(Vertex &srcVert, Vertex &dstVert);
  bool visitVertex(Vertex &srcVert, Vertex &dstVert);
  bool bfs();
  void increment(int node);
  bool dfsDeadEdge();
  bool dfs();

  inline void addEdge(const Vertex &start, const Vertex &end, int cap) {
    // Max prevent capacity override if there are 2 node cycles
    int cap_value = std::max(neighbors[start.index][end.index].cap, cap);
    neighbors[start.index][end.index] = {cap_value, cap_value};
    cap_value = std::max(neighbors[end.index][start.index].cap, 0);
    neighbors[end.index][start.index] = {cap_value, cap_value};
    // neighbors[end.index][start.index].cap =
    //     std::max(neighbors[end.index][start.index].cap, 0);
  }

  void reset();
  void dinicsAlgo();
  friend std::ostream &operator<<(std::ostream &os, const Graph &graph);
};

bool operator==(const Vertex &a, const Vertex &b);
// {
//   return a.index == b.index && a.current_edge == b.current_edge &&
//          a.layer == b.layer;
// }
bool operator!=(const Vertex &a, const Vertex &b); /* { return !(a == b); } */
std::ostream &operator<<(std::ostream &os, const Graph &graph);
#endif // !DINICS_SEQ_HPP
