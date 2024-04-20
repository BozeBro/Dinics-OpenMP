#pragma once
#include <limits>
#include <unordered_map>
#include <vector>
#define UNSET std::numeric_limits<int>::max()
#define SOURCE 0
#define SINK 1
struct Edge {
  int cap;
  int initial_cap;
};
struct Vertex {
  int index;
  int current_edge = 0;
  int layer = UNSET;
  int parent = -1;
  std::vector<int> layered_dst;

  void reset() {
    this->current_edge = 0;
    this->layer = UNSET;
    this->parent = -1;
    this->layered_dst.clear();
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
  void addEdge(const Vertex &start, const Vertex &end, int cap);
  void reset();
};
