
#ifndef MGRAPH_HPP
#define MGRAPH_HPP

#include "dinics-seq.hpp"
#include <unordered_map>
#include <vector>

struct MGraph {
  double bfs_time = 0;
  double dfs_time = 0;
  std::vector<Vertex> vertices;
  std::vector<std::unordered_map<int, Edge>> neighbors;
  std::vector<std::vector<int>> adj_list;
  MGraph(int size);
  bool bfsParallel();
  bool dfsDeadEdge();
  bool isLayerReachable(const Vertex &srcVert, const Vertex &dstVert);
  bool visitVertexParallel(Vertex &srcVert, Vertex &dstVert);
  void increment(int node);
  void dinicsAlgo();
  void reset();
  void addEdge(const Vertex &start, const Vertex &end, int cap);
  void printEdgesVisualized();
  bool bfsStep(std::vector<bool> &visited, std::vector<int> &sizes, int step);
};
#endif // !MGRAPH_HPP
