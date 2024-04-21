#include "mgraph.hpp"
#include "dinics-seq.hpp"
#include <cassert>
#include <chrono>
// #include <omp.h>
#include <stack>
#include <vector>
void MGraph::printEdgesVisualized() {
  PRINTF("digraph G {\n");
  for (int i = 0; i < this->vertices.size(); i++) {
    for (int neigh = 0; neigh < this->vertices.size(); neigh++) {
      Edge edge = this->neighbors[i][neigh];
      if (edge.initial_cap > 0)
        PRINTF("%d -> %d [label = \"%d/%d\"];\n", i, neigh,
               std::max(0, edge.initial_cap - edge.cap), edge.initial_cap);
    }
  }
  PRINTF("}\n");
}

MGraph::MGraph(int size)
    : vertices(size), neighbors(size), adj_list(size) {
  for (int i = 0; i < size; i++) {
    vertices[i].index = i;
  }
}
void MGraph::increment(int node) {
  if (node >= 0)
    this->vertices[node].current_edge++;
}

bool MGraph::isLayerReachable(const Vertex &srcVert, const Vertex &dstVert) {
  if (dstVert.layer <= srcVert.layer ||
      this->neighbors[srcVert.index][dstVert.index].cap == 0) {
    PRINTF("dst layer: %d src layer: %d capacity %d\n", dstVert.layer,
           srcVert.layer, this->neighbors[srcVert.index][dstVert.index].cap);
    // PRINTF("capacity %d\n",
    // this->neighbors[srcVert.index][dstVert.index].cap);
    PRINTF("Fail here %d->%d\n", srcVert.index, dstVert.index);
    return false;
  }
  return true;
}

bool MGraph::bfsStep(std::vector<bool> &visited, int step) {
  bool sz = false;
  #pragma omp parallel for reduction(| : sz)
  for (int i = 0; i < visited.size(); i++) {
    if (this->vertices[i].layer != step)
      continue;
    for (int neigh = 0; neigh < this->adj_list[i].size(); neigh++) {

      int neighInd = this->adj_list[i][neigh];
      Vertex &dstVert = this->vertices[neighInd];
      Vertex &srcVert = this->vertices[i];
      if (visited[neighInd] || !isLayerReachable(srcVert, dstVert)) {
        continue;
      }
      
      if (visited[neighInd]) {
        continue;
      }
      visited[neighInd] = true;
      dstVert.layer = step + 1;
      sz |= true;
    }
  }
  return sz;
}
bool MGraph::bfsParallel() {
  auto start = std::chrono::steady_clock::now();
  bool foundSink = false;
  std::vector<bool> visited(this->vertices.size(), false);
  visited[SOURCE] = true;
  this->vertices[SOURCE].layer = 0;
  bool seen = true;
  int step = 0;
  while (seen && !visited[SINK]) {

    seen = bfsStep(visited, step++);
  }
  auto end = std::chrono::steady_clock::now();
  bfs_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();

  return visited[SINK];
}
bool MGraph::visitVertexParallel(Vertex &srcVert, Vertex &dstVert) {

  if (dstVert.layer == UNSET) {
    // frontier.push(dstVert.index);
    return true;
  } else {
    if (dstVert.layer != srcVert.layer + 1) {
      PRINTF("index: %d %d\n", dstVert.index, dstVert.layer);
      abort();
    }
  }
  return false;
}
bool MGraph::dfsDeadEdge() {
  auto start = std::chrono::steady_clock::now();
  std::vector<bool> visited(this->vertices.size(), false);
  std::stack<int> stack;
  stack.push(SOURCE);
  while (!stack.empty()) {
    int nodeInd = stack.top();
    Vertex &srcVert = this->vertices[nodeInd];
    visited[nodeInd] = true;
    auto &neighborEdges = this->adj_list[nodeInd];
    if (srcVert.current_edge == neighborEdges.size()) {
      stack.pop();
      increment(srcVert.parent);
      continue;
    }
    assert(srcVert.current_edge < neighborEdges.size());
    int neigh = this->adj_list[nodeInd][srcVert.current_edge];
    Vertex &dstVert = this->vertices[neigh];
    // Vertex &dstVert = this->vertices[neigh];
    if (dstVert.layer <= srcVert.layer || visited[neigh] ||
        this->neighbors[nodeInd][neigh].cap == 0) {
      increment(nodeInd);
      continue;
    }
    this->vertices[neigh].parent = nodeInd;
    if (neigh == SINK) {
      auto end = std::chrono::steady_clock::now();
      dfs_time +=
          std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
              .count();
      return true;
    }
    stack.push(neigh);
  }
  auto end = std::chrono::steady_clock::now();
  dfs_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();
  return false;
}
void MGraph::dinicsAlgo() {
  while (bfsParallel()) {
    PRINTF("Did a BFS\n");

    // printEdges();
    while (dfsDeadEdge()) {
      // PRINTF("Finished dfs iteration\n");
      Vertex &dstVert = vertices[SINK];
      int minCapacity = UNSET;
      for (Vertex cur = dstVert; cur != vertices[SOURCE];
           cur = vertices[cur.parent]) {
        minCapacity =
            std::min(minCapacity, neighbors[cur.parent][cur.index].cap);
      }
      for (Vertex cur = dstVert; cur != vertices[SOURCE];
           cur = vertices[cur.parent]) {
        neighbors[cur.parent][cur.index].cap -= minCapacity;
        neighbors[cur.index][cur.parent].cap += minCapacity;
      }
    }
    reset();
  }
}
void MGraph::reset() {
  for (Vertex &v : this->vertices) {
    v.reset();
  }
}
void MGraph::addEdge(const Vertex &start, const Vertex &end, int cap) {
  // Max prevent capacity override if there are 2 node cycles
  this->adj_list[start.index].push_back(end.index);
  int cap_value = std::max(neighbors[start.index][end.index].cap, cap);
  neighbors[start.index][end.index] = {cap_value, cap_value};
  cap_value = std::max(neighbors[end.index][start.index].cap, 0);
  neighbors[end.index][start.index] = {cap_value, cap_value};
  // neighbors[end.index][start.index].cap =
  //     std::max(neighbors[end.index][start.index].cap, 0);
}
