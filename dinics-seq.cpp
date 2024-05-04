#include "dinics-seq.hpp"
// create aug graph

// find a layered graph (BFS)

// run lots of DFS to find blocking flow in layered graph

/*
Be able to change flow of Forward and Backward Edge Easily
Be able to run BFS and DFS through it in a simple way

*/
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <omp.h>
#include <ostream>
#include <queue>
#include <stack>
#include <unistd.h>
#include <unordered_map>
#include <vector>

bool operator==(const Vertex &a, const Vertex &b) {
  return a.index == b.index && a.current_edge == b.current_edge &&
         a.layer == b.layer;
}
bool operator!=(const Vertex &a, const Vertex &b) { return !(a == b); }

Graph::Graph(int size) : vertices(size), neighbors(size) {
  for (int i = 0; i < size; i++) {
    vertices[i].index = i;
  }
}

void Graph::printEdges() {
  for (int i = 0; i < this->vertices.size(); i++) {
    for (auto [neigh, edge] : this->neighbors[i]) {
      PRINTF("(src: %d dst: %d, cap: %d) ", i, neigh, edge.cap);
    }
    PRINTF("\n");
  }
}

void Graph::printEdgesVisualized() {
  PRINTF("digraph G {\n");
  for (int i = 0; i < this->vertices.size(); i++) {
    for (auto [neigh, edge] : this->neighbors[i]) {
      if (edge.initial_cap > 0)
        PRINTF("%d -> %d [label = \"%d/%d\"];\n", i, neigh,
               std::max(0, edge.initial_cap - edge.cap), edge.initial_cap);
    }
  }
  PRINTF("}\n");
}

void Graph::validate() {
  std::vector<int> vertFlows(this->vertices.size());
  for (int vert = 0; vert < this->vertices.size(); vert++) {
    // PRINTF("outer vert %d\n", vert);
    for (auto [i, edge] : this->neighbors[vert]) {
      vertFlows[vert] += std::max(0, edge.initial_cap - edge.cap);
      vertFlows[i] -= std::max(0, edge.initial_cap - edge.cap);
    }
  }
  for (int vert = 0; vert < this->vertices.size(); vert++) {
    if (vert == 0 || vert == 1)
      assert(vertFlows[0] == -vertFlows[1]);
    else
      assert(vertFlows[vert] == 0);
  }
  // PRINTF("validation success!\n");
  // return false;
}

inline bool Graph::isLayerReachable(const Vertex &srcVert, const Vertex &dstVert) {
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

bool Graph::visitVertex(Vertex &srcVert, Vertex &dstVert) {
  srcVert.layered_dst.push_back(dstVert.index);
  PRINTF("is layer set? %d\n", dstVert.layer == UNSET);
  if (dstVert.layer == UNSET) {
    // frontier.push(dstVert.index);
    dstVert.layer = srcVert.layer + 1;
    return true;
  } else {
    if (dstVert.layer != srcVert.layer + 1) {
      PRINTF("index: %d %d\n", dstVert.index, dstVert.layer);
      abort();
    }
  }
  return false;
}

bool Graph::bfs() {
  auto start = std::chrono::steady_clock::now();
  bool foundSink = false;
  std::queue<int> frontier;
  std::vector<bool> visited(this->vertices.size(), false);
  this->vertices[SOURCE].layer = 0;
  frontier.push(SOURCE);
  while (!frontier.empty()) {
    int src = frontier.front();
    if (src == SINK)
      foundSink = true;
    assert(!visited[src]);
    frontier.pop();
    Vertex &srcVert = this->vertices[src];
    for (auto &[dst, edge] : this->neighbors[src]) {
      Vertex &dstVert = this->vertices[dst];
      if (!visited[dstVert.index] && isLayerReachable(srcVert, dstVert) &&
          visitVertex(srcVert, dstVert)) {
        frontier.push(dst);
      }
    }
    visited[srcVert.index] = true;
  }

  auto end = std::chrono::steady_clock::now();
  bfs_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();
  return foundSink;
}
void Graph::increment(int node) {
  if (node >= 0)
    this->vertices[node].current_edge++;
}
bool Graph::dfsDeadEdge() {
  auto start = std::chrono::steady_clock::now();
  std::vector<bool> visited(this->vertices.size(), false);
  std::stack<int> stack;
  stack.push(SOURCE);
  while (!stack.empty()) {
    int nodeInd = stack.top();
    Vertex &srcVert = this->vertices[nodeInd];
    visited[nodeInd] = true;

    int s = this->edgesStart[nodeInd];
    int c = this->layeredEdgesCount[nodeInd];

    if (srcVert.current_edge == c) {
      stack.pop();
      increment(srcVert.parent);
      continue;
    }
    assert(srcVert.current_edge < c);
    int neigh = layeredEdges[s + srcVert.current_edge];
    Edge e = this->neighbors[nodeInd][neigh];

    // Vertex &dstVert = this->vertices[neigh];
    if (visited[neigh] || this->edgeCapacities[e.index] == 0) {
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

void Graph::reset() {
  for (Vertex &v : this->vertices) {
    v.reset();
  }
}

std::ostream &operator<<(std::ostream &os, const Graph &graph) {
  for (int i = 0; i < graph.neighbors.size(); i++) {
    os << "(" << i << " " << graph.vertices[i].layer << ") ";
  }
  os << "\n";
  return os;
}

bool Graph::bfsCudaTimed() {
  auto start = std::chrono::steady_clock::now();
  bool result = bfsCuda();
  auto end = std::chrono::steady_clock::now();
  bfs_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
              .count();
  return result;
}

void Graph::dinicsAlgo() {
  while (bfsCudaTimed()) {
    // PRINTF("Did a BFS\n");

    // printEdges();
    int i = 0;
    while (dfsDeadEdge()) {
      i++;
      // PRINTF("Finished dfs iteration\n");
      Vertex &dstVert = vertices[SINK];
      int minCapacity = UNSET;
      for (Vertex cur = dstVert; cur != vertices[SOURCE];
           cur = vertices[cur.parent]) {
        minCapacity =
            std::min(minCapacity, edgeCapacities[neighbors[cur.parent][cur.index].index]);
      }
      for (Vertex cur = dstVert; cur != vertices[SOURCE];
           cur = vertices[cur.parent]) {
        edgeCapacities[neighbors[cur.parent][cur.index].index] -= minCapacity;
        edgeCapacities[neighbors[cur.index][cur.parent].index] += minCapacity;
      }
    }
    printf("dfs count = %d\n", i);

    reset();
  }
}
