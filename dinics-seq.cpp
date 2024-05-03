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

#ifdef USE_OPEN_MP
bool Graph::bfsParallel() {
  auto start = std::chrono::steady_clock::now();
  // initialize visited array
  std::vector<int> visited(this->vertices.size(), 0);
  // initialize current frontier
  std::vector<int> frontier{SOURCE};
  this->vertices[SOURCE].layer = 0;
  int size;
  bool found = false;
  while (!frontier.empty() && !found) {
    size = 0;  
    std::vector<int> threadFrontiers[MAX_THREADS];
    #pragma omp parallel for
    for (int i = 0; i < frontier.size(); i++) {
      int threadIndex = omp_get_thread_num();
      std::vector<int>& localFrontier = threadFrontiers[threadIndex];
      int index = frontier[i];
      // if (!(0 <= index && index < this->vertices.size())) {
      //   PRINTF("index is %d\n", index);
      //   abort();
      // }
      Vertex &srcVert = this->vertices[index];
      for (const auto [neigh, edge] : this->neighbors[index]) {
        Vertex &dstVert = this->vertices[neigh];
        if (dstVert.layer <= srcVert.layer || edge.cap == 0)
          continue;
        if (neigh == SINK) {
          PRINTF("done\n");
          found = true;
        }

//             assert(visited[dstVert.index]);
        if (visitVertexParallel(srcVert, dstVert)) {
//               PRINTF("adding vertex %d\n", dstVert.index);
          localFrontier.push_back(neigh);
        }
      }
    }

    int startIndices[MAX_THREADS];
    startIndices[0] = 0;
    int size = threadFrontiers[0].size();
    for (int i = 1; i < MAX_THREADS; i++) {
      startIndices[i] = size;
      size += threadFrontiers[i].size();
    }

    frontier.resize(size);
    // The time this takes is negligible
    #pragma omp parallel
    {
      int num = omp_get_thread_num();
      int startIndex = startIndices[num];
      std::vector<int>& localFrontier = threadFrontiers[num];
      // assert(startIndex >= 0);
      // assert(startIndex < this->vertices.size());

      // if (localFrontier.size() > 0)
      //   PRINTF("Start localFrontier print %d\n", omp_get_thread_num());
      for (int k = 0; k < localFrontier.size(); k++) {
        // PRINTF("%d ", startIndex + k);
        // assert(k >= 0 && k < this->vertices.size());
        frontier[startIndex + k] = localFrontier[k];
      }
    //   if (localFrontier.size() > 0)
    //     PRINTF("\n End local Frontier print\n");
    //   memcpy(newFrontier + startIndex * sizeof(int),
    //          localFrontier.data(), sizeof(int) * localFrontier.size());
    }
    // PRINTF("newFrontier %d\n", lastWrite);
    // for (int i = 0; i < this->vertices.size(); i++) {
    //   PRINTF("%d ", newFrontier[i]);
    // }
    // PRINTF("\n");
    // PRINTF("Size is %d\n", lastWrite);
  }
  auto end = std::chrono::steady_clock::now();
  bfs_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();
  return found;
}
#endif

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

#ifdef USE_OPEN_MP
inline bool Graph::visitVertexParallel(Vertex &srcVert, Vertex &dstVert) {
  int num = omp_get_thread_num();
  srcVert.layered_dst_array[num].push_back(dstVert.index);
  PRINTF("is layer set? %d\n", dstVert.layer == UNSET);
  // maybe make this atomic?
  // on second thought, maybe don't... dfs already checks for duplicates
  if (dstVert.layer == UNSET) {
    // frontier.push(dstVert.index);
    dstVert.layer = srcVert.layer + 1;
    return true;
  } 
  // else {
  //   if (dstVert.layer != srcVert.layer + 1) {
  //     PRINTF("index: %d %d\n", dstVert.index, dstVert.layer);
  //     abort();
  //   }
  // }
  return false;
}
#endif

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
    auto &neighborEdges = srcVert.layered_dst_array[srcVert.current_edge_array];
    if (srcVert.current_edge == neighborEdges.size()) {
      if (srcVert.current_edge_array == MAX_THREADS - 1) {
        stack.pop();
        increment(srcVert.parent);
      } else {
        srcVert.current_edge = 0;
        srcVert.current_edge_array++;
      }
      continue;
    }
    assert(srcVert.current_edge < neighborEdges.size());
    int neigh = neighborEdges[srcVert.current_edge];
    // Vertex &dstVert = this->vertices[neigh];
    if (visited[neigh] || this->neighbors[nodeInd][neigh].cap == 0) {
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

bool Graph::dfs() {
  std::vector<bool> visited(this->vertices.size(), false);
  std::stack<int> stack;
  stack.push(SOURCE);
  while (!stack.empty()) {
    int nodeInd = stack.top();
    stack.pop();
    if (visited[nodeInd])
      continue;
    visited[nodeInd] = true;
    for (auto neigh : this->vertices[nodeInd].layered_dst) {
      assert(this->neighbors[nodeInd][neigh].cap >= 0);
      if (visited[neigh] || this->neighbors[nodeInd][neigh].cap == 0)
        continue;
      this->vertices[neigh].parent = nodeInd;
      stack.push(neigh);
      if (neigh == SINK)
        return true;
    }
  }
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

void Graph::dinicsAlgo() {
  while (bfsParallel()) {
    // PRINTF("Did a BFS\n");

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
