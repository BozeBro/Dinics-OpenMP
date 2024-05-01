#include "dinics-seq.hpp"
// create aug graph

// find a layered graph (BFS)

// run lots of DFS to find blocking flow in layered graph

/*
Be able to change flow of Forward and Backward Edge Easily
Be able to run BFS and DFS through it in a simple way

*/
#include <algorithm>
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
      ASSERT(vertFlows[0] == -vertFlows[1]);
    else
      ASSERT(vertFlows[vert] == 0);
  }
  // PRINTF("validation success!\n");
  // return false;
}
bool Graph::bfsParallel() {
  auto start = std::chrono::steady_clock::now();
  // initialize visited array
  std::vector<int> visited(this->vertices.size(), 0);
  // initialize current frontier
  std::vector<int> frontier{SOURCE};
  this->vertices[SOURCE].layer = 0;
  int *newFrontier = new int[this->vertices.size()];
  for (int i = 0; i < this->vertices.size(); i++) {
    newFrontier[i] = -1;
  }
  int size;
  int lastWrite = 0;
  // // initialize threads
  // // Per iteration:
  // PRINTF("size: %d\n", this->neighbors[SOURCE].size());
  // for (auto [ind, dst] : this->neighbors[SOURCE]) {
  //   PRINTF("%d ", ind);
  // }
  // PRINTF("\n");
  bool found = false;
  while (!frontier.empty() && !found) {
    size = 0;
    lastWrite = 0;

    std::vector<int> localFrontier;
    if (frontier.size() != 10)
      PRINTF("%zu\n", frontier.size());
#pragma omp parallel shared(newFrontier, frontier, lastWrite, found, visited)  \
    firstprivate(localFrontier)
    {

#pragma omp parallel for
      for (int i = 0; i < frontier.size(); i++) {
        int index = frontier[i];
        if (!(0 <= index && index < this->vertices.size())) {
          PRINTF("index is %d\n", index);
          abort();
        }
        for (const auto [neigh, edge] : this->neighbors[index]) {
          Vertex &srcVert = this->vertices[index];
          Vertex &dstVert = this->vertices[neigh];
          if (!isLayerReachable(srcVert, dstVert)) {
            continue;
          }
          bool edited = true;
#pragma omp atomic capture
          {
            edited = visited[neigh];
            visited[neigh] = true;
          }
          if (!edited) {
            if (neigh == SINK) {
              PRINTF("done\n");
              found = true;
            }

            ASSERT(visited[dstVert.index]);
            if (visitVertexParallel(srcVert, dstVert)) {
              PRINTF("adding vertex %d\n", dstVert.index);
              localFrontier.push_back(neigh);
            }
          }
        }
      }
      int startIndex = -1;
#pragma omp atomic capture
      {
        startIndex = lastWrite;
        lastWrite += localFrontier.size();
      }
      ASSERT(startIndex >= 0);
      ASSERT(startIndex < this->vertices.size());

      if (localFrontier.size() > 0)
        PRINTF("Start localFrontier print %d\n", omp_get_thread_num());
      for (int k = 0; k < localFrontier.size(); k++) {
        PRINTF("%d ", startIndex + k);
        ASSERT(k >= 0 && k < this->vertices.size());
        newFrontier[startIndex + k] = localFrontier[k];
      }
      if (localFrontier.size() > 0)
        PRINTF("\n End local Frontier print\n");
      // memcpy(newFrontier + startIndex * sizeof(int),
      //        localFrontier.data(), sizeof(int) * localFrontier.size());
    }
    PRINTF("newFrontier %d\n", lastWrite);
    for (int i = 0; i < this->vertices.size(); i++) {
      PRINTF("%d ", newFrontier[i]);
    }
    PRINTF("\n");
    PRINTF("Size is %d\n", lastWrite);
    frontier.resize(lastWrite);
    for (int i = 0; i < lastWrite; i++) {
      frontier[i] = newFrontier[i];
    }
  }
  delete[] newFrontier;
  auto end = std::chrono::steady_clock::now();
  bfs_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();
  return found;
}

bool Graph::isLayerReachable(const Vertex &srcVert, const Vertex &dstVert) {
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
bool Graph::visitVertexParallel(Vertex &srcVert, Vertex &dstVert) {
#pragma omp critical
  { srcVert.layered_dst.push_back(dstVert.index); }

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
    ASSERT(!visited[src]);
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
    auto &neighborEdges = this->vertices[nodeInd].layered_dst;
    if (srcVert.current_edge == neighborEdges.size()) {
      stack.pop();
      increment(srcVert.parent);
      continue;
    }
    ASSERT(srcVert.current_edge < neighborEdges.size());
    int neigh = this->vertices[nodeInd].layered_dst[srcVert.current_edge];
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
      ASSERT(this->neighbors[nodeInd][neigh].cap >= 0);
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

void Graph::addEdge(const Vertex &start, const Vertex &end, int cap) {
  // Max prevent capacity override if there are 2 node cycles
  int cap_value = std::max(neighbors[start.index][end.index].cap, cap);
  neighbors[start.index][end.index] = {cap_value, cap_value};
  cap_value = std::max(neighbors[end.index][start.index].cap, 0);
  neighbors[end.index][start.index] = {cap_value, cap_value};
  // neighbors[end.index][start.index].cap =
  //     std::max(neighbors[end.index][start.index].cap, 0);
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
  while (bfs()) {
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
