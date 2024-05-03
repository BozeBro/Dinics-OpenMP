// #include "dinics-seq.hpp"
// create aug graph

// find a layered graph (BFS)

// run lots of DFS to find blocking flow in layered graph

/*
Be able to change flow of Forward and Backward Edge Easily
Be able to run BFS and DFS through it in a simple way

*/
#include <cstdlib>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <ostream>
#include <queue>
#include <stack>
#include <string>
#include <unistd.h>
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

static double bfs_time = 0;
static double dfs_time = 0;

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

            assert(visited[dstVert.index]);
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
      assert(startIndex >= 0);
      assert(startIndex < this->vertices.size());

      if (localFrontier.size() > 0)
        PRINTF("Start localFrontier print %d\n", omp_get_thread_num());
      for (int k = 0; k < localFrontier.size(); k++) {
        PRINTF("%d ", startIndex + k);
        assert(k >= 0 && k < this->vertices.size());
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
    assert(!visited[src]);
    frontier.pop();
    Vertex &srcVert = this->vertices[src];
    for (auto &[dst, edge] : this->neighbors[src]) {
      Vertex &dstVert = this->vertices[dst];
      if (isLayerReachable(srcVert, dstVert) && visitVertex(srcVert, dstVert)) {
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
    assert(srcVert.current_edge < neighborEdges.size());
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

int main(int argc, char **argv) {

  // Generate graph from input
  // first two are src, dst
  // src: 0
  // dst: 1
  // nodes 2,3
  // src -5-> [2,3] -5-> dstVert
  const auto init_start = std::chrono::steady_clock::now();
  std::string input_filename;
  char mode = '\0';
  int num_threads = 0;
  int opt;
  while ((opt = getopt(argc, argv, "f:n:m:")) != -1) {
    switch (opt) {
    case 'f':
      input_filename = optarg;
      break;
    case 'n':
      num_threads = atoi(optarg);
      break;
    case 'm':
      mode = *optarg;
      break;

    default:
      std::cerr << "Usage: " << argv[0]
                << " -f input_filename -n num_threads -m parallel_mode -b "
                   "batch_size\n";
      exit(EXIT_FAILURE);
    }
  }
  if (empty(input_filename)) {
    std::cerr << "Usage: " << argv[0]
              << " -f input_filename -n num_threads "
                 "-m parallel_mode\n";
    exit(EXIT_FAILURE);
  }

  std::ifstream fin(input_filename);
  if (!fin) {
    std::cerr << "Unable to open file: " << input_filename << ".\n";
    exit(EXIT_FAILURE);
  }
  int n;
  fin >> n;
  Graph graph(n);
  for (int i = 0; i < n; i++) {
    int cnt;
    fin >> cnt;
    if (cnt == 0)
      continue;
    graph.neighbors[i].reserve(cnt);
    int neigh, cap;
    for (int j = 0; j < cnt; j++) {
      fin >> neigh >> cap;

      graph.addEdge({i}, {neigh}, cap);
    }
  }

  std::unordered_map<Vertex, Edge> src_map;

  /*
    Make the layer graph BFS
    1. Run BFS on the graph
      a. notate level of each vertex in the graph.
      b. construct lists for "next vertex" for Blocking Flow
  */
  const double init_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - init_start)
          .count();
  const auto compute_start = std::chrono::steady_clock::now();
  // graph.printEdgesVisualized();
  // exit(0);
  while (graph.bfs()) {
    // PRINTF("Did a BFS\n");

    // graph.printEdges();
    while (graph.dfsDeadEdge()) {
      // PRINTF("Finished dfs iteration\n");
      Vertex &dstVert = graph.vertices[SINK];
      int minCapacity = UNSET;
      for (Vertex cur = dstVert; cur != graph.vertices[SOURCE];
           cur = graph.vertices[cur.parent]) {
        minCapacity =
            std::min(minCapacity, graph.neighbors[cur.parent][cur.index].cap);
      }
      for (Vertex cur = dstVert; cur != graph.vertices[SOURCE];
           cur = graph.vertices[cur.parent]) {
        graph.neighbors[cur.parent][cur.index].cap -= minCapacity;
        graph.neighbors[cur.index][cur.parent].cap += minCapacity;
      }
    }
    graph.reset();
  }
  int flow = 0;

  for (auto [i, edge] : graph.neighbors[SOURCE]) {
    // std::cout << i << ": " << edge.initial_cap << " " << edge.cap << '\n';
    flow += edge.initial_cap - edge.cap;
  }
  for (int src = 1; src < n; src++) {
    auto edge = graph.neighbors[src][SOURCE];
    flow -= std::max(0, edge.initial_cap - edge.cap);
  }
  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  printf("Compute time %f\n", compute_time);
  printf("Init time %f\n", init_time);
  printf("Total time %f\n", compute_time + init_time);
  printf("BFS time %f\n", bfs_time);
  printf("DFS time %f\n", dfs_time);
  printf("Flow value %d\n", flow);

  return 0;
}
