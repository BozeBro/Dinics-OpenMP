
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
#include <omp.h>

#define UNSET std::numeric_limits<int>::max()
#define SOURCE 0
#define SINK 1
#define MY_PRINT_ENABLED
#ifdef MY_PRINT_ENABLED
// If it's defined, define PRINTF as printf
#define PRINTF(format, ...) printf(format, ##__VA_ARGS__)
#else
// If it's not defined, define PRINTF as an empty macro
#define PRINTF(format, ...) ((void)0)
#endif

static double bfs_time = 0;
static double dfs_time = 0;

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
struct Edge {
  int cap;
  int initial_cap;
};
namespace std {
template <> struct hash<Vertex> {
  std::size_t operator()(const Vertex &obj) const {
    // index is the unique identifier anyway
    return std::hash<int>()(obj.index);
  }
};
} // namespace std
bool operator==(const Vertex &a, const Vertex &b) {
  return a.index == b.index && a.current_edge == b.current_edge &&
         a.layer == b.layer;
}
bool operator!=(const Vertex &a, const Vertex &b) { return !(a == b); }

struct Graph {

  std::vector<Vertex> vertices;
  std::vector<std::unordered_map<int, Edge>> neighbors;

  Graph(int size) : vertices(size), neighbors(size) {
    for (int i = 0; i < size; i++) {
      vertices[i].index = i;
    }
  }

  void printEdges() {
    for (int i = 0; i < this->vertices.size(); i++) {
      for (auto [neigh, edge] : this->neighbors[i]) {
        printf("(src: %d dst: %d, cap: %d) ", i, neigh, edge.cap);
      }
      printf("\n");
    }
  }

  void printEdgesVisualized() {
    printf("digraph G {\n");
    for (int i = 0; i < this->vertices.size(); i++) {
      for (auto [neigh, edge] : this->neighbors[i]) {
        if (edge.initial_cap > 0)
          printf("%d -> %d [label = \"%d/%d\"];\n", i, neigh,
                 std::max(0, edge.initial_cap - edge.cap), edge.initial_cap);
      }
    }
    printf("}\n");
  }

  void validate() {
    std::vector<int> vertFlows(this->vertices.size());
    for (int vert = 0; vert < this->vertices.size(); vert++) {
      // printf("outer vert %d\n", vert);
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
    // printf("validation success!\n");
    // return false;
  }
  bool bfsParallel() {
    // initialize visited array
    std::vector<bool> visited(this->vertices.size(), false);
    // initialize current frontier
    std::vector<int> frontier{SOURCE};
    std::vector<int> newFrontier(this->vertices.size());
    int size;
    int lastWrite;
    // initialize threads
    // Per iteration:
    while (!frontier.empty()) {
      // newFrontier.reserve(this->vertices.size());
      size = 0;
      lastWrite = 0;
      bool found = false;
      #pragma omp parallel shared(newFrontier, frontier, lastWrite, found) threadprivate(localFrontier, localLastWrite) {
        // initialize local frontier
        std::vector<int> localFrontier;
      
        #pragma omp parallel for reduction(+: size)
        for (int i = 0; i < frontier.size(); i++) {
          int index = frontier[i];
          for (auto [neigh, edge] : this->neighbors[index]) {
            bool edited = true;
            #pragma omp atomic capture 
            { edited = visited[neigh]; visited[neigh] = true; }
            if (!edited) {
              if (index == SINK)
                found = true;
              Vertex& srcVert = this->neighbors[index];
              Vertex& dstVert = this->vertices[neigh];
              visitVertex(srcVert, dstVert);
              localFrontier.push_back(neigh);
              size += 1;
            }
          }
        }
        int startIndex = -1;
        #pragma omp atomic capture 
        { startIndex = lastWrite; lastWrite += localFrontier.size(); }
        memcpy(newFrontier.data() + startIndex * sizeof(int), localFrontier.data(), sizeof(int) * localFrontier.size());
      }
      if (found)
        return true;
      frontier = newFrontier;
    }
  
    // Per thread:
    // Pick out vert from old frontier
    // For each neighbor:
    // Test and set visited 
    // If visited already, skip
    // Else, add it to the new frontier and update accordingly
      // change depth
      // set visited
      // add to local frontier
    // upon finish: 
    /*
    #pragma parallel threadprivate(tfront, last_write) {
    #program omp for  {
        // do work --> local tfront edited, last_write edited
    }
    atomic {
      i = glast_write; i += last_write;
    }
    mempcy(frontier, )
    }

    */

    // Merge all thread new frontiers, set old frontier to that
    // Exit if we found the sink
    return false;
  }
  
  void visitVertex(Vertex& srcVert, Vertex& dstVert) {
    if (dstVert.layer <= srcVert.layer || visited[dstVert.index] ||
        this->neighbors[srcVert.index][dst].cap == 0)
      continue;
    srcVert.layered_dst.push_back(dst);
    if (dstVert.layer == UNSET) {
      frontier.push(dstVert.index);
      dstVert.layer = srcVert.layer + 1;
    } else {
      if (dstVert.layer != srcVert.layer + 1) {
        printf("index: %d %d\n", dstVert.index, dstVert.layer);
        abort();
      }
    }
  }

  bool bfs() {
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
        visitVertex(srcVert, dstVert);
      }
      visited[srcVert.index] = true;
    }

    auto end = std::chrono::steady_clock::now();
    bfs_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    return foundSink;
  }
  void increment(int node) {
    if (node >= 0)
      this->vertices[node].current_edge++;
  }
  bool dfsDeadEdge() {
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
        dfs_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
        return true;
      }
      stack.push(neigh);
    }
    auto end = std::chrono::steady_clock::now();
    dfs_time += std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    return false;
  }

  bool dfs() {
    std::vector<bool> visited(this->vertices.size(), false);
    std::stack<int> stack;
    stack.push(SOURCE);
    while (!stack.empty()) {
      int nodeInd = stack.top();
      // printf("%d\n", nodeInd);
      stack.pop();
      if (visited[nodeInd])
        continue;
      visited[nodeInd] = true;
      for (auto neigh : this->vertices[nodeInd].layered_dst) {
        assert(this->neighbors[nodeInd][neigh].cap >= 0);
        if (visited[neigh] || this->neighbors[nodeInd][neigh].cap == 0)
          continue;
        // printf("src: %d, dst: %d, cap: %d\n", nodeInd, neigh,
        //        this->neighbors[nodeInd][neigh].cap);
        this->vertices[neigh].parent = nodeInd;
        stack.push(neigh);
        if (neigh == SINK)
          return true;
      }
    }
    return false;
  }

  void addEdge(const Vertex &start, const Vertex &end, int cap) {
    // Max prevent capacity override if there are 2 node cycles
    int cap_value = std::max(neighbors[start.index][end.index].cap, cap);
    neighbors[start.index][end.index] = {cap_value, cap_value};
    cap_value = std::max(neighbors[end.index][start.index].cap, 0);
    neighbors[end.index][start.index] = {cap_value, cap_value};
    // neighbors[end.index][start.index].cap =
    //     std::max(neighbors[end.index][start.index].cap, 0);
  }

  void reset() {
    for (Vertex &v : this->vertices) {
      v.reset();
    }
  }

  friend std::ostream &operator<<(std::ostream &os, const Graph &graph) {
    for (int i = 0; i < graph.neighbors.size(); i++) {
      os << "(" << i << " " << graph.vertices[i].layer << ") ";
    }
    os << "\n";
    return os;
  }
};

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
  // int n = 5;
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
  // graph.addEdge({0}, {2}, 5);
  // graph.addEdge({0}, {3}, 3);

  // graph.addEdge({2}, {3}, 10);
  // graph.addEdge({2}, {4}, 3);

  // graph.addEdge({3}, {4}, 10);

  // graph.addEdge({4}, {1}, 10);
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
  while (graph.bfs()) {
    // printf("Did a BFS\n");
    // graph.printEdges();
    while (graph.dfsDeadEdge()) {
      // printf("Finished dfs iteration\n");
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
  // graph.printEdges();
  // printf("Start validating\n");
  // graph.validate();
  // printf("End validating\n");
  int flow = 0;

  for (auto [i, edge] : graph.neighbors[SOURCE]) {
    // std::cout << i << ": " << edge.initial_cap << " " << edge.cap << '\n';
    flow += edge.initial_cap - edge.cap;
  }
  // std::cout << "\n\n";
  for (int src = 1; src < n; src++) {
    auto edge = graph.neighbors[src][SOURCE];
    // std::cout << src << ": " << edge.initial_cap << " " << edge.cap << '\n';
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

  // graph.printEdgesVisualized();

  /*

    s -> a

  */
  /*
  2. Run DFS on the level graph
    a. have pointer/iterator to edge for each vertex, increment on dead edge
    b. loop and update aug graph accordingly
      - dead edges happen when DFS does not reach destination
      - when going back, increment pointer for edge and try again
      - when find a path to destination, subtract the minimum flow of any edge
  on the path
  */

  // Get blocking flow (recursive)
  return 0;
}
