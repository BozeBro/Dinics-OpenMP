
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
          printf("%d -> %d [label = \"%d/%d\"];\n", i, neigh, std::max(0, edge.initial_cap - edge.cap), edge.initial_cap);
      }
    }
    printf("}\n");
  }

  void validate() {
    int srcFlow = -1;
    for (int vert = 0; vert < this->vertices.size(); vert++) {
      int flow = 0;
      for (auto [i, edge] : this->neighbors[vert]) {
        flow += std::max(0, edge.initial_cap - edge.cap);
      }
      for (int src = 0; src < this->vertices.size(); src++) {
        auto edge = this->neighbors[src][vert];
        flow -= std::max(0, edge.initial_cap - edge.cap);
      }
      if (vert == 0)
        srcFlow = flow;
      else if (vert == 1)
        assert(flow == -srcFlow);
      else
        assert(flow == 0);
    }
    printf("validation success!\n");
  }

  bool bfs() {
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
      visited[srcVert.index] = true;
    }
    return foundSink;
  }
  void increment(int node) {
    if (node >= 0)
      this->vertices[node].current_edge++;
  }
  bool dfsDeadEdge() {
    std::vector<bool> visited(this->vertices.size(), false);
    std::stack<int> stack;
    stack.push(SOURCE);
    while (!stack.empty()) {
      int nodeInd = stack.top();
      Vertex &srcVert = this->vertices[nodeInd];
      if (srcVert.parent != -1)
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
      if (neigh == SINK)
        return true;
      stack.push(neigh);
    }
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
  while (graph.bfs()) {
    printf("Did a BFS\n");
    // graph.printEdges();
    while (graph.dfsDeadEdge()) {
      printf("Finished dfs iteration\n");
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
  int flow = 0;
  for (auto [i, edge] : graph.neighbors[SOURCE]) {
    std::cout << i << ": " << edge.initial_cap << " " << edge.cap << '\n';
    flow += edge.initial_cap - edge.cap;
  }
<<<<<<< HEAD
  printf("%d\n", flow);
=======
  std::cout << "\n\n";
  for (int src = 1; src < n; src++) {
    auto edge = graph.neighbors[src][SOURCE];
    std::cout << src << ": " << edge.initial_cap << " " << edge.cap << '\n';
    flow -= std::max(0, edge.initial_cap - edge.cap);
  }
  graph.printEdgesVisualized();
  std::cout << flow << std::endl;
  graph.validate();
>>>>>>> refs/remotes/origin/main
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
