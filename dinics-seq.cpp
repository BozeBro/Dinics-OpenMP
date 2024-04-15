
// create aug graph

// find a layered graph (BFS)

// run lots of DFS to find blocking flow in layered graph

/*
Be able to change flow of Forward and Backward Edge Easily
Be able to run BFS and DFS through it in a simple way

*/
#include <unordered_map>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <limits>
#include <cassert>
#include <iostream>
#include <ostream>
#include <cstdio>

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
    this->layered_dst.clear();
  }
};
struct Edge {
  int cap;
};
namespace std {
template <>
    struct hash<Vertex> {
        std::size_t operator()(const Vertex& obj) const {
            // Combine hashes of individual members to form the hash value
            return std::hash<int>()(obj.index);
        }
    };
}
bool operator==(const Vertex& a, const Vertex& b) {
  return a.index == b.index && a.current_edge == b.current_edge && a.layer == b.layer;
}
bool operator!=(const Vertex& a, const Vertex& b) {
  return !(a == b);
}

struct Graph {

  std::vector<Vertex> vertices;
  std::vector<std::unordered_map<int, Edge > > neighbors;

  Graph(int size) 
  : vertices(size)
  , neighbors(size) {
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
      for (auto& [dst, edge] : this->neighbors[src]) {
        Vertex &dstVert = this->vertices[dst];
        if (dstVert.layer <= srcVert.layer || visited[dstVert.index] || this->neighbors[srcVert.index][dst].cap == 0)
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
            
          // assert(dstVert.layer == srcVert.layer + 1);
        }
      }
      visited[srcVert.index] = true;
    }
    return foundSink;
  }

  bool dfs() {
    std::vector<bool> visited(this->vertices.size(), false);
    std::stack<int> stack;
    stack.push(SOURCE);
    while (!stack.empty()) {
      int nodeInd = stack.top();
      printf("%d\n", nodeInd);
      stack.pop();
      if (visited[nodeInd])
        continue;
      visited[nodeInd] = true;
      for (auto neigh : this->vertices[nodeInd].layered_dst) {
        assert(this->neighbors[nodeInd][neigh].cap >= 0);
        if (visited[neigh] || this->neighbors[nodeInd][neigh].cap == 0)
          continue;
        printf("src: %d, dst: %d, cap: %d\n", nodeInd, neigh, this->neighbors[nodeInd][neigh].cap);
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
    neighbors[start.index][end.index].cap = std::max(neighbors[start.index][end.index].cap, cap);
    neighbors[end.index][start.index].cap = std::max(neighbors[end.index][start.index].cap, 0);
  }

  void reset() {
    for (Vertex &v: this->vertices) {
      v.reset();
    }
  }

  friend std::ostream& operator<<(std::ostream& os, const Graph& graph) {
    for (int i = 0; i < graph.neighbors.size(); i++) {
      os << "(" << i << " " << graph.vertices[i].layer << ") ";
    }
    os << "\n";
    return os;
  }
};

int main(int argc, char const *argv[]) {
  
  // Generate graph from input
  // first two are src, dst
  // src: 0
  // dst: 1
  // nodes 2,3
  // src -5-> [2,3] -5-> dst
  int n = 5;
  Graph graph(n);
  graph.addEdge({0}, {2}, 5);
  graph.addEdge({0}, {3}, 3);

  graph.addEdge({2}, {3}, 10);
  graph.addEdge({2}, {4}, 3);

  graph.addEdge({3}, {4}, 10);

  graph.addEdge({4}, {1}, 10);
  std::unordered_map<Vertex, Edge> src_map;

  /*
    Make the layer graph BFS
    1. Run BFS on the graph
      a. notate level of each vertex in the graph.  
      b. construct lists for "next vertex" for Blocking Flow
  */
  while (graph.bfs()) {
    printf("Did a BFS\n");
    graph.printEdges();
    while (graph.dfs()) {
      printf("FInished dfs iteration\n");
      Vertex& dstVert = graph.vertices[SINK];
      int minCapacity = UNSET;
      for (Vertex cur = dstVert; cur != graph.vertices[SOURCE]; cur = graph.vertices[cur.parent]) {
        minCapacity = std::min(minCapacity, graph.neighbors[cur.parent][cur.index].cap);
      }
      for (Vertex cur = dstVert; cur != graph.vertices[SOURCE]; cur = graph.vertices[cur.parent]) {
        graph.neighbors[cur.parent][cur.index].cap -= minCapacity;
        graph.neighbors[cur.index][cur.parent].cap += minCapacity;
      }
    }
    graph.reset();
  }
  graph.printEdges();
  /*

    s -> a
      
  */
  /* 
  2. Run DFS on the level graph
    a. have pointer/iterator to edge for each vertex, increment on dead edge
    b. loop and update aug graph accordingly
      - dead edges happen when DFS does not reach destination
      - when going back, increment pointer for edge and try again
      - when find a path to destination, subtract the minimum flow of any edge on the path
  */
 

  // Get blocking flow (recursive)
  return 0;
}

