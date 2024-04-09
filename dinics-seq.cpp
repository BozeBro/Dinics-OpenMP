
// create aug graph

// find a layered graph (BFS)

// run lots of DFS to find blocking flow in layered graph

/*
Be able to change flow of Forward and Backward Edge Easily
Be able to run BFS and DFS through it in a simple way

*/
#include <unordered_map>
#include <vector>

struct Edge {
  int cap;
};

struct Vertex {
  int index;
  int current_edge = 0;
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

int main(int argc, char const *argv[])
{
  // Generate graph from input
  std::vector< std::unordered_map< Vertex, Edge > > graph;
  // Make the layer graph BFS
  /* 
  1. Run BFS on the graph
    a. notate level of each vertex in the graph.  
    b. construct lists for "next vertex" for Blocking Flow
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

