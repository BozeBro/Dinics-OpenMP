#include "dinics-seq.hpp"
#include "mgraph.hpp"
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <omp.h>
#include <ostream>
#include <unistd.h>

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
  omp_set_num_threads(num_threads);

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
    // graph.adj_list[i].reserve(cnt);
    int neigh, cap;
    for (int j = 0; j < cnt; j++) {
      fin >> neigh >> cap;

      Vertex a{i};
      Vertex b{neigh};
      graph.addEdge(a, b, cap);
    }
  }
  // graph.printEdgesVisualized();
  // exit(0);

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
  graph.dinicsAlgo();

  int flow = 0;

  for (auto [_, edge] : graph.neighbors[SOURCE]) {
    // std::cout << i << ": " << edge.initial_cap << " " << edge.cap << '\n';
    flow += edge.initial_cap - edge.cap;
    // flow += i - edge;
  }
  for (int src = 1; src < n; src++) {
    Edge edge = graph.neighbors[src][SOURCE];
    flow -= std::max(0, edge.initial_cap - edge.cap);
  }
  const double compute_time =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          std::chrono::steady_clock::now() - compute_start)
          .count();
  printf("Compute time %f\n", compute_time);
  printf("Init time %f\n", init_time);
  printf("Total time %f\n", compute_time + init_time);
  printf("BFS time %f\n", graph.bfs_time);
  printf("DFS time %f\n", graph.dfs_time);
  printf("Flow value %d\n", flow);

  return 0;
}
