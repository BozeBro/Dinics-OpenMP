#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "dinics-seq.hpp"
#include <chrono>

__global__ void bfsKernel(char* frontier, char* newFrontier, int level, 
    int* edgesStart, int* edgesCount, int* edges, int* layeredEdgesCount, 
    int* layeredEdges, int* edgeCapacities, unsigned int* vertDists, 
    int count, char* foundSink, char* progressed) {
  int vert = blockIdx.x * blockDim.x + threadIdx.x;

  if (!frontier[vert] || vert >= count)
    return;

  int es = edgesStart[vert];
  int ec = edgesCount[vert];
  int dist = level;
  int lc = 0;

  for (int i = 0; i < ec; i++) {
    int dest = edges[es + i];
    
    if (edgeCapacities[es + i] > 0 && vertDists[dest] >= dist + 1) {
      if (dest == 1)
        *foundSink = true;

      vertDists[dest] = dist + 1;
      newFrontier[dest] = true;
      layeredEdges[es + lc] = dest;
      *progressed = true;
      lc++;
    }
  }

  layeredEdgesCount[vert] = lc;
}

void Graph::initCuda() {
  cudaMalloc(&cudaFrontier, sizeof(char) * vertices.size());
  cudaMalloc(&cudaNewFrontier, sizeof(char) * vertices.size());
  cudaMalloc(&cudaEdgesStart, sizeof(int) * vertices.size());
  cudaMalloc(&cudaEdgesCount, sizeof(int) * vertices.size());
  cudaMalloc(&cudaEdges, sizeof(int) * edgeCount);
  cudaMalloc(&cudaLayeredEdgesCount, sizeof(int) * vertices.size());
  cudaMalloc(&cudaLayeredEdges, sizeof(int) * edgeCount);
  cudaMalloc(&cudaEdgeCapacities, sizeof(int) * edgeCount);
  cudaMalloc(&cudaVertDists, sizeof(int) * vertices.size());
  cudaMalloc(&cudaFoundSink, sizeof(char));
  cudaMalloc(&cudaProgressed, sizeof(char));

  cudaMemcpy(cudaEdges, edges.data(), edgeCount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(cudaEdgesStart, edgesStart.data(), vertices.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(cudaEdgesCount, edgesCount.data(), vertices.size() * sizeof(int), cudaMemcpyHostToDevice);
}

void Graph::destroyCuda() {
  cudaFree(cudaFrontier);
  cudaFree(cudaNewFrontier);
  cudaFree(cudaEdgesStart);
  cudaFree(cudaEdgesCount);
  cudaFree(cudaEdges);
  cudaFree(cudaLayeredEdgesCount);
  cudaFree(cudaLayeredEdges);
  cudaFree(cudaEdgeCapacities);
  cudaFree(cudaVertDists);
  cudaFree(cudaFoundSink);
}

bool Graph::bfsCuda() {
  cudaMemset(cudaFrontier, 0, sizeof(char) * vertices.size());
  cudaMemset(cudaFrontier, 1, sizeof(char));
  cudaMemset(cudaNewFrontier, 0, sizeof(char) * vertices.size());
  cudaMemset(cudaVertDists, 255, sizeof(int) * vertices.size());
  cudaMemset(cudaVertDists, 0, sizeof(int));
  cudaMemset(cudaFoundSink, 0, sizeof(char));

  // std::vector<int> edges(edgeCount);
  // std::vector<int> edgeCapacities(edgeCount);
  // std::vector<int> edgesStart(vertices.size());
  // std::vector<int> edgesCount(vertices.size());

  auto start = std::chrono::steady_clock::now();
  // int ec = 0;
  // for (int i = 0; i < vertices.size(); i++) {
  //   edgesStart[i] = ec;
  //   int numNeighbors = neighbors[i].size();
  //   edgesCount[i] = numNeighbors;

  //   int e = 0;
  //   for (auto &[dst, edge] : this->neighbors[i]) {
  //     edgeCapacities[ec + e] = edge.cap;
  //     edges[ec + e] = dst;
  //     e++;
  //   }
  //   ec += e;
  // }

  int level = 0;
  cudaMemcpy(cudaEdgeCapacities, edgeCapacities.data(), edgeCount * sizeof(int), cudaMemcpyHostToDevice);
  auto end = std::chrono::steady_clock::now();

  const int threadsPerBlock = 512;
  const int blocks = (vertices.size() + threadsPerBlock - 1) / threadsPerBlock;

  char foundSink = false;
  char progressed = false;
  bfs_aux_time +=
      std::chrono::duration_cast<std::chrono::duration<double>>(end - start)
          .count();

  do {
    cudaMemset(cudaProgressed, 0, sizeof(char));
    bfsKernel<<<blocks, threadsPerBlock>>>(cudaFrontier, cudaNewFrontier, level, 
      cudaEdgesStart, cudaEdgesCount, cudaEdges, cudaLayeredEdgesCount, 
      cudaLayeredEdges, cudaEdgeCapacities, cudaVertDists, 
      vertices.size(), cudaFoundSink, cudaProgressed);

    cudaMemcpy(cudaFrontier, cudaNewFrontier, sizeof(char) * vertices.size(), cudaMemcpyDeviceToDevice);
    cudaMemset(cudaNewFrontier, 0, sizeof(char) * vertices.size());

    cudaMemcpy(&foundSink, cudaFoundSink, sizeof(char), cudaMemcpyDeviceToHost);
    cudaMemcpy(&progressed, cudaProgressed, sizeof(char), cudaMemcpyDeviceToHost);

    level++;
  } while (!foundSink && progressed);

  cudaMemcpy(layeredEdges.data(), cudaLayeredEdges, edgeCount * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(layeredEdgesCount.data(), cudaLayeredEdgesCount, vertices.size() * sizeof(int), cudaMemcpyDeviceToHost);
  // for (int i = 0; i < vertices.size(); i++) {
  //   vertices[i].layered_dst.resize(layeredEdgesCount[i]);
  //   for (int j = 0; j < layeredEdgesCount[i]; j++) {
  //     vertices[i].layered_dst[j] = layeredEdges[edgesStart[i] + j];
  //   }
  // }

  return foundSink;
}

