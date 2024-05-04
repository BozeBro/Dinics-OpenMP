#include <stdio.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#include "dinics-seq.hpp"

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

bool Graph::bfsCuda() {
  char* cudaFrontier;
  cudaMalloc(&cudaFrontier, sizeof(char) * vertices.size());
  cudaMemset(cudaFrontier, 0, sizeof(char) * vertices.size());
  cudaMemset(cudaFrontier, 1, sizeof(char));

  char* cudaNewFrontier;
  cudaMalloc(&cudaNewFrontier, sizeof(char) * vertices.size());
  cudaMemset(cudaNewFrontier, 0, sizeof(char) * vertices.size());

  int* cudaEdgesStart;
  cudaMalloc(&cudaEdgesStart, sizeof(int) * vertices.size());

  int* cudaEdgesCount;
  cudaMalloc(&cudaEdgesCount, sizeof(int) * vertices.size());

  int* cudaEdges;
  cudaMalloc(&cudaEdges, sizeof(int) * edgeCount);

  int* cudaLayeredEdgesCount;
  cudaMalloc(&cudaLayeredEdgesCount, sizeof(int) * vertices.size());

  int* cudaLayeredEdges;
  cudaMalloc(&cudaLayeredEdges, sizeof(int) * edgeCount);

  int* cudaEdgeCapacities;
  cudaMalloc(&cudaEdgeCapacities, sizeof(int) * edgeCount);

  unsigned int* cudaVertDists;
  cudaMalloc(&cudaVertDists, sizeof(int) * vertices.size());
  cudaMemset(cudaVertDists, 255, sizeof(int) * vertices.size());
  cudaMemset(cudaVertDists, 0, sizeof(int));

  char* cudaFoundSink;
  cudaMalloc(&cudaFoundSink, sizeof(char));
  cudaMemset(cudaFoundSink, 0, sizeof(char));

  char* cudaProgressed;
  cudaMalloc(&cudaProgressed, sizeof(char));

  std::vector<int> edges(edgeCount);
  std::vector<int> edgeCapacities(edgeCount);
  std::vector<int> edgesStart(vertices.size());
  std::vector<int> edgesCount(vertices.size());

  int ec = 0;
  for (int i = 0; i < vertices.size(); i++) {
    edgesStart[i] = ec;
    int numNeighbors = neighbors[i].size();
    edgesCount[i] = numNeighbors;

    int e = 0;
    for (auto &[dst, edge] : this->neighbors[i]) {
      edgeCapacities[ec + e] = edge.cap;
      edges[ec + e] = dst;
      e++;
    }
    ec += e;
  }

  int level = 0;
  cudaMemcpy(cudaEdges, edges.data(), edgeCount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(cudaEdgeCapacities, edgeCapacities.data(), edgeCount * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(cudaEdgesStart, edgesStart.data(), vertices.size() * sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(cudaEdgesCount, edgesCount.data(), vertices.size() * sizeof(int), cudaMemcpyHostToDevice);

  const int threadsPerBlock = 512;
  const int blocks = (vertices.size() + threadsPerBlock - 1) / threadsPerBlock;

  char foundSink = false;
  char progressed = false;

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

  std::vector<int> layeredEdgesCount(vertices.size());
  std::vector<int> layeredEdges(edgeCount);
  cudaMemcpy(layeredEdges.data(), cudaLayeredEdges, edgeCount * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(layeredEdgesCount.data(), cudaLayeredEdgesCount, vertices.size() * sizeof(int), cudaMemcpyDeviceToHost);
  for (int i = 0; i < vertices.size(); i++) {
    vertices[i].layered_dst_array[0].resize(layeredEdgesCount[i]);
    for (int j = 0; j < layeredEdgesCount[i]; j++) {
      vertices[i].layered_dst_array[0][j] = layeredEdges[edgesStart[i] + j];
    }
  }

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

  return foundSink;
}

