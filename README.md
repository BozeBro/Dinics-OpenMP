### TITLE
Parallelizing and Evaluating Dinic’s Max Flow Algorithm using OpenMP
### URL
https://github.com/BozeBro/Dinics-OpenMP

### SUMMARY:
We are going to parallelize the Dinic’s Max Flow Algorithm on the CPU using the OpenMP library and measure the scalability of our algorithm on large graph sizes.

### BACKGROUND:
Flow networks are directed graphs with a source and sink node. Every edge has a nonnegative capacity assigned to it. A flow in a flow network is an assignment of values to each edge such that all values are between 0 and the edge’s capacity, and such that for all nodes (except the source and sink) the sum of values for in-edges equals that for out edges. The total capacity of the flow is defined as the sum of the out edge values from the source minus the sum of the in edge values of the sink.
Dinic’s algorithm is an algorithm to find the maximum flow of a flow network. It involves repeating two steps: computing a layered graph inside of a flow network (using BFS to get the distances of each node from the source), and finding a blocking flow in the layered graph (repeatedly invoking DFS to find paths in the graph to add to the flow value). 

### THE CHALLENGE:
In general, graph algorithms can be tricky to parallelize because of shared edges/vertices, poor cache locality, and inconsistent access patterns.
Two potential areas for parallelism in Dinic’s algorithm are with the BFS and with the many DFS calls. Breadth-first search is in theory parallelizable, with each vertex in the frontier finding its neighbors in parallel. In practice, however, there are several challenges to a practical implementation. There is a significant component of synchronization and reduction - as each iteration of advancing the frontier is dependent on the previous iteration’s frontier, and as reduction is needed to compute a frontier if the computation is done in parallel. There could be other ways to compute this level graph that don’t involve a fully synchronized BFS too that could be worth testing.
Dinic’s classically relies on running many DFS sequentially to find blocking flows, but there are probably ways to parallelize these searches, likely involving making them operate on disjoint edges if possible.

### RESOURCES:
The plan is to use the GHC machines for most of our development and preliminary performance measurements. There will be usage of the PSC machines to measure parallelizing on larger thread counts. 

### GOALS AND DELIVERABLES
#### PLAN TO ACHIEVE
Gain interesting insights into how different graph geometries affect algorithm performance
Working Parallel implementation of Dinic’s algorithm
Parallelize BFS
Parallelize blocking flow (groups of DFS calls)
Test on sparse and dense graphs
Predict better performance on dense graphs compared to sparse graphs
Raw dinic’s runs in O(n^2m)
Estimate speedup is better for graphs with shorter distance from source to sink
Parallel BFS is dependent on distance from source to sink.
Vary flow numbers (unit networks, huge flows)
We expect a lot of variation in speedup based on the graph
We want to try several approaches for each part of the algorithm and analyze how they perform on various graphs to get insight about what works for what
#### HOPE TO ACHIEVE
CUDA implementation of Parallel Dinic’s algorithm
Push-Relabel parallelization & comparison

### PLATFORM CHOICE
The machines we are using have 8 processor CPU cores for GHC and at least 128 processor cores on the PSC machines. We plan to use C++ for parallelization as we will have access to OpenMP. Also, for a large computer like PSC, we will get to see how this algorithm handles graphs that have high memory requirements.

### SCHEDULE
By April 1 - Finalize design of code and get initial layer of code for sequential Dinic’s
By April 8 - Create sequential dinic’s algorithm and create some inputs with varied graph characteristics
By April 15 - Finish sequential algorithm, input generation, and generate results for sequential algorithm.
By April 22 - Explore parallelizing BFS (layered graph)
By April 29 - Explore parallelizing blocking flow generation (DFS calls)
Final - Measure different versions of the project compared to the sequential version.

# Milestone Report

## Sequential results on Large Graphs
Device Information

| System Info | Value                  |
| ----------- | ---------------------- |
| OS          | Mac OS (Sonoma 14.1.1) |
| Chip        | M1 (2020)              |
| Memory      | 16 GB                  | 


Node Count: 50k
Edge Count: 50M
Max Capacity available: 10k
`java generator/GraphGenerator.java 0 random 50000 50000000 10000`
| Data Keys           | Values     |
| ------------------- | ---------- |
| Initialization time | 131.585684 |
| Compute time        | 63.106117  |
| Total time          | 194.691801 |
| BFS time            | 60.385660  |
| DFS time            | 0.763157   |
| Flow Value          | 4715076    | 



Results when run on a very dense graph with 4000 vertices and 1000 max flow:
(M1 Max Macbook Pro, 2021)
java GraphGenerator.java 0 random 4000 16000000 1000
| Data Keys           | Values    |
| ------------------- | --------- |
| Initialization time | 10.669528 |
| Compute time        | 4.623483  |
| Total time          | 15.293011 |
| BFS time            | 3.525797  |
| DFS time            | 0.688907  |
| Flow Value          | 1253545   | 

## Work Completed so far
Currently, we have implemented the sequential algorithm for Dinic’s algorithm. We took measurements of the computation time of the BFS and DFS usage within the algorithm, as well as initialization time and reference flow value. We created a tool to export information about the graph into GraphViz format, which can be used to visualize the graph - useful for debugging. We created a graph generator that supports generating two kinds of graphs in various sizes - ones with random edges between vertices, and ones with guaranteed connectivity. We created initial tests that contained a large number of vertices, edges, and edge capacities. We started work on parallelizing the BFS portion of the Dinic’s algorithm. From previous benchmarking, we found that BFS takes most of the computation time and shows to be the bottleneck, taking between 80% of the time (on very dense graphs) to 99% (on most graphs) of the time.

## Change of Plans
Given how much of a bottleneck BFS is in our measurements, we plan to spend more time finding faster parallel BFS algorithms to speed up Dinic’s instead of equally exploring DFS and BFS approaches. Our current approach is to use atomic operations with OpenMP to create a lock-free BFS, but we are also considering exploring using CUDA and/or a per-vertex approach (as opposed to a per-frontier approach). If we are able to speed up BFS by 5x, we may still consider optimizing DFS (because in that case it would be more of a bottleneck), but otherwise, it may be more productive for us to explore other BFS approaches or other algorithms such as Push-Relabel.

## Schedule Updated
April 21 - Get Parallel (lock free queue) BFS working
April 23 - Make progress for CUDA BFS and per-vertex BFS approach
April 30 - Complete CUDA and per-vertex approach, perhaps start attempting push-relabel
May 6 (final) - Finish Creating Presentation Poster, finish results for different BFS implementations, visualize results and speedup graphs.

## Deliverables
Visualization of different implementation results and speedups of Dinic’s algorithm using different BFS implementations running with Mac OS, Gates machines, and PSC machines.
Show time taken of BFS vs DFS for each implementation.
Plots for computation time over different graph inputs and techniques

## Concerns
Not too sure if it is more interesting to examine different types of BFS implementations vs. making a parallel push-relabel algorithm. 
