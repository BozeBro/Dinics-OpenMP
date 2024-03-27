TITLE: Parallelizing and Evaluating Dinic’s Max Flow Algorithm using OpenMP
URL: https://github.com/BozeBro/Dinics-OpenMP

SUMMARY:
We are going to parallelize the Dinic’s Max Flow Algorithm on the CPU using the OpenMP library and measure the scalability of our algorithm on large graph sizes.

BACKGROUND:
Flow networks are directed graphs with a source and sink node. Every edge has a nonnegative capacity assigned to it. A flow in a flow network is an assignment of values to each edge such that all values are between 0 and the edge’s capacity, and such that for all nodes (except the source and sink) the sum of values for in-edges equals that for out edges. The total capacity of the flow is defined as the sum of the out edge values from the source minus the sum of the in edge values of the sink.
Dinic’s algorithm is an algorithm to find the maximum flow of a flow network. It involves repeating two steps: computing a layered graph inside of a flow network (using BFS to get the distances of each node from the source), and finding a blocking flow in the layered graph (repeatedly invoking DFS to find paths in the graph to add to the flow value). 

THE CHALLENGE:
In general, graph algorithms can be tricky to parallelize because of shared edges/vertices, poor cache locality, and inconsistent access patterns.
Two potential areas for parallelism in Dinic’s algorithm are with the BFS and with the many DFS calls. Breadth-first search is in theory parallelizable, with each vertex in the frontier finding its neighbors in parallel. In practice, however, there are several challenges to a practical implementation. There is a significant component of synchronization and reduction - as each iteration of advancing the frontier is dependent on the previous iteration’s frontier, and as reduction is needed to compute a frontier if the computation is done in parallel. There could be other ways to compute this level graph that don’t involve a fully synchronized BFS too that could be worth testing.
Dinic’s classically relies on running many DFS sequentially to find blocking flows, but there are probably ways to parallelize these searches, likely involving making them operate on disjoint edges if possible.

RESOURCES:
The plan is to use the GHC machines for most of our development and preliminary performance measurements. There will be usage of the PSC machines to measure parallelizing on larger thread counts. 

GOALS AND DELIVERABLES
PLAN TO ACHIEVE
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
HOPE TO ACHIEVE
CUDA implementation of Parallel Dinic’s algorithm
Push-Relabel parallelization & comparison

PLATFORM CHOICE
The machines we are using have 8 processor CPU cores for GHC and at least 128 processor cores on the PSC machines. We plan to use C++ for parallelization as we will have access to OpenMP. Also, for a large computer like PSC, we will get to see how this algorithm handles graphs that have high memory requirements.

SCHEDULE
By April 1 - Finalize design of code and get initial layer of code for sequential Dinic’s
By April 8 - Create sequential dinic’s algorithm and create some inputs with varied graph characteristics
By April 15 - Finish sequential algorithm, input generation, and generate results for sequential algorithm.
By April 22 - Explore parallelizing BFS (layered graph)
By April 29 - Explore parallelizing blocking flow generation (DFS calls)
Final - Measure different versions of the project compared to the sequential version.
