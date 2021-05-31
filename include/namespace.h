#ifndef NAMESPACE_H
#define NAMESPACE_H

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <map>
#include <unordered_map>
#include <random>
#include <sstream>
#include <algorithm>
// #include <utility>

using namespace std;

using int32     = int;
using ui        = unsigned int;
using i64 		= int64_t;
using VertexIdx = ui;
using EdgeIdx   = ui;
using Count     = ui;
using ePair 	= pair<VertexIdx, VertexIdx>;
using VecMat 	= vector< vector<VertexIdx> >;


const EdgeIdx invalidEdge = -1;
// A structure to store an edge with all the relavant information
struct OrderedEdge {
    VertexIdx u;  // the lower degree end point
    VertexIdx v; // the higher degree end point
    // EdgeIdx index;   // nbors[index] represents (src, dest)
    VertexIdx degree; // deg(u)
};

struct sampledGraphletInfo {
    Count graphletId;
    VertexIdx sampledNode;
    vector<OrderedEdge> extendedGraphlet;
};

struct pivotNeighborsAndSizes_3_2{
    vector<VertexIdx> pivotNeighbors;
    vector<unsigned int> neighborSizes;
};

#endif      // NAMESPACE_H
