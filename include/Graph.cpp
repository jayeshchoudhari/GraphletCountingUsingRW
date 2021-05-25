#include "include/Graph.h"
// #include "include/JointSort.h"
#include <algorithm>

using namespace rwNameSpace;


//Basic binary search procedure
// Input: pointer array, index of last entry end, and val to search for
// Output: index if val is found, and -1 otherwise

// VertexIdx rwNameSpace::binarySearch(vector<EdgeIdx> array, VertexIdx end, EdgeIdx val) 
// {
//     VertexIdx low = 0;
//     VertexIdx high = end - 1;
//     VertexIdx mid;

//     while (low <= high) {
//         mid = (low + high) / 2;
//         if (array[mid] == val)
//             return mid;
//         if (array[mid] > val)
//             high = mid - 1;
//         if (array[mid] < val)
//             low = mid + 1;
//     }
//     return -1;
// }


// comparator that only compares the first in pair

// bool rwNameSpace::pairCompareFirst(Pair firstPair, Pair nextPair) {
//     return firstPair.first < nextPair.first;
// }

// int rwNameSpace :: Graph :: printGraphDetails()
// {
// 	std::cout << "num of nodes = " << nVertices << " num of edges = " << nEdges << "\n";
// 	return 0;
// }


// rwNameSpace :: Graph :: Graph(vector <VertexIdx> &v, vector<ePair> &e)
// // Graph::Graph(vector <VertexIdx> &v, vector<ePair> &edList)
// {
// 	nVertices = v.size();
// 	nEdges = e.size();

// 	// std::cout << "num of nodes = " << nVertices << " num of edges = " << nEdges << "\n";

// 	vector<VertexIdx>::iterator vIt; 

// 	for(vIt = v.begin(); vIt != v.end(); vIt++)
// 	{
// 		nodeList.push_back(*vIt);
// 	}
	
// 	// std::cout << "got vertices \n"; 

// 	adjList.resize(nVertices);

// 	// vector<ePair>::iterator edgeIt;

// 	for(const ePair &edgeIt : e)
// 	{
// 		VertexIdx src = edgeIt.first;
// 		VertexIdx dest = edgeIt.second;

// 		adjList[src].push_back(dest);
// 		adjList[dest].push_back(src);
// 	}
// }

// int rwNameSpace::printGraphDetails()