#ifndef GRAPH_IO_H
#define GRAPH_IO_H


class Graph 
{
    public:

    	Count nVertices; 						//number of vertices in the graph
        Count nEdges;     						//number of edges in this list
        int smoothNessParam; 					// value of q = i + 1 - |cut(seg_i)|
        vector <VertexIdx> nodeList;
        vector <vector<VertexIdx>> adjList;		// adj List
        vector <vector<VertexIdx>> edgeList;
        // unordered_map <VertexIdx, Count> nodeDeg;
        vector <Count> nodeDeg;

        vector<VecMat> seg3IsoGraphs, seg4IsoGraphs, seg5IsoGraphs; 
		unsigned int seg3NumEdges, seg4NumEdges, seg5NumEdges;

        /* Count Lists... */
        // vector<int> g3s(3);
        vector <Count> g3s = vector<Count>(3);

		Graph(string fileName);
		int initializeAdjList(Count numV, Count numE);

		Count getNumVertices();
		Count getNumEdges();
        Count getDegree(VertexIdx u);
        vector<VertexIdx> getNeighbors(VertexIdx u);
        VertexIdx getKthNeighbor(VertexIdx u, int k);
    	bool checkEdgeInAdjList(VertexIdx v1, VertexIdx v2);
		bool checkConnectionOfXToAny2OfUVW(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode);

		int checkEdgeInAdjListInt(VertexIdx v1, VertexIdx v2);
		int checkConnectionOfXToOtherThree(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode);
        vector<VertexIdx> checkConnectionXNotConnectedTo(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode);
        vector<VertexIdx> checkLastNodeNotConnectedTo(vector<VertexIdx> setOfNodes);


		// struct pivotNeighborsAndSizes_3_2 getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q);
		struct pivotNeighborsAndSizes_X_2 getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q);
		Count getCombinedNeighborSize(VertexIdx uNode, VertexIdx vNode);
        
        int printGraphDetails();
        
		// RW Functions...
		vector<OrderedEdge> getAllEdgesFromRStepRandomWalk(Count numSteps, VertexIdx startNode);
		// ePair lStepRandomWalk(Count, VertexIdx);

		// UAR Edge Sample
		vector<OrderedEdge> getUniformRandomEdges(Count numEdgesToSample);



    	int setQ(int qVal);
    	int CountG3s();

		// vector<VertexIdx> getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q);
		

		vector<VertexIdx> createSetOfNodesFromEdges(vector<OrderedEdge> subGraphlet);
		vector<VertexIdx> getMinNeighborhoodSet(vector<vector<int>> allCombsOfQNodes, vector<VertexIdx> subGraphletNodes);
		// vector<vector<VertexIdx>> getPivotNeighborSetsForAllSubGraphlets(vector<vector<ePair>> R_j);
		vector<vector<VertexIdx>> getPivotNeighborSetsForAllSubGraphlets(vector<vector<OrderedEdge>> R_j);
		// vector<ePair> sampleGraphletAndNeighbor(vector<i64> distValsPrev, vector<vector<VertexIdx>> prevPivotNeighborSets, vector<vector<ePair>> prevRj);
		pair<VertexIdx, vector<OrderedEdge>> sampleGraphletAndNeighbor(vector<i64> distValsPrev, vector<vector<VertexIdx>> prevPivotNeighborSets, vector<vector<OrderedEdge>> prevRj);
		sampledGraphletInfo sampleGraphletInfoAndNeighbor(vector<i64> distValsPrev, vector<vector<VertexIdx>> prevPivotNeighborSets, vector<vector<OrderedEdge>> prevRj);
		// int checkPivotSetContainsSampledNode(vector<ePair> subGraphlet, VertexIdx sampledNeighborNode);
		int checkPivotSetContainsSampledNode(vector<OrderedEdge> subGraphlet, VertexIdx sampledNeighborNode);
		int checkTrianglePivotSetContainsSampledNode(vector<OrderedEdge> subGraphlet, VertexIdx sampledNeighborNode);
		int checkSampledNodeDegOrdering(vector<OrderedEdge> subGraphlet, VertexIdx sampledNeighborNode, VertexIdx lastNode);
		pair<vector<VertexIdx>, int> getMinNeighborhoodSetForAssignment(vector<vector<int>> allCombsOfQNodes, vector<VertexIdx> subGraphletNodes, VertexIdx sampledNeighborNode);

		// utilities for building graphlet R_j...
};


#endif 		// GRAPH_IO_H
