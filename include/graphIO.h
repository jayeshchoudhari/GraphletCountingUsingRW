class Graph 
{
    private:
        Count nVertices; 				//number of vertices in the graph
        Count nEdges;     				//number of edges in this list
        int smoothNessParam; 			// value of q = i + 1 - |cut(seg_i)|
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


    public:

		Graph(string fileName);
		int initializeAdjList(Count numV, Count numE);

		Count getNumVertices();
		Count getNumEdges();
        Count getDegree(VertexIdx u);
        vector<VertexIdx> getNeighbors(VertexIdx u);
        VertexIdx getKthNeighbor(VertexIdx u, int k);
    	bool checkEdgeInAdjList(VertexIdx v1, VertexIdx v2);
        int printGraphDetails();
        
        // RW Functions...
		vector<OrderedEdge> getAllEdgesFromRStepRandomWalk(Count numSteps, VertexIdx startNode);
		// ePair lStepRandomWalk(Count, VertexIdx);

    	int setQ(int qVal);
    	int CountG3s();

        
		// vector<VertexIdx> getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q);
        struct pivotNeighborsAndSizes_3_2 getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q);
		Count getCombinedNeighborSize(VertexIdx uNode, VertexIdx vNode);
		
		bool checkConnectionOfXToAny2OfUVW(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode);

		double buildCliqueGraphlet(vector<OrderedEdge> rwEdges, int maxGraphLetSize);
		// vector<VertexIdx> createSetOfNodesFromEdges(vector<ePair> subGraphlet);
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


Graph :: Graph(string graphFileName)		// constructor
{
	// nVertices = v.size();
	vector<ePair> allEdges;

	ifstream graphFile;

	int lineNum = 0;
	Count numVertices, numEdges;
	VertexIdx u, v;

	stringstream ss;
	string line;

	graphFile.open(graphFileName, ifstream::in);

	if(graphFile.is_open())
	{
		if (getline(graphFile, line))
		{
			ss.clear();
			ss.str("");
			ss << line;
			ss >> numVertices  >> numEdges;

			initializeAdjList(numVertices, numEdges);
			
			while(getline(graphFile, line))
			{
				ss.clear();
				ss.str("");
				ss << line;
				ss >> u >> v;
				adjList[u].push_back(v);
				adjList[v].push_back(u);

				// populate node degress
				nodeDeg[u] += 1;
				nodeDeg[v] += 1;

				vector <VertexIdx> tempEdge;
				tempEdge.push_back(u);
				tempEdge.push_back(v);
				edgeList.push_back(tempEdge);
				// allEdges.push_back(make_pair(u, v));
				lineNum += 1;
			}
		}

		// store sorted vectors....
		for(unsigned int i = 0; i < adjList.size(); i++)
		{
			sort(adjList[i].begin(), adjList[i].end());
		}
		cout << "Done with Initializing and graph in memory now...\n";
	}
	else
	{
		cout << "Cannot open the said file...\n";
		exit(0);
	}
	
}

int Graph :: initializeAdjList(Count numVertices, Count numEdges)
{
	nVertices = numVertices + 1;
	nEdges = numEdges;

	// for(vIt = v.begin(); vIt != v.end(); vIt++)
	for(Count i = 0; i < nVertices; i++)
	{
		// nodeList.push_back(*vIt);
		nodeList.push_back((VertexIdx)i);
	}
	
	// std::cout << "got vertices \n"; 
	adjList.resize(nVertices);
    nodeDeg.resize(nVertices);

	return 0;
}

Count Graph :: getNumVertices()
{
	return nVertices;
}

Count Graph :: getNumEdges()
{
	return nEdges;
}

Count  Graph :: getDegree(VertexIdx u)
{
    return nodeDeg[u];
}

vector<VertexIdx> Graph :: getNeighbors(VertexIdx u)
{
    return adjList[u];
}

VertexIdx Graph :: getKthNeighbor(VertexIdx u, int k)
{
    return adjList[u][k];
}

int Graph :: printGraphDetails()
{
	std::cout << "num of nodes = " << nVertices << " num of edges = " << nEdges << "\n";
    map<int, int> degdist;
    Count maxDeg = 0;
	for(unsigned int i = 0; i < nVertices; ++i)
	{
		// cout << i << " : " << nodeDeg[i] << "\n";
        Count deg = nodeDeg[i];
        degdist[deg/10] += 1;
        if(deg>maxDeg) 
        {
            maxDeg = deg;
        }
	}

    map<int, int>::iterator mIt;
	for(mIt = degdist.begin(); mIt != degdist.end(); ++mIt)
	{
		std::cout << mIt->first << " : " << mIt->second << "\n";
	}

	cout << "Maxdeg = " << maxDeg << endl;
    /*
	// Print graph  //
	for (int i = 0; i < nVertices; i++)
	{
		// print current vertex number
		cout << i << " --> ";

		// print all neighboring vertices of vertex i
		for (int v : adjList[i])
			cout << v << " ";
		cout << endl;
	}
	//
    */ 
	return 0;
}

vector<OrderedEdge> Graph :: getAllEdgesFromRStepRandomWalk(Count numSteps, VertexIdx startNode)
{
	vector<ePair> localRwEdges;
    vector<OrderedEdge> localOrdRwEdges;

	vector<VertexIdx> rWalk;
	rWalk.push_back(startNode);
	VertexIdx currentNode = startNode;
	VertexIdx nextNode = startNode;

	cout << "Random Walk: " << startNode << "\t";
	ePair lastEdge;

	for(unsigned int i = 0; i < numSteps; i++)
	{
		Count numNbs = adjList[currentNode].size();
		int randNbr = rand() % numNbs;
		nextNode = adjList[currentNode][randNbr];
		rWalk.push_back(nextNode);
		lastEdge = make_pair(currentNode, nextNode);
		localRwEdges.push_back(lastEdge);
        VertexIdx u, v, low_deg;
        if((nodeDeg[currentNode] <= nodeDeg[nextNode]) || (nodeDeg[currentNode] <= nodeDeg[nextNode] && currentNode < nextNode))
        {
            u = currentNode;
            v = nextNode;
            low_deg = nodeDeg[currentNode];
        }
        else
        {
            u = nextNode;
            v = currentNode;
            low_deg = nodeDeg[nextNode];
        } 
        // cout << u << " " << v << " " << nodeDeg[u] << " " << nodeDeg[v] << " " << low_deg << "\t";
        localOrdRwEdges.push_back(OrderedEdge{u, v, low_deg});

		currentNode = nextNode;
		// cout << nextNode << "\t";
	}
	// cout << "\n";
	// return localRwEdges;
	return localOrdRwEdges;
}

bool Graph :: checkEdgeInAdjList(VertexIdx v1, VertexIdx v2)
{
    bool boolVar; 
    vector<VertexIdx> nbrs = adjList[v1];
    boolVar = binary_search (nbrs.begin(), nbrs.end(), v2);
    return boolVar;
}