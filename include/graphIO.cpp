#include "namespace.h"
#include "graphIO.h"
#include "utilities.h"
// #include <bits/stdc++.h>

using namespace std;

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

		graphFile.close();
		
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


vector<OrderedEdge> Graph :: getUniformRandomEdges(Count numEdgesToSample)
{
	std::default_random_engine generator;
  	std::uniform_int_distribution<int> distribution(0, nEdges-1);
    vector<OrderedEdge> localOrdRwEdges;

	for(unsigned int i = 0; i < numEdgesToSample; i++)
	{
		int edgeIdToSample = distribution(generator);
		vector<VertexIdx> selectedEdge = edgeList[edgeIdToSample];

        VertexIdx u, v, low_deg;
        if((nodeDeg[selectedEdge[0]] <= nodeDeg[selectedEdge[1]]) || (nodeDeg[selectedEdge[0]] <= nodeDeg[selectedEdge[1]] && selectedEdge[0] < selectedEdge[1]))
        {
            u = selectedEdge[0];
            v = selectedEdge[1];
            low_deg = nodeDeg[selectedEdge[0]];
        }
        else
        {
            u = selectedEdge[1];
            v = selectedEdge[0];
            low_deg = nodeDeg[selectedEdge[1]];
        } 
        // cout << u << " " << v << " " << nodeDeg[u] << " " << nodeDeg[v] << " " << low_deg << "\t";
        localOrdRwEdges.push_back(OrderedEdge{u, v, low_deg});

	}
	// cout << "\n";
	// return localRwEdges;
	return localOrdRwEdges;
}


bool Graph :: checkEdgeInAdjList(VertexIdx v1, VertexIdx v2)
{
    bool boolVar; 
    // vector<VertexIdx> nbrs = adjList[v1];
    boolVar = binary_search (adjList[v1].begin(), adjList[v1].end(), v2);
    return boolVar;
}

int Graph :: checkEdgeInAdjListInt(VertexIdx v1, VertexIdx v2)
{
    bool boolVar; 
    // vector<VertexIdx> nbrs = adjList[v1];
    boolVar = binary_search(adjList[v1].begin(), adjList[v1].end(), v2);
    if(boolVar)
    {
    	return 1;
    }
    else
    {
    	return 0;
    }
}


// check if 5 clique
bool Graph :: checkClique(set<VertexIdx> setOf5Nodes)
{
	bool boolVar;
	int innerFlag = 0;

	for(auto it1 = setOf5Nodes.begin(); it1 != setOf5Nodes.end(); it1++)
    {
    	innerFlag = 0;
        for(auto it2 = std::next(it1); it2 != setOf5Nodes.end(); it2++)
        {
        	bool checkEdgeBool = checkEdgeInAdjList(*it1, *it2);
			if(!checkEdgeBool)
			{
				innerFlag = 1;
				break;
			}
            // std::cout << *it1 << " " << *it2 << " ";
        }
        // std::cout << "\n";
        if (innerFlag == 1)
			break;
    }

	if(innerFlag)
		return 0;
	else
		return 1;
}


// check if the last vertex in the list of arguments is connected to at least 2 of other 3 vertices...
bool Graph :: checkConnectionOfXToAny2OfUVW(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode)
{
	bool isConnected = 0;
	bool uxEdge = checkEdgeInAdjList(uNode, xNode);

	if(uxEdge)
	{
		bool vxEdge = checkEdgeInAdjList(vNode, xNode);

		if(vxEdge)
		{
			isConnected = 1;				// uxEdge and vxEdge
		}
		else
		{
			bool wxEdge = checkEdgeInAdjList(wNode, xNode);
			if(wxEdge)
			{
				isConnected = 1;			// uxEdge and wxEdge
			}
		}
	}
	else
	{
		bool vxEdge = checkEdgeInAdjList(vNode, xNode);

		if(vxEdge)
		{
			bool wxEdge = checkEdgeInAdjList(wNode, xNode);
			if(wxEdge)
			{
				isConnected = 1;			// vxEdge and wxEdge
			}
		}
	}

	return isConnected;
}

int Graph :: checkConnectionOfXToOtherThree(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode)
{
	int conns = 0;

	int uxconn = checkEdgeInAdjListInt(uNode, xNode);
	int vxconn = checkEdgeInAdjListInt(vNode, xNode);
	int wxconn = checkEdgeInAdjListInt(wNode, xNode);

	conns = uxconn + vxconn + wxconn;

	return conns;
}


vector<VertexIdx> Graph :: checkConnectionXNotConnectedTo(VertexIdx uNode, VertexIdx vNode, VertexIdx wNode, VertexIdx xNode)
{
	// int conns = 0;
	vector<VertexIdx> notConnectedTo;

	int uxconn = checkEdgeInAdjListInt(uNode, xNode);
	int vxconn = checkEdgeInAdjListInt(vNode, xNode);
	int wxconn = checkEdgeInAdjListInt(wNode, xNode);

	if(uxconn == 0)
	{
		notConnectedTo.push_back(uNode);
	}

	if(vxconn == 0)
	{
		notConnectedTo.push_back(vNode);
	}

	if(wxconn == 0)
	{
		notConnectedTo.push_back(wNode);
	}

	// conns = uxconn + vxconn + wxconn;

	return notConnectedTo;
}


vector<VertexIdx> Graph :: checkLastNodeNotConnectedTo(vector<VertexIdx> setOfNodes)
{
	// int conns = 0;
	vector<VertexIdx> notConnectedTo;

	VertexIdx lastNode = setOfNodes.back();


	for(int i = 0; i < setOfNodes.size() - 1; i++)
	{
		int uxconn = checkEdgeInAdjListInt(setOfNodes[i], lastNode);
		if(uxconn == 0)
		{
			notConnectedTo.push_back(setOfNodes[i]);
		}

		if(notConnectedTo.size() > 1)
		{
			break;
		}
	}

	return notConnectedTo;
}

Count Graph :: getCombinedNeighborSize(VertexIdx uNode, VertexIdx vNode)
{
	vector<VertexIdx> uNeighbors = getNeighbors(uNode);
	vector<VertexIdx> vNeighbors = getNeighbors(vNode);

	vector<VertexIdx> uvNeighborsVector = getUnionOfNeighbors(uNeighbors, vNeighbors);
	Count numNeighbors = (Count)uvNeighborsVector.size();
	return numNeighbors;
}
/*
struct pivotNeighborsAndSizes_3_2 Graph :: getDegreeAndNeighborsOf2Qset(vector<VertexIdx>tempComponent, int q)
{
	
	VertexIdx uNode, vNode, wNode;
	uNode = tempComponent[0];
	vNode = tempComponent[1];
	wNode = tempComponent[2];

	vector<VertexIdx> uNeighbors = getNeighbors(uNode);
	vector<VertexIdx> vNeighbors = getNeighbors(vNode);
	vector<VertexIdx> wNeighbors = getNeighbors(wNode);

	vector<VertexIdx> uvNeighbors, vwNeighbors, uwNeighbors;
	uvNeighbors = getUnionOfNeighbors(uNeighbors, vNeighbors);
	vwNeighbors = getUnionOfNeighbors(vNeighbors, wNeighbors);
	uwNeighbors = getUnionOfNeighbors(uNeighbors, wNeighbors);

	unsigned int uvSize, vwSize, uwSize;
	uvSize = uvNeighbors.size();
	uwSize = uwNeighbors.size();
	vwSize = vwNeighbors.size();

	vector<unsigned int> sizeVector;
	sizeVector.push_back(uvSize); 
	sizeVector.push_back(uwSize); 
	sizeVector.push_back(vwSize);

	// cout << "returning -- " << sizeVector[0] << " " << sizeVector[1] << " " << sizeVector[2] << "\n";
	// cout << "created vector of sizes -- " << uvSize << " " << uwSize << " " << vwSize << "\n";

	struct pivotNeighborsAndSizes_3_2 neigborsAndSizes;

	if((uvSize <= vwSize) && (uvSize <= uwSize))
	{
		neigborsAndSizes.pivotNeighbors = uvNeighbors;
		neigborsAndSizes.neighborSizes = sizeVector;
		// struct pivotNeighborsAndSizes_3_2 neigborsAndSizes = {uvNeighbors, sizeVector};
		// return neigborsAndSizes;
	}
	else if((vwSize <= uwSize) && (vwSize <= uvSize))
	{
		neigborsAndSizes.pivotNeighbors = vwNeighbors;
		neigborsAndSizes.neighborSizes = sizeVector;
		// struct pivotNeighborsAndSizes_3_2 neigborsAndSizes = {vwNeighbors, sizeVector};
		// return neigborsAndSizes;
	}
	else if((uwSize <= vwSize) && (uwSize <= uvSize))
	{
		neigborsAndSizes.pivotNeighbors = uwNeighbors;
		neigborsAndSizes.neighborSizes = sizeVector;
		// struct pivotNeighborsAndSizes_3_2 neigborsAndSizes = {uwNeighbors, sizeVector};
		// return neigborsAndSizes;
	}

	return neigborsAndSizes;
}
*/


struct pivotNeighborsAndSizes_X_2 Graph :: getDegreeAndNeighborsOf2Qset(vector<VertexIdx> tempComponent, int q)
{
	struct pivotNeighborsAndSizes_X_2 neigborsAndSizes;

	vector<vector<VertexIdx>> neighborsOfPairOfVertices;

	int minNeighborSize = getNumVertices();

	vector<unsigned int> sizeVector;

	for(int i = 0; i < tempComponent.size(); i++)
	{
		vector<VertexIdx> uNeighbors = getNeighbors(tempComponent[i]);
		for(int j = i+1; j < tempComponent.size(); j++)
		{
			// cout << i << j << endl;
			vector<VertexIdx> vNeighbors = getNeighbors(tempComponent[j]);
			vector<VertexIdx> neighborsSet = getUnionOfNeighbors(uNeighbors, vNeighbors);
			sizeVector.push_back(neighborsSet.size());

			if(minNeighborSize > neighborsSet.size())
			{
				minNeighborSize = neighborsSet.size();
				neigborsAndSizes.pivotNeighbors = neighborsSet;
			}
		}
	}

	neigborsAndSizes.neighborSizes = sizeVector;

	return neigborsAndSizes;
}