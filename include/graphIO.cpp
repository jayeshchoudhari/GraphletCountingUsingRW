#include "namespace.h"
#include "graphIO.h"
#include <bits/stdc++.h>

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