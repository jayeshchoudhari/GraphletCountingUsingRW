/***
* SRW2CSS for 4 and 5 vertex graphlets
* Algorithm:
* Random Walk on G2, which is defined as follows.
* Each node in G2 is an edge in G and two nodes
* e1 and e2 are connected in G2 if they share an
* endpoint in G.
* Simulating a random walk in G2 is quite straight forward.

* We give an algorithm for counting 5 cliques using SRW2.
* A major requirement for this algorithm is to count the
* number of edges in G2. This is quite non-trivial in the
* random walk model. For now, we assume this quantity is 
* given to us.
* 
* 
* 1. Start a random walk from a vertex u and take a step to vertex v.
* 2. Now the current edge (which is a vertex in G2) is e = {u,v}
* 2. To take another step in the random walk:
*    (i) select u w.p. d_u/(d_u+d_v), and (ii) select a neighbor u.a.r
*
* 3. Now take 4 steps of the walk. Let X^4 = (e1,e2,e3,e4)
*  Note that each two consecutive edges in X^4 share an edge.
*  We assume these four edges consists of exactly five distinct
*  vertices. If not, continue till we find four distinct vertices.
*
* 4. Now perform a random walk of given length.
* 4.a. Assume e_i is the last edge in the collection.
*      Then, to take a step, choose a neighbor of e_i
*      uniformly at random. 
*
* 5a. If the last four edges has exactly five distinct vertices,
*     and they form a five clique, then add m * d(e2) * d(e3) * d(e4) / 240
*     where m is the number of edges in G2. 
* 5b. There is a CSS option in the algorithm. If this option
*     is set, the the algorithm remains the same but the normalization
*     factor changes. 
*      If the last four edges has exactly five distinct vertices,
*     and they form a five clique, then add m /(1/d(e1)+1/d(e_2)+...+1/d(e_4))
*     where m is the number of edges in G2. 
* 
* 6. Finally, take the average over the length of the random walk.
* Note that the factor m/24 can be multiplied at the end of the algorithm 

* 5b. If it forms a triangle, and CSS flag is set, then add m/(1/d(v1)+1/(dv2)+1/d(v3)).
*/


#include "namespace.h"
#include "graphIO.h"
#include "count5Graphlets.h"
#include "utilities.h"

using namespace std;


double rwCount5Graphlets :: SRWCount5CliqueGraphlet(Graph &G, int numSteps, Count numEdgesInG2) 
{
	Count numVertices = G.getNumVertices();
	double final5CliqueCounter = 0.0;

	// Random starting point...
	// 1. Start a random walk from a vertex u and take a step to vertex v.
	VertexIdx startNode = getRandomStartPoint(numVertices);

	// create an initial 4 graphlet first...
	vector<vector<VertexIdx>> listOfRWEdges = getX3Edges(G, startNode);

	// get the degree of each edge... (this can be optimized later...)
	vector<Count> edgeDegreeList = getDegreeOfEdges(G, listOfRWEdges);

	// get set of nodes in the edges accumulated till now... specifically the last 3 edges...
	set<VertexIdx> currentSetOfNodes = getSetOfNodesFromEdges(listOfRWEdges);
	set<VertexIdx>::iterator setIt;

	for(int i = 0; i < numSteps + 1; i++)
	{
		// get the last edge to sample a new one...
		vector<VertexIdx> lastEdge = listOfRWEdges.back();
		
		// Till a new node is found...
		while(true)
		{
			// get a new uar edge... using the last edge to sample from...
			vector<VertexIdx> newEdge = sampleUARNeighboringEdge(G, lastEdge[0], lastEdge[1]);
			VertexIdx newNode = newEdge[1];
			
			// check if the new node exists in the already existing nodes...
			// we need 5 unique nodes...
			setIt = currentSetOfNodes.find(newNode);

			if(setIt == currentSetOfNodes.end())
			{
				currentSetOfNodes.insert(newNode);
				listOfRWEdges.push_back(newEdge);
				edgeDegreeList.push_back(G.getDegree(newEdge[0]) + G.getDegree(newEdge[1]) - 2);
				break;
			}
		}

		// degree of the last 3 edges...
		vector<Count> last3DegreeValues = vector<Count>(edgeDegreeList.end()-3, edgeDegreeList.end());

		// check if the current set of nodes form a clique...
		bool cliqueIdentifier = G.checkClique(currentSetOfNodes);
		if (cliqueIdentifier)
		{
			// add m * d(e2) * d(e3) * d(e4) / 24
			double qtyToAdd = numEdgesInG2;
			for(int i = 0; i < last3DegreeValues.size(); i++)
			{
				qtyToAdd = qtyToAdd * last3DegreeValues[i];
			}
			qtyToAdd = qtyToAdd/24;
			final5CliqueCounter += qtyToAdd;
		}

		// update new current set of nodes to the nodes from the last 3 edges...
		vector<vector<VertexIdx>> last3Edges = vector<vector<VertexIdx>>(listOfRWEdges.end() - 3, listOfRWEdges.end());
		currentSetOfNodes = getSetOfNodesFromEdges(last3Edges);
	}
	
	return final5CliqueCounter;
}


vector<VertexIdx> rwCount5Graphlets :: sampleUARNeighboringEdge(Graph &G, VertexIdx uNode, VertexIdx vNode)
{
	VertexIdx samplerNode, newNode;

	Count deg_of_u = G.getDegree(uNode);
	Count deg_of_v = G.getDegree(vNode);
	Count deg_of_uv = deg_of_u + deg_of_v;

	int nextSamplerNodeId = rand() % deg_of_uv;

	if(nextSamplerNodeId < deg_of_u)
		samplerNode = uNode;
	else
		samplerNode = vNode;

	Count numNbs = G.getDegree(samplerNode);
	int randNbr = rand() % numNbs;
	newNode = G.getKthNeighbor(samplerNode, randNbr);

	vector<VertexIdx> nextEdge = {samplerNode, newNode};

	return nextEdge;
}


vector<vector<VertexIdx>> rwCount5Graphlets :: getX3Edges(Graph &G, VertexIdx startNode)
{
	VertexIdx uNode = startNode;

	vector<VertexIdx> currentWalkNodes;
	currentWalkNodes.push_back(uNode);

	Count numNbs = G.getDegree(uNode);
	int randNbr = rand() % numNbs;
	VertexIdx vNode = G.getKthNeighbor(uNode, randNbr);

	currentWalkNodes.push_back(vNode);
	// 2. Now the current edge (which is a vertex in G2) is e = {u,v}
	// ePair lastPair = make_pair(uNode, vNode);
	vector<vector<VertexIdx>> XkEdges;
	vector<VertexIdx> edgeVec;
	edgeVec.push_back(uNode);
	edgeVec.push_back(vNode);
	XkEdges.push_back(edgeVec);
	VertexIdx nextNode;

	while(true)
	{
		vector<VertexIdx> uarEdge = sampleUARNeighboringEdge(G, uNode, vNode);
		nextNode = uarEdge[1];

		if(notInList(currentWalkNodes, nextNode))
		{
			uNode = uarEdge[0];
			vNode = nextNode;
			currentWalkNodes.push_back(vNode);

			// add the last edge to the list of edges...
			XkEdges.push_back(uarEdge);

			if(XkEdges.size() == 3)
				break;
		}
		else
		{
			continue;	
		}
	}

	return XkEdges;
} 

vector<Count> rwCount5Graphlets :: getDegreeOfEdges(Graph &G, vector<vector<VertexIdx>> edgeList)
{
	vector<Count> degreeList;
	for(int i = 0; i < edgeList.size(); i++)
	{
		degreeList.push_back(G.getDegree(edgeList[i][0]) + G.getDegree(edgeList[i][1]) - 2);
	}

	return degreeList;
}
