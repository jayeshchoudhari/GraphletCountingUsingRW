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
#include "include/utilities.h"

using namespace std;

double rwCount5Graphlets :: chooseNeighborOfEdge(Graph &G, VertexIdx u, VertexIdx v)
{
	vector<VertexIdx> nextEdge;

	Count deg_of_u = G.getDegree(uNode);
	Count deg_of_v = G.getDegree(vNode);
	Count deg_of_uv = deg_of_u + deg_of_v;

	int nextSamplerNodeId = rand() % deg_of_uv;

	if(nextSamplerNodeId < deg_of_u)
		nextSamplerNode = uNode;
	else
		nextSamplerNode = vNode;

	Count numNbs = adjList[nextSamplerNode].size();
	randNbr = rand() % numNbs;
	nextNode = adjList[nextSamplerNode][randNbr];



	return 
}


double rwCount5Graphlets :: getX4Edges(Graph &G, VertexIdx startNode)
{
	VertexIdx uNode = startNode;

	vector<VertexIdx> initialRandomWalkNodes;
	initialRandomWalkNodes.push_back(uNode);

	Count numNbs = adjList[uNode].size();
	int randNbr = rand() % numNbs;
	VertexIdx vNode = adjList[uNode][randNbr];

	initialRandomWalkNodes.push_back(vNode);
	// 2. Now the current edge (which is a vertex in G2) is e = {u,v}
	// ePair lastPair = make_pair(uNode, vNode);
	vector<vector<VertexIdx>> list_of_edges;
	vector<VertexIdx> edgeVec;
	edgeVec.push_back(uNode);
	edgeVec.push_back(vNode);
	listOfEdgesX4.push_back(edgeVec);
	VertexIdx nextNode;

	while(True)
	{
		Count deg_of_u = G.getDegree(uNode);
		Count deg_of_v = G.getDegree(vNode);
		Count deg_of_uv = deg_of_u + deg_of_v;

		int nextSamplerNodeId = rand() % deg_of_uv;

		if(nextSamplerNodeId < deg_of_u)
			nextSamplerNode = uNode;
		else
			nextSamplerNode = vNode;
		
		Count numNbs = adjList[nextSamplerNode].size();
		randNbr = rand() % numNbs;
		nextNode = adjList[nextSamplerNode][randNbr];
		if(nextNode notInList(initialRandomWalkNodes, nextNode))
		{
			uNode = nextSamplerNode;
			vNode = nexNode;
			initialRandomWalkNodes.push_back(vNode);

			// add the last edge to the list of edges...
			vector<VertexIdx> edgeVec;
			edgeVec.push_back(uNode);
			edgeVec.push_back(vNode);
			listOfEdgesX4.push_back(edgeVec);

			if(listOfEdgesX4.size() == 4)
				break
		}
		else
		{
			continue;	
		}

	}

	return listOfEdgesX4;
} 


double rwCount5Graphlets :: SRWCount5CliqueGraphlet(Graph &G, int numSteps) 
{
	Count numVertices = G.getNumVertices();

	// Random starting point...
	// 1. Start a random walk from a vertex u and take a step to vertex v.
	VertexIdx startNode = getRandomStartPoint(numVertices);

	// create an initial 4 graphlet first...
	vector<vector<VertexIdx>> listOfRWEdges = getX4Edges(G, startNode);
	vector<VertexIdx> lastSetOfNodes;
	set<VertexIdx> last5Nodes;
	for(int i = 0; i < listOfRWEdges; i++)
	{
		last5Nodes.insert(listOfRWEdges[0][0]);
		last5Nodes.insert(listOfRWEdges[0][1]);
	}

	for(int i = 0; i < numSteps; i++)
	{
		vector<VertexIdx> lastEdge = listOfRWEdges.back();

	}
	

}
