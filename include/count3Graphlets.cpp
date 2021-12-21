#include "namespace.h"
#include "graphIO.h"
#include "count3Graphlets.h"
// #include <bits/stdc++.h>

using namespace std;

// double rwCount3Graphlets :: countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges)
prevSubGraphletSet rwCount3Graphlets :: countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc)
{
	vector<double> dRVals = {0, 0};
	vector<double> ljVals = {0, 0, rwEdges.size() * 1.0};
	Count numEdges = G.getNumEdges();

	random_device rd;
	mt19937 gen(rd());

	Count l2 = rwEdges.size();

	vector<VertexIdx> edge_degree_list(l2);
	double dR2 = 0.0;
	for (unsigned int i = 0; i < l2; i++) 
	{
		edge_degree_list[i] = rwEdges[i].degree;
		dR2 += rwEdges[i].degree;
	}

	vector<vector<VertexIdx>> nextLevelComponents;
	vector<VertexIdx> nextLevelDegrees;
	double dR3 = 0.0;

	// int subsample_size = l2/20;
	int subsample_size = l2 * l3Perc/100.0;
	ljVals.push_back(subsample_size);

	double X = 0, Y = 0, Z = 0;

	discrete_distribution<int> distribution (edge_degree_list.begin(), edge_degree_list.end());
	
	for(int s = 0; s < subsample_size; s++)
	{
		int sampledId = distribution(gen);
		VertexIdx uNode = rwEdges[sampledId].u;
		Count deg_of_u = rwEdges[sampledId].degree;
		uniform_int_distribution<int> distNbor(0, deg_of_u - 1);
		int rdNbrId = distNbor(gen);
		VertexIdx wNode = G.getKthNeighbor(uNode, rdNbrId);

		VertexIdx vNode = rwEdges[sampledId].v;
		Count deg_of_v = G.getDegree(vNode);
		Count deg_of_w = G.getDegree(wNode);

		if((deg_of_w > deg_of_v || (deg_of_w == deg_of_v && wNode > vNode)))
		{
			bool vwEdge = G.checkEdgeInAdjList(vNode, wNode);
			if(vwEdge)
			{
				Z = 1;
				vector<VertexIdx> tempComponent {uNode, vNode, wNode};
				nextLevelComponents.push_back(tempComponent);
				nextLevelDegrees.push_back(deg_of_u);
				dR3 += deg_of_u;
			}
			else
			{
				Z = 0;
			}
		}
		else
		{
			Z = 0;
		}
		Y += Z;
	}

	X =  (numEdges / ljVals[2]) * (dR2/ljVals[3]) * Y;

	struct prevSubGraphletSet threeSubGraphlets;
	threeSubGraphlets.graphletEstimate = X;
	threeSubGraphlets.graphletsForNextLevel = nextLevelComponents;
	threeSubGraphlets.degreesForNextLevel = nextLevelDegrees;
	threeSubGraphlets.totalDegree = dR3;

	// return X;
	return threeSubGraphlets;
}
