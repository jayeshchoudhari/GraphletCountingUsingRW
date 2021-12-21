#include "namespace.h"
#include "graphIO.h"
#include "countCliques.h"
#include <bits/stdc++.h>

using namespace std;

prevSubGraphletSet  rwCountXCliqueGraphlets :: countXCliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges, vector<int> percSamples, int k)
{
	vector<double> dRVals = {0.0, 0.0};
	vector<double> ljVals = {0.0, 0.0, rwEdges.size() * 1.0};
	Count numEdges = G.getNumEdges();

	random_device rd;
	mt19937 gen(rd());

	Count l2 = rwEdges.size();
	dRVals.push_back(0.0);		// index 2

	double finalEstimate;

	vector<vector<VertexIdx>> nextLevelComponents, newNextLevelComponents;
	vector<VertexIdx> nextLevelDegrees, newNextLevelDegrees;

	for (unsigned int i = 0; i < l2; i++)
	{
		nextLevelDegrees.push_back(rwEdges[i].degree);
		dRVals[2] += rwEdges[i].degree;
	}

	vector<double> allEstimates;
	// First two estimates for clique size 1 and 2...
	// Ideally should be number of nodes and number of edges..
	// For now keeping 0 as placeholder... and thus shouldn't be used...
	allEstimates.push_back(0.0);
	allEstimates.push_back(0.0);


	for(int i = 2; i < k; i++)
	{
		dRVals.push_back(0.0);		// index i+1

		int subsample_size = ljVals[i] * percSamples[i-2]/100.0;
		ljVals.push_back(subsample_size);

		double X = 1.0, Y = 0.0, Z = 0.0;

		discrete_distribution<int> distribution (nextLevelDegrees.begin(), nextLevelDegrees.end());

		for(int s = 0; s < subsample_size; s++)
		{
			int sampledId = distribution(gen);
			vector<VertexIdx> seqNodesFromComponent;
			vector<Count> seqNodeDegFromComponent;

			if(i == 2)
			{
				VertexIdx uNode = rwEdges[sampledId].u;
				Count deg_of_u = rwEdges[sampledId].degree;
				VertexIdx vNode = rwEdges[sampledId].v;
				Count deg_of_v = G.getDegree(vNode);

				seqNodesFromComponent.push_back(uNode);
				seqNodesFromComponent.push_back(vNode);

				seqNodeDegFromComponent.push_back(deg_of_u);
				seqNodeDegFromComponent.push_back(deg_of_v);
			}
			else
			{
				seqNodesFromComponent = nextLevelComponents[sampledId];
				for(int k = 0; k < seqNodesFromComponent.size(); ++k)
				{
					Count deg_of_node = G.getDegree(seqNodesFromComponent[k]);
					seqNodeDegFromComponent.push_back(deg_of_node);
				}
			}

			uniform_int_distribution<int> distNbor(0, seqNodeDegFromComponent[0] - 1);
			int rdNbrId = distNbor(gen);
			
			VertexIdx lastNode = G.getKthNeighbor(seqNodesFromComponent[0], rdNbrId);
			Count degOfLastNode = G.getDegree(lastNode);


			if((degOfLastNode > seqNodeDegFromComponent.back() || (degOfLastNode == seqNodeDegFromComponent.back() && lastNode > seqNodesFromComponent.back())))
			{

				int edgeCheckCount = 1;
				// starting k = 1, as at index 0 we have u and the current node is sampled from its neighborhood... 
				// and thus edgeCheckCount is 1 already...
				for(int k = 1; k < seqNodesFromComponent.size(); ++k)
				{
					bool edgeCheck = G.checkEdgeInAdjList(seqNodesFromComponent[k], lastNode);
					if(edgeCheck)
					{
						edgeCheckCount += 1;
						continue;
					}
					else
					{
						break;
					}
				}

				if(edgeCheckCount == seqNodesFromComponent.size())
				{
					Z = 1;
					seqNodesFromComponent.push_back(lastNode);	
					newNextLevelComponents.push_back(seqNodesFromComponent);
					newNextLevelDegrees.push_back(seqNodeDegFromComponent[0]);
					dRVals[i+1] += seqNodeDegFromComponent[0];
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

		// X = (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) *  Y;
		X = X * (numEdges / ljVals[2]);

		for(int j = 2; j <= i; j++)
		{
			X = X * (dRVals[j]/ljVals[j+1]);
		}

		X = X * Y;
		finalEstimate = X;
		allEstimates.push_back(finalEstimate);

		cout << i+1 << " Clique est = " << finalEstimate << endl;

		nextLevelDegrees.clear();
		nextLevelComponents.clear();

		nextLevelDegrees = newNextLevelDegrees;
		nextLevelComponents = newNextLevelComponents;

		newNextLevelDegrees.clear();
		newNextLevelComponents.clear();
	}

	struct prevSubGraphletSet lastSubGraphlets;
	lastSubGraphlets.graphletEstimate = finalEstimate;
	lastSubGraphlets.graphletsForNextLevel = nextLevelComponents;
	lastSubGraphlets.degreesForNextLevel = nextLevelDegrees;
	lastSubGraphlets.totalDegree = dRVals.back();
	lastSubGraphlets.allSizeCliqueEstimates = allEstimates;

	return lastSubGraphlets;
}