#include "namespace.h"
#include "graphIO.h"
#include "count4Graphlets.h"
#include <bits/stdc++.h>

using namespace std;

double rwCount4Graphlets :: count4CliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges) 
{
    vector<double> dRVals = {0, 0};
	vector<double> ljVals = {0, 0, rwEdges.size() * 1.0};
    Count numEdges = G.getNumEdges();

	vector<ePair> graphletEdges;
	// int k = maxGraphLetSize;

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

    int subsample_size = l2/20;
    ljVals.push_back(subsample_size);

    double X = 0, Y = 0, Z = 0;

    discrete_distribution<int> distribution (edge_degree_list.begin(), edge_degree_list.end());
    
    // omp_set_num_threads(1);
	// 18060671
    // #pragma omp parallel for
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
        // if(vwEdge && (deg_of_w > deg_of_u || (deg_of_w == deg_of_u && wNode > uNode)))
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
    // X = Y / subsample_size;

	X =  (numEdges / ljVals[2]) * (dR2/ljVals[3]) * Y;

	cout << "Triangle Est = " << X << endl; 

	Y = 0;
	X = 0;
	Count l3 = nextLevelComponents.size();

	cout << "size of triangles -- " << l3 << endl;
	int l3_subsample_size = l3/5;
    ljVals.push_back(l3_subsample_size);

	discrete_distribution<int> nextLevelDist (nextLevelDegrees.begin(), nextLevelDegrees.end());
	// omp_set_num_threads(1);
	// #pragma omp parallel for
    for(int s = 0; s < l3_subsample_size; s++)
	{
        int sampledId = nextLevelDist(gen);

        VertexIdx uNode = nextLevelComponents[sampledId][0];
        VertexIdx vNode = nextLevelComponents[sampledId][1];
        VertexIdx wNode = nextLevelComponents[sampledId][2];

		Count deg_of_u = nextLevelDegrees[sampledId];
		// Count deg_of_v = G.getDegree(vNode);
		Count deg_of_w = G.getDegree(wNode);

		uniform_int_distribution<int> distNbor(0, deg_of_u - 1);
        int rdNbrId = distNbor(gen);

        VertexIdx xNode = G.getKthNeighbor(uNode, rdNbrId);
		Count deg_of_x = G.getDegree(xNode);

		bool wdegreeOrder = (deg_of_x > deg_of_w) || (deg_of_x == deg_of_w && xNode > wNode);

		bool finalInd = 0;
        if(wdegreeOrder)
        {
		    // check if x is connected to every other vertex...
            bool vxEdge = G.checkEdgeInAdjList(vNode, xNode);
            if(vxEdge)
            {
                bool wxEdge = G.checkEdgeInAdjList(wNode, xNode);
                if(wxEdge)
                {
                    finalInd = 1;
                }
            }
        }

		Y += finalInd;
	}

	cout << "samples found -- " << Y << endl;

	X = (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) *  Y;
	cout << "4 Clique est = " << X << endl;

    return X;
}



double rwCount4Graphlets :: count4ChordCycle(Graph &G, vector<OrderedEdge> rwEdges)
{
	vector<double> dRVals = {0, 0};
	vector<double> ljVals = {0, 0, rwEdges.size() * 1.0};
    Count numEdges = G.getNumEdges();

	vector<ePair> graphletEdges;
	// int k = maxGraphLetSize;

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
	cout << "Computed total degree at l2\n";

	vector<vector<VertexIdx>> nextLevelComponents;
	vector<vector<VertexIdx>> nextLevelNeighbors;
	vector<vector<VertexIdx>> nextLevelAllPivotSizes;
	vector<VertexIdx> nextLevelDegrees;
	double dR3 = 0.0;

    int subsample_size = l2/10;
    ljVals.push_back(subsample_size);
	// cout << "subsample_size = " << subsample_size << "\n";
	// cout << "edge_degree_list_size = " << edge_degree_list.size() << "\n";

    double X = 0, Y = 0, Z = 0;

	discrete_distribution<int> distribution (edge_degree_list.begin(), edge_degree_list.end());
    // omp_set_num_threads(1);
	// 18060671
    // #pragma omp parallel for
    for(int s = 0; s < subsample_size; s++)
	{
		int sampledId = distribution(gen);
		
		VertexIdx uNode = rwEdges[sampledId].u;
		Count deg_of_u = rwEdges[sampledId].degree;
		uniform_int_distribution<int> distNbor(0, deg_of_u - 1);
		int rdNbrId = distNbor(gen);
		VertexIdx wNode = G.getKthNeighbor(uNode, rdNbrId);
		// cout << "Got the next node w = "<< wNode << "--deg of u = " << deg_of_u << "...\n";

		VertexIdx vNode = rwEdges[sampledId].v;
		Count deg_of_v = G.getDegree(vNode);
		Count deg_of_w = G.getDegree(wNode);

        // Degree order check -- assignment rule for triangle...
        if((deg_of_w > deg_of_v || (deg_of_w == deg_of_v && wNode > vNode)))
        // if(vwEdge && (deg_of_w > deg_of_u || (deg_of_w == deg_of_u && wNode > uNode)))
        {   
            // Check for triangle...
        	bool vwEdge = G.checkEdgeInAdjList(vNode, wNode);
			if(vwEdge)
			{
				Z = 1;
				vector<VertexIdx> tempComponent {uNode, vNode, wNode};
				nextLevelComponents.push_back(tempComponent);

				// Get  pivot set of (2 out of 3)nodes and their neighbors... 
				struct pivotNeighborsAndSizes_3_2 pivotNeighborsAndSizes = G.getDegreeAndNeighborsOf2Qset(tempComponent, 2);
				vector<VertexIdx> pivotNbrs = pivotNeighborsAndSizes.pivotNeighbors;

				nextLevelNeighbors.push_back(pivotNbrs);
				nextLevelDegrees.push_back(pivotNbrs.size());
				dR3 += pivotNbrs.size();
				
				nextLevelAllPivotSizes.push_back(pivotNeighborsAndSizes.neighborSizes);
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
		// cout << "Got the indicator for triangle -- " << "\n";
        Y += Z;
    }
    // X = Y / subsample_size;
 
	X =  (numEdges / ljVals[2]) * (dR2/ljVals[3]) * Y;
	cout << "Triangle Est = " << X << endl; 


	Y = 0;
	X = 0;
	// Count l3 = nextLevelComponents.size();
	Count l3 = nextLevelNeighbors.size();

	cout << "size of triangles -- " << l3 << endl;
	int l3_subsample_size = l3/5;
    ljVals.push_back(l3_subsample_size);

	discrete_distribution<int> nextLevelDist (nextLevelDegrees.begin(), nextLevelDegrees.end());
	// omp_set_num_threads(1);
	// #pragma omp parallel for
	for(int s = 0; s < l3_subsample_size; s++)
	{
		int sampledId = nextLevelDist(gen);

		VertexIdx uNode = nextLevelComponents[sampledId][0];
		VertexIdx vNode = nextLevelComponents[sampledId][1];
		VertexIdx wNode = nextLevelComponents[sampledId][2];

		Count componentDeg = nextLevelDegrees[sampledId];
		uniform_int_distribution<int> distNbor(0, componentDeg - 1);
		int rdNbrId = distNbor(gen);

		VertexIdx xNode = nextLevelNeighbors[sampledId][rdNbrId];

		bool finalInd = 0;
		bool connectionCheck = G.checkConnectionOfXToAny2OfUVW(uNode, vNode, wNode, xNode);

		// degree order check...
		if(connectionCheck)
		{
			// check if assignment to seg_{k-1}
			Count deg_of_uv = nextLevelAllPivotSizes[sampledId][0];
			Count deg_of_uw = nextLevelAllPivotSizes[sampledId][1];
			Count deg_of_vw = nextLevelAllPivotSizes[sampledId][2];

			Count minDeg_uvw =  deg_of_uv <= deg_of_uw ? ((deg_of_uv <= deg_of_vw) ? deg_of_uv : deg_of_vw) : ((deg_of_uw <= deg_of_vw) ? deg_of_uw : deg_of_vw);
			// (one < ((two < three) ? two:three)) ? one:((two < three) ? two:three)

			Count deg_of_ux = G.getCombinedNeighborSize(uNode, xNode);
			Count deg_of_vx = G.getCombinedNeighborSize(vNode, xNode);
			Count deg_of_wx = G.getCombinedNeighborSize(wNode, xNode);

			// bool checkuv = (deg_of_uv <= deg_of_ux) && (deg_of_uv <= deg_of_vx) && (deg_of_uv <= deg_of_wx);
			// bool checkuw = (deg_of_uw <= deg_of_ux) && (deg_of_uw <= deg_of_vx) && (deg_of_uw <= deg_of_wx);
			// bool checkvw = (deg_of_vw <= deg_of_ux) && (deg_of_vw <= deg_of_vx) && (deg_of_vw <= deg_of_wx);

            // check if the assignment is to {k-1}-subgraphlet... 
			bool checkuvw = (minDeg_uvw <= deg_of_ux) && (minDeg_uvw <= deg_of_vx) && (minDeg_uvw <= deg_of_wx);

			finalInd = checkuvw;
		}
		
		if(finalInd)
		{
			Z = 1;
		}
		else
		{
			Z = 0;
		}
		Y += Z;
	}

	cout << "samples found -- " << Y << endl;

	X = (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) *  Y;
	cout << "4 Chord Cycle Est = " << X << endl;

    return X;
}
