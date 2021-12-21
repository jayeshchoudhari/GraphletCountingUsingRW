#include "namespace.h"
#include "graphIO.h"
#include "countCliques.h"
#include "count5Graphlets.h"

using namespace std;


double rwCount5Graphlets :: count5CliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc, int l5Perc) 
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

	vector<vector<VertexIdx>> nextLevelComponents, fourthLevelComponents;
	vector<VertexIdx> nextLevelDegrees, fourthLevelDegrees;
	double dR3 = 0.0, dR4 = 0.0;

    // int subsample_size = l2/20;
    int subsample_size = l2 * l3Perc/100.0;
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
	
	int l3_subsample_size = l3 * l4Perc/100.0;
    ljVals.push_back(l3_subsample_size);

	discrete_distribution<int> nextLevelDist (nextLevelDegrees.begin(), nextLevelDegrees.end());
	// omp_set_num_threads(1);
	// #pragma omp parallel for
    for(int s = 0; s < l3_subsample_size; s++)
	{
        int sampledId = nextLevelDist(gen);

        VertexIdx uNode = nextLevelComponents[sampledId][0];
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
        	VertexIdx vNode = nextLevelComponents[sampledId][1];
		    // check if x is connected to every other vertex...
            bool vxEdge = G.checkEdgeInAdjList(vNode, xNode);
            if(vxEdge)
            {
                bool wxEdge = G.checkEdgeInAdjList(wNode, xNode);
                if(wxEdge)
                {
                    finalInd = 1;
                    vector<VertexIdx> tempComponent {uNode, vNode, wNode, xNode};
					fourthLevelComponents.push_back(tempComponent);
					fourthLevelDegrees.push_back(deg_of_u);
					dR4 += deg_of_u;
                }
            }
        }

		Y += finalInd;
	}

	cout << "samples found -- " << Y << endl;

	X = (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) *  Y;
	cout << "4 Clique est = " << X << endl;

	Y = 0;
	X = 0;

	Count l4 = fourthLevelComponents.size();
	cout << "size of 4-cliques -- " << l4 << endl;
	
	int l4_subsample_size = l4 * l5Perc/100.0;
    ljVals.push_back(l4_subsample_size);

	discrete_distribution<int> fourthLevelDist (fourthLevelDegrees.begin(), fourthLevelDegrees.end());

	for(int s = 0; s < l4_subsample_size; s++)
	{
        int sampledId = fourthLevelDist(gen);

        VertexIdx uNode = fourthLevelComponents[sampledId][0];
        VertexIdx xNode = fourthLevelComponents[sampledId][3];

		Count deg_of_u = fourthLevelDegrees[sampledId];
		Count deg_of_x = G.getDegree(xNode);

		uniform_int_distribution<int> distNbor(0, deg_of_u - 1);
        int rdNbrId = distNbor(gen);

        VertexIdx yNode = G.getKthNeighbor(uNode, rdNbrId);
		Count deg_of_y = G.getDegree(yNode);

		bool xdegreeOrder = (deg_of_y > deg_of_x) || (deg_of_y == deg_of_x && yNode > xNode);

		bool finalInd = 0;
        if(xdegreeOrder)
        {
			VertexIdx vNode = fourthLevelComponents[sampledId][1];
        	VertexIdx wNode = fourthLevelComponents[sampledId][2];
		    // check if y is connected to every other vertex...
            bool vyEdge = G.checkEdgeInAdjList(vNode, yNode);
            if(vyEdge)
            {
                bool wyEdge = G.checkEdgeInAdjList(wNode, yNode);
                if(wyEdge)
                {
                	bool xyEdge = G.checkEdgeInAdjList(xNode, yNode);
	                if(xyEdge)
	                {
                    	finalInd = 1;
                        // cout << "Hit... \t";
                    }
                }
            }
        }

		Y += finalInd;
	}

	cout << "samples found -- " << Y << endl;
	cout <<  ljVals[2] << " " <<  dR2 << " " << ljVals[3]  << " " <<  dR3 << " " << ljVals[4] << " " << dR4 << " " << ljVals[5] << " " << Y << "\n";
	X = (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) * (dR4/ljVals[5]) * Y;
	cout << "5 Clique est = " << X << endl;

    return X;
}


double rwCount5Graphlets :: count5CliqueBut1Edge(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc, int l5Perc) 
{

    rwCountXCliqueGraphlets C4;
    struct prevSubGraphletSet x4SubGraphlet;

    vector<int> percVals;
    percVals.push_back(l3Perc);
    percVals.push_back(l4Perc);
    x4SubGraphlet = C4.countXCliqueGraphlet(G, rwEdges, percVals, 4);

    vector<vector<VertexIdx>> nextLevelComponents;
    vector<vector<VertexIdx>> nextLevelNeighbors;
    vector<vector<VertexIdx>> nextLevelAllPivotSizes;
    vector<VertexIdx> nextLevelDegrees;
    double degreeSum = 0.0;

    nextLevelComponents = x4SubGraphlet.graphletsForNextLevel;
    int maxCliqueSize = x4SubGraphlet.allSizeCliqueEstimates.size();
    double clique4Est = x4SubGraphlet.allSizeCliqueEstimates[maxCliqueSize-1];

    for(int s = 0; s < nextLevelComponents.size(); s++)
    {
       vector<VertexIdx> currComponent = nextLevelComponents[s];
       
        struct pivotNeighborsAndSizes_X_2 pivotNeighborsAndSizes = G.getDegreeAndNeighborsOf2Qset(currComponent, 2);
        vector<VertexIdx> pivotNbrs = pivotNeighborsAndSizes.pivotNeighbors;

        nextLevelNeighbors.push_back(pivotNbrs);
        nextLevelDegrees.push_back(pivotNbrs.size());
        degreeSum += pivotNbrs.size();
    }

    Count l4 = nextLevelComponents.size();
    int subsample_size = l4 * l5Perc/100.0;
    double ljVals_5 = subsample_size;

    double prevY = l4;
    double X = 0.0, Y = 0.0, Z = 0.0;

    random_device rd;
    mt19937 gen(rd());

    discrete_distribution<int> nextLevelDist(nextLevelDegrees.begin(), nextLevelDegrees.end());

    for(int s = 0; s < subsample_size; s++)
    {
        int sampledId = nextLevelDist(gen);
        vector<VertexIdx> seqNodesFromComponent;
        vector<Count> seqNodeDegFromComponent;

        seqNodesFromComponent = nextLevelComponents[sampledId];

        Count componentDeg = nextLevelDegrees[sampledId];
        uniform_int_distribution<int> distNbor(0, componentDeg - 1);
        int rdNbrId = distNbor(gen);

        VertexIdx lastNode = nextLevelNeighbors[sampledId][rdNbrId];
        seqNodesFromComponent.push_back(lastNode);

        vector<VertexIdx> lastNotConnnectedTo = G.checkLastNodeNotConnectedTo(seqNodesFromComponent);

        bool finalInd = 0;

        if(lastNotConnnectedTo.size() == 1)
        {
            VertexIdx checkNode = lastNotConnnectedTo[0];
            Count deg_of_checkNode = G.getDegree(checkNode);
            Count deg_of_lastNode = G.getDegree(lastNode);

            // if((deg_of_checkNode <= deg_of_x) && (deg_of_w <= deg_of_x))
            if((deg_of_checkNode < deg_of_lastNode) || (deg_of_checkNode == deg_of_lastNode && checkNode < lastNode))
            {
                finalInd = 1;
            }
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
    // (numEdges / ljVals[2]) * (dR2/ljVals[3]) * (dR3/ljVals[4]) * (dR4/ljVals[5]) * Y;
    X = (clique4Est / prevY) * (degreeSum/ljVals_5) * Y;
    cout << clique4Est << " " << prevY << " " << degreeSum << " " << ljVals_5 << " " << Y << endl; 
    cout << "5 Clique -1 est = " << X << endl;
    return X;
}