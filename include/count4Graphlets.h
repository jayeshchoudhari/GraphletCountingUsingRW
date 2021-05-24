class rwCount4Graphlets
{
    public:
        double count4CliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges);

        // double count4ChordCycle(vector<OrderedEdge> rwEdges);
        // double count4ChordCycleNew(vector<OrderedEdge> rwEdges);

        // double count4Cycle(vector<OrderedEdge> rwEdges);
        // double count4CycleNew(vector<OrderedEdge> rwEdges);
};

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

    int subsample_size = l2/50;
    ljVals.push_back(subsample_size);

    double X = 0, Y = 0, Z = 0;

    // omp_set_num_threads(1);
	// 18060671
    // #pragma omp parallel for
    for(int s = 0; s < subsample_size; s++)
	{
        discrete_distribution<int> distribution (edge_degree_list.begin(), edge_degree_list.end());
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
				vector<VertexIdx> tempComponent;
				tempComponent.push_back(uNode);
				tempComponent.push_back(vNode);
				tempComponent.push_back(wNode);
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
	int l3_subsample_size = l3;
    ljVals.push_back(l3_subsample_size);

	// omp_set_num_threads(1);
	// #pragma omp parallel for
    for(int s = 0; s < l3_subsample_size; s++)
	{
		discrete_distribution<int> distribution (nextLevelDegrees.begin(), nextLevelDegrees.end());
        int sampledId = distribution(gen);

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