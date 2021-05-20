class rwCount3Graphlets
{
    public:
        double countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int maxGraphLetSize);
};

double rwCount3Graphlets :: countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int maxGraphLetSize)
{
    vector<double> dRVals = {0, 0};
	vector<double> ljVals = {0, 0, rwEdges.size() * 1.0};

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

    int subsample_size = l2/50;
    ljVals.push_back(subsample_size);

    double X = 0, Y = 0, Z = 0;

    // omp_set_num_threads(2);
	// 18060671
    // #pragma omp parallel for
    for(int s = 0; s < subsample_size; s++)
	{
        discrete_distribution<int> distribution (edge_degree_list.begin(), edge_degree_list.end());
        int sampledId = distribution(gen);

        VertexIdx uNode = rwEdges[sampledId].u;
        VertexIdx vNode = rwEdges[sampledId].v;

        uniform_int_distribution<int> distNbor(0, rwEdges[sampledId].degree - 1);
        int rdNbrId = distNbor(gen);
        VertexIdx wNode = G.getKthNeighbor(uNode, rdNbrId);

        Count deg_of_v = G.getDegree(vNode);
        Count deg_of_w = G.getDegree(wNode);
        bool vwEdge = G.checkEdgeInAdjList(vNode, wNode);

        if(vwEdge && (deg_of_w > deg_of_v || (deg_of_w == deg_of_v && wNode > vNode)))
        {
            Z = 1;
        }
        else
        {
            Z = 0;
        }
        Y += Z;
		if(s % 1000 == 0)
		{
			cout << "samples seen = " << s << endl;
		}
    }
    // X = Y / subsample_size;

	X =  (G.getNumEdges() / ljVals[2]) * (dR2/ljVals[3]) * Y;
    return X;
}