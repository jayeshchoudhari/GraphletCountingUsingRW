#ifndef COUNT_5_GRAPHLETS_H 
#define COUNT_5_GRAPHLETS_H 

class rwCount5Graphlets
{
    public:
        double count5CliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc, int l5Perc);
        double count5CliqueBut1Edge(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc, int l5Perc);
        double count5CliqueBut2Edges(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc, int l5Perc);


        double SRWCount5CliqueGraphlet(Graph &G, int numSteps);
};

#endif      //COUNT_5_GRAPHLETS_H