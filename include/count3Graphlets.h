#ifndef COUNT_3_GRAPHLETS_H 
#define COUNT_3_GRAPHLETS_H 

class rwCount3Graphlets
{
    public:
        // double countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges);
        prevSubGraphletSet countTriangleGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc);
};

#endif      // COUNT_3_GRAPHLETS_H