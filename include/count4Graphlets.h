#ifndef COUNT_4_GRAPHLETS_H 
#define COUNT_4_GRAPHLETS_H 

class rwCount4Graphlets
{
    public:
        double count4CliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc);
        double count4ChordCycle(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc);
        double count4Cycle(Graph &G, vector<OrderedEdge> rwEdges, int l3Perc, int l4Perc);

        // double count4Cycle(vector<OrderedEdge> rwEdges);
        // double count4CycleNew(vector<OrderedEdge> rwEdges);

};

#endif      //COUNT_4_GRAPHLETS_H