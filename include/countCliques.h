#ifndef COUNT_X_CLIQUE_GRAPHLETS_H 
#define COUNT_X_CLIQUE_GRAPHLETS_H

class rwCountXCliqueGraphlets
{
    public:
        prevSubGraphletSet countXCliqueGraphlet(Graph &G, vector<OrderedEdge> rwEdges, vector<int> percSamples, int k);
};

#endif      //COUNT_X_CLIQUE_GRAPHLETS_H