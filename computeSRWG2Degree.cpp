#include "include/namespace.h"
#include "include/graphIO.h"
#include "include/utilities.h"
#include <iomanip>

using namespace std;

int main(int argc, char *argv[])
{
	string inputFileName = argv[1];
	Graph G(inputFileName);
	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	cout << "Graph loaded in the memory... " << numVertices << " " << numEdges << endl;

	double G2Degree = 0.0;

	vector <vector<VertexIdx>>::iterator edgeListIt;

	for(edgeListIt = G.edgeList.begin(); edgeListIt != G.edgeList.end(); ++edgeListIt)
	{
		vector<VertexIdx> eachEdge = *edgeListIt;

		Count deg_of_u = G.getDegree(eachEdge[0]);
		Count deg_of_v = G.getDegree(eachEdge[1]);
		Count eachEdgeDegree =  deg_of_u + deg_of_v;
		G2Degree = G2Degree + eachEdgeDegree;
	}


	cout << "G2 Degree = " << setprecision(20) << G2Degree << endl;

	return G2Degree;
}