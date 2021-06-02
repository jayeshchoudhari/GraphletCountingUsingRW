#include "include/namespace.h"
#include "include/graphIO.h"
#include "include/count3Graphlets.h"
#include "include/count4Graphlets.h"
#include "include/utilities.h"
// #include <bits/stdc++.h>
#include <iomanip>
#include <ctime>
#include <chrono>
#include <omp.h>


int main(int argc, char *argv[])
{
    string inputFileName = argv[1];

    cout << "initializing graph...\n";
	Graph G(inputFileName);				// constructor
	// Graph G(numV, e);
	// e.clear();

	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	// G.printGraphDetails();
    cout << "n =  " << numVertices << endl;
    cout << "E =  " << numEdges << endl;

    // int maxK;
	int numRandomWalks = 10;
	
    int lStep = stoi(argv[2]);      /// not using this argument currently

    string outFileName = argv[3];
	// lStep = numEdges / 20;

    ofstream outFile;
    outFile.open(outFileName, ofstream::out);

    // for(int i = 0; i < 1; ++i)
    // {
    vector<OrderedEdge> rWEdges;
    vector<OrderedEdge>::iterator rWEdgesIt;

    double avgVal = 0;
    cout.precision(20);
    std::chrono::steady_clock::time_point beginClock;
    std::chrono::steady_clock::time_point endClock;
    long long int perIterationTime = 0, totalTimePerEdgePerc = 0;

    // cout << "No floating point exceptions here\n";

    srand(10);
    // int randStartPoint = rand() % numVertices;
    cout << "before the start point -- \n";
    VertexIdx randStartPoint = next() % numVertices;
    cout << "Got start point -- "<< randStartPoint << "\n";

    rwCount3Graphlets C3;
    // rwCount4Graphlets C4;

    // vector<int> percEdges = {1, 5, 7, 10};
    // vector<double> percEdges = {0.1, 0.3, 0.5, 0.7};
    // vector<double> percEdges = {1.5, 2.0, 2.5, 3.0};
    vector<double> percEdges = {3.0};

    cout << "Going for random walks....\n";
    for(unsigned int j = 0; j < percEdges.size(); j++)
    {
        lStep = percEdges[j] * (numEdges / 100);
        outFile << "% Edges = " << percEdges[j] << endl;
        
        vector<double> allRunEsts;
        totalTimePerEdgePerc = 0;
        for(int k = 0; k < numRandomWalks; k++)
        {
            beginClock = chrono::steady_clock::now();
            rWEdges = G.getAllEdgesFromRStepRandomWalk(lStep, randStartPoint);
            cout << k << "th Random Walk -- Got random walk edges.... -- " << rWEdges.size() << "---" << percEdges[j] << "\n";
            
            double kGraphletCount =  C3.countTriangleGraphlet(G, rWEdges);	// passing all seg_2's
            // double kGraphletCount =  C4.count4CliqueGraphlet(G, rWEdges);	// passing all seg_2's
            // double kGraphletCount =  C4.count4ChordCycle(G, rWEdges);	// passing all seg_2's
            endClock = chrono::steady_clock::now();
            perIterationTime = chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
            totalTimePerEdgePerc += perIterationTime;
            cout << k << "-th Estimate = " << kGraphletCount << endl;
            outFile << setprecision(20) << kGraphletCount << " " << setprecision(20) <<  perIterationTime << endl;
            
            allRunEsts.push_back(kGraphletCount);
            avgVal = avgVal + (kGraphletCount - avgVal)/(k+1);

            cout << "Running avg = " << avgVal << endl;
            cout << "Time per iteration = " << perIterationTime << endl;
        }
    }
    // }

    return 0;
}
