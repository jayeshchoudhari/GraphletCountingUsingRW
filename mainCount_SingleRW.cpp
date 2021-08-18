#include "include/namespace.h"
#include "include/graphIO.h"
#include "include/count3Graphlets.h"
#include "include/count4Graphlets.h"
#include "include/utilities.h"
#include <iomanip>
#include <ctime>
#include <chrono>
#include <omp.h>


int main(int argc, char *argv[])
{
    if (argc != 4) 
	{
        std::cerr << "ERROR -- Requires 3 parameters: 1) Input Graph filename; 2) CodeWord for Graphlet To Count; and 3) Output filename;\n";
        std::cout << "Codeword for Graphlet to Count -- :\n" 
                    << "1) g32: Triangle (3-Clique)\n" 
                    << "2) g43: 4-Cycle\n"
                    << "3) g45: 4-Chord-Cycle\n" 
                    << "4) g46: 4-Clique\n" 
                    << "Run code to count 3-Clique or Triangles: ./main-organized-make InputGraphFile g32 OutputFile";  
        exit(1);
    }

    string inputFileName = argv[1];
    string whatCount = argv[2];
    string outFileName = argv[3];

    cout << "Initializing graph... Populating Edges in Memory...\n";
	Graph G(inputFileName);				// constructor

	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	// G.printGraphDetails();
    cout << "n =  " << numVertices << endl;
    cout << "E =  " << numEdges << endl;

	int numRandomWalks = 100;
	
    ofstream outFile;
    outFile.open(outFileName, ofstream::out);

    // if(whatCount.find("g3") != std::string::npos)
    //     rwCount3Graphlets C3;
    // else if(whatCount.find("g4") != std::string::npos)
    //     rwCount4Graphlets C4;

    vector<OrderedEdge> rWEdges;
    vector<OrderedEdge>::iterator rWEdgesIt;

    double avgVal = 0;
    cout.precision(20);

    std::chrono::steady_clock::time_point beginClock;
    std::chrono::steady_clock::time_point endClock;
    long long int perIterationTime = 0, totalTimePerEdgePerc = 0;

    // srand(10000);
    // int randStartPoint = rand() % numVertices;
    cout << "before the start point -- \n";
    VertexIdx randStartPoint = getRandomStartPoint(numVertices);
    cout << "Got start point -- "<< randStartPoint << "\n";

    // vector<int> percEdges = {1, 5, 7, 10};
    // vector<double> percEdges = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0};
    vector<float> percEdges = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
    // vector<double> percEdges = {3.0};

    float maxPercentageValue = 10.0;
    int maxPercentageEdges = (maxPercentageValue * (numEdges / 100)) + 1;

    vector<vector<OrderedEdge>> allIncrementalRWEdges;
    
    for(int k = 0; k < numRandomWalks; k++)
    {
        rWEdges = G.getAllEdgesFromRStepRandomWalk(maxPercentageEdges, randStartPoint);
        allIncrementalRWEdges.push_back(rWEdges);
        cout << k << "-th Random walk -- Got random walk edges for k-th random walk" << rWEdges.size()  << "\n";
    }

    cout << "Got all the random walks each of max size...\n";

    for(unsigned int j = 0; j < percEdges.size(); j++)
    {
        int lStep = percEdges[j] * (numEdges / 100);
        outFile << "% Edges = " << percEdges[j] << endl;
        
        vector<double> allRunEsts;
        totalTimePerEdgePerc = 0;
        for(int k = 0; k < numRandomWalks; k++)
        {
            beginClock = chrono::steady_clock::now();
            // rWEdges = G.getAllEdgesFromRStepRandomWalk(lStep, randStartPoint);
            vector<OrderedEdge> slicedRwEdges(&allIncrementalRWEdges[k][0], &allIncrementalRWEdges[k][lStep]);
            cout << k << "th Random Walk -- Got random walk edges.... -- " << slicedRwEdges.size() << "---" << percEdges[j] << "\n";
            
            double kGraphletCount = 0;

            if(whatCount.compare("g32") == 0)
            {
                rwCount3Graphlets C3;
                kGraphletCount =  C3.countTriangleGraphlet(G, slicedRwEdges);	// passing all seg_2's
            }
            else if(whatCount.compare("g43") == 0)
            {
                rwCount4Graphlets C4;
                // kGraphletCount =  C4.count4Cycle(G, slicedRwEdges);	// passing all seg_2's
                std::cout << "To be added...\n";
            }
            else if(whatCount.compare("g45") == 0)
            {
                rwCount4Graphlets C4;
                kGraphletCount =  C4.count4ChordCycle(G, slicedRwEdges);	// passing all seg_2's
            }
            else if(whatCount.compare("g46") == 0)
            {
                rwCount4Graphlets C4;
                kGraphletCount =  C4.count4CliqueGraphlet(G, slicedRwEdges, 5, 20);	// passing all seg_2's
            }
            endClock = chrono::steady_clock::now();

            perIterationTime = chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
            totalTimePerEdgePerc += perIterationTime;

            allRunEsts.push_back(kGraphletCount);
            avgVal = avgVal + (kGraphletCount - avgVal)/(k+1);
            // cout << "Running avg = " << avgVal << endl;
            
            // cout << k << "-th Estimate = " << kGraphletCount << endl;
            cout << k << " " << percEdges[j] << " " << setprecision(20) << kGraphletCount << " " << setprecision(20) <<  perIterationTime << endl;
            outFile << k << " " << percEdges[j] << " " << setprecision(20) << kGraphletCount << " " << setprecision(20) <<  perIterationTime << endl;
        }
    }

    outFile.close();
    return 0;
}
