#include "include/namespace.h"
#include "include/graphIO.h"
#include "include/countCliques.h"
#include "include/count3Graphlets.h"
#include "include/count4Graphlets.h"
#include "include/count5Graphlets.h"
#include "include/utilities.h"
#include <iomanip>
#include <ctime>
#include <chrono>
#include <omp.h>


int main(int argc, char *argv[])
{
	if (argc != 5) 
	{
		std::cerr << "ERROR -- Requires 4 parameters: 1) Input Graph filename; 2) CodeWord for Graphlet To Count; 3) Graph name; and 4) Random walk edges(0) or U.A.R. edges (1)\n";
		std::cout << "Codeword for Graphlet to Count -- :\n" 
					<< "1) g32: Triangle (3-Clique)\n" 
					<< "2) g43: 4-Cycle\n"
					<< "3) g45: 4-Chord-Cycle\n" 
					<< "4) g46: 4-Clique\n" 
					<< "5) g58: 5-Clique-But-2-Edges-(Hatted-4-Clique)\n" 
					<< "6) g59: 5-Clique-But-1-Edge-(Almost-5-Clique)\n"
					<< "7) g510: 5-Clique\n"
					<< "8) g615: 6-Clique\n"
					<< "9) SRWg510: 5-Clique using SRW2\n"
					<< "Run code to count 3-Clique or Triangles using random walk edges: ./main-DiffStartPoint InputGraphFile g32 GraphName 0";  
		exit(1);
	}

	// sinaweibo G2Degree = 1611189745086
	// orkut G2Degree = 73948630564

	string inputFileName = argv[1];
	string whatCount = argv[2];
	string outFileGraphName = argv[3];
	int rwOrUar = atoi(argv[4]);

	double numEdgesInG2;

	if(outFileGraphName.compare("orkut") == 0 || outFileGraphName.compare("Orkut") == 0)
	{
		numEdgesInG2 = 73948630564.0;
	}
	else if(outFileGraphName.compare("sinaweibo") == 0 || outFileGraphName.compare("Sinaweibo") == 0)
	{
		numEdgesInG2 = 1611189745086.0;
	}


	cout << "Initializing graph... Populating Edges in Memory...\n";
	Graph G(inputFileName);				// constructor

	Count numVertices = G.getNumVertices();
	Count numEdges = G.getNumEdges();

	// G.printGraphDetails();
	cout << "n =  " << numVertices << endl;
	cout << "E =  " << numEdges << endl;

	int numRandomWalks = 100;
	
	ofstream outFile;

	vector<OrderedEdge> rWEdges;
	vector<OrderedEdge>::iterator rWEdgesIt;

	double avgVal = 0;
	cout.precision(20);

	std::chrono::steady_clock::time_point beginClock;
	std::chrono::steady_clock::time_point endClock;
	long long int perIterationTime = 0, totalTimePerEdgePerc = 0;

	srand(10000);
	// int randStartPoint = rand() % numVertices;
	cout << "before the start point -- \n";
	// VertexIdx randStartPoint = next() % numVertices;
	// VertexIdx randStartPoint = getRandomStartPoint(numVertices);
	// cout << "Got random start points -- "<< randStartPoint << "\n";

	// vector<float> percEdges = {5, 7, 10, 12, 15};
	// vector<double> percEdges = {0.1, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.50, 10.0};
	// vector<double> percEdges = {0.1, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
	vector<double> percEdges = {1.5, 2.0, 2.5, 3.0};
	// vector<double> percEdges = {3.0};

	vector<int> l3SamplesPerc = {80};
	vector<int> l4SamplesPerc = {100};
	vector<int> l5SamplesPerc = {100};
	vector<int> l6SamplesPerc = {100};


	// for(unsigned int l3i = 0; l3i < l3SamplesPerc.size(); l3i++)
	// {
	// int l3Perc = l3SamplesPerc[l3i];
	int l3Perc = l3SamplesPerc[0];
	// for(unsigned int l4i = 0; l4i < l4SamplesPerc.size(); l4i++)
	// {
	// int l4Perc = l4SamplesPerc[l4i];
	int l4Perc = l4SamplesPerc[0];
	int l5Perc = l5SamplesPerc[0];
	int l6Perc = l6SamplesPerc[0];

	string outFileName;

	if(rwOrUar == 0)
	{
		outFileName = "./allOutput/demet-" + outFileGraphName + "_" + whatCount + "_IntCliqueEstimates_RW.out";
		cout << "Going for random walks....\n";
	}
	else
	{
		outFileName = "./allOutput/demet-" + outFileGraphName + "_" + whatCount + "_IntCliqueEstimates_UAR.out";
		cout << "Getting UAR Edges....\n";
	}
	
	outFile.open(outFileName, ofstream::out);

	for(unsigned int j = 0; j < percEdges.size(); j++)
	{
		int lStep = percEdges[j] * (numEdges / 100);
		// outFile << "% Edges = " << percEdges[j] << " " << l3Perc << " " << l4Perc << " " << l5Perc << endl;
		outFile << "% Edges = " << percEdges[j] << endl;
		
		vector<double> allRunEsts;
		totalTimePerEdgePerc = 0;
		for(int k = 0; k < numRandomWalks; k++)
		{
			vector<double> allSizeEstimates;
			beginClock = chrono::steady_clock::now();
			
			if(rwOrUar == 0)
			{
				VertexIdx randStartPoint = getRandomStartPoint(numVertices);
				cout << "Got random start points -- "<< randStartPoint << "\n";

				// get edges from lstep random walk... (graphIO instance, graphIO.h, graphIO.cpp)
				rWEdges = G.getAllEdgesFromRStepRandomWalk(lStep, randStartPoint);
				cout << k << "th Random Walk -- Got random walk edges.... -- " << rWEdges.size() << "---" << percEdges[j] << "\n";
			}
			else
			{
				rWEdges = G.getUniformRandomEdges(lStep);
				cout << k << "th Random sample of edges.... -- " << rWEdges.size() << "---" << percEdges[j] << "\n";
			}

			double kGraphletCount = 0;
			struct prevSubGraphletSet xSubGraphlet;

			if(whatCount.compare("g32") == 0)
			{
				// rwCount3Graphlets C3;
				rwCountXCliqueGraphlets C3;
				// kGraphletCount =  C3.countTriangleGraphlet(G, rWEdges); // passing all seg_2's

				// xSubGraphlet = C3.countTriangleGraphlet(G, rWEdges, l3Perc);
				vector<int> percVals;
				percVals.push_back(l3Perc);
				xSubGraphlet = C3.countXCliqueGraphlet(G, rWEdges, percVals, 3);
				kGraphletCount = xSubGraphlet.graphletEstimate;
				allSizeEstimates = xSubGraphlet.allSizeCliqueEstimates;
			}
			else if(whatCount.compare("g43") == 0)
			{
				rwCount4Graphlets C4;
				kGraphletCount =  C4.count4Cycle(G, rWEdges, l3Perc, l4Perc);  // passing all seg_2's
			}
			else if(whatCount.compare("g45") == 0)
			{
				// cout << "I am coming here....\n";
				rwCount4Graphlets C4;
				kGraphletCount =  C4.count4ChordCycle(G, rWEdges, l3Perc, l4Perc);  // passing all seg_2's
				allSizeEstimates.push_back(kGraphletCount);
			}
			else if(whatCount.compare("g46") == 0)
			{
				// rwCount4Graphlets C4;
				// kGraphletCount =  C4.count4CliqueGraphlet(G, rWEdges, l3Perc, l4Perc);  // passing all seg_2's

				rwCountXCliqueGraphlets C4;

				vector<int> percVals;
				percVals.push_back(l3Perc);
				percVals.push_back(l4Perc);
				xSubGraphlet = C4.countXCliqueGraphlet(G, rWEdges, percVals, 4);
				kGraphletCount = xSubGraphlet.graphletEstimate;
				allSizeEstimates = xSubGraphlet.allSizeCliqueEstimates;
			}
			else if(whatCount.compare("g510") == 0)
			{
				// rwCount5Graphlets C5;
				// kGraphletCount =  C5.count5CliqueGraphlet(G, rWEdges, l3Perc, l4Perc, l5Perc);  // passing all seg_2's

				rwCountXCliqueGraphlets C5;

				vector<int> percVals;
				percVals.push_back(l3Perc);
				percVals.push_back(l4Perc);
				percVals.push_back(l5Perc);
				xSubGraphlet = C5.countXCliqueGraphlet(G, rWEdges, percVals, 5);
				kGraphletCount = xSubGraphlet.graphletEstimate;
				allSizeEstimates = xSubGraphlet.allSizeCliqueEstimates;
			}
			else if(whatCount.compare("g59") == 0)
			{
				rwCount5Graphlets C5;
				kGraphletCount =  C5.count5CliqueBut1Edge(G, rWEdges, l3Perc, l4Perc, l5Perc);  // passing all seg_2's
				allSizeEstimates.push_back(kGraphletCount);
			}
			else if(whatCount.compare("g615") == 0)
			{
				// rwCount5Graphlets C5;
				// kGraphletCount =  C5.count5CliqueGraphlet(G, rWEdges, l3Perc, l4Perc, l5Perc);  // passing all seg_2's

				rwCountXCliqueGraphlets C6;

				vector<int> percVals;
				percVals.push_back(l3Perc);
				percVals.push_back(l4Perc);
				percVals.push_back(l5Perc);
				percVals.push_back(l6Perc);
				xSubGraphlet = C6.countXCliqueGraphlet(G, rWEdges, percVals, 6);
				kGraphletCount = xSubGraphlet.graphletEstimate;
				allSizeEstimates = xSubGraphlet.allSizeCliqueEstimates;
			}
			else if(whatCount.compare("SRWg510") == 0)
			{
				rwCount5Graphlets C5;
				int numSteps = lStep;
				kGraphletCount =  C5.SRWCount5CliqueGraphlet(G, numSteps, numEdgesInG2);  // passing all seg_2's
				allSizeEstimates.push_back(kGraphletCount);
			} 

			endClock = chrono::steady_clock::now();

			perIterationTime = chrono::duration_cast<std::chrono::microseconds> (endClock - beginClock).count();
			totalTimePerEdgePerc += perIterationTime;

			allRunEsts.push_back(kGraphletCount);
			avgVal = avgVal + (kGraphletCount - avgVal)/(k+1);
			// cout << "Running avg = " << avgVal << endl;
			
			// cout << k << "-th Estimate = " << kGraphletCount << endl;
			// cout << k << " " << percEdges[j] << " " << setprecision(20) << kGraphletCount << " " << setprecision(20) <<  perIterationTime << endl;

			cout << k << " " << percEdges[j] << " ";

			if(allSizeEstimates.size() > 1)
			{
				for(int q = 2; q < allSizeEstimates.size(); q++)
					cout << setprecision(20) << allSizeEstimates[q] << " ";
			}
			else
			{
				cout << setprecision(20) << allSizeEstimates[0] << " ";
			}

			cout << setprecision(20) <<  perIterationTime << endl;

			// outFile << k << " " << percEdges[j] << " " << setprecision(20) << kGraphletCount << " " << setprecision(20) <<  perIterationTime << endl;

			outFile << k << " " << percEdges[j] << " ";
			if(allSizeEstimates.size() > 1)
			{
				for(int q = 2; q < allSizeEstimates.size(); q++)
					outFile << setprecision(20) << allSizeEstimates[q] << " ";
			}
			else
			{
				outFile << setprecision(20) << allSizeEstimates[0] << " ";
			}

			outFile << setprecision(20) <<  perIterationTime << endl;

			outFile.flush();
		}
	}

	outFile.close();
		// }
	// }
	return 0;
}
