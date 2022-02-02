#include "namespace.h"
#include "utilities.h"
#include <random>

using namespace std;

VertexIdx getRandomStartPoint(Count numVertices)
{
	std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> distVert(0, numVertices-1); // distribution in range [1, 6]

	VertexIdx randVertId = distVert(rng);
    // std::cout << distVert(rng) << std::endl;
	return randVertId;
}

int notInList(vector<VertexIdx> initialRandomWalkNodes, VertexIdx nextNode)
{
	for(int i = 0; i < initialRandomWalkNodes.size(); i++)
	{
		if(nextNode == initialRandomWalkNodes[i])
			return 0;
	}
	return 1;
}


uint64_t next()
{
    uint64_t shuffle_table[4] = {23673, 34793, 5690, 4673};
    uint64_t s1 = shuffle_table[0];
    uint64_t s0 = shuffle_table[1];
    uint64_t result = s0 + s1;
    shuffle_table[0] = s0;
    s1 ^= s1 << 23;
    shuffle_table[1] = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5);
    return result;
}


// Read file and get number of vertices and all the edges as a vector of vectors...
pair<Count, vector<ePair>> getSetOfEdges(string graphFileName)
{
	vector<ePair> allEdges;

	ifstream graphFile;

	int lineNum = 0;
	Count numVertices, numEdges, dummyVal;
	VertexIdx u, v;

	stringstream ss;
	string line;

	graphFile.open(graphFileName, ifstream::in);
	if(graphFile.is_open())
	{
		if (getline(graphFile, line))
		{
			ss.clear();
			ss.str("");
			ss << line;
			ss >> numVertices >> dummyVal >> numEdges;

			while(getline(graphFile, line))
			{
				ss.clear();
				ss.str("");
				ss << line;
				ss >> u >> v;

				allEdges.push_back(make_pair(u, v));
				lineNum += 1;
			}
		}
	}
	else
	{
		cout << "Cannot open the said file...\n";
		exit(0);
	}

	return make_pair(numVertices, allEdges);
}

vector<VertexIdx> getUnionOfNeighbors(vector<VertexIdx> set1, vector<VertexIdx> set2)
{
	vector<VertexIdx> combinedSet;
	set_union(set1.begin(),set1.end(), set2.begin(),set2.end(),back_inserter(combinedSet));
	return combinedSet;
}

vector<VertexIdx> getIntersectionOfTwoVectors(vector<VertexIdx> set1, vector<VertexIdx> set2)
{
	vector<VertexIdx> set3;
    set_intersection(set1.begin(),set1.end(), set2.begin(),set2.end(),back_inserter(set3));
	return set3;
}

int sampleFromDiscreteDistOfR_prev(vector<i64> weightVector)	
{	
  	// default_random_engine generator;
	random_device rd;
    mt19937 gen(rd());
	discrete_distribution<int> distribution (weightVector.begin(), weightVector.end());

  	int sampledId = distribution(gen);

	return sampledId;
}

int sampleFromUniformDist(int maxSize)
{
  	// default_random_engine generator;
	random_device rd;
    mt19937 gen(rd());
	uniform_int_distribution<int> distribution (0, maxSize);

  	int sampledId = distribution(gen);

	return sampledId;
}

vector<VertexIdx> unorderedMapToVec(unordered_map<VertexIdx, int> localMap)
{
	vector<VertexIdx> combinedNeighborsVector;
	unordered_map<VertexIdx, int> :: iterator mapIt;

	for(mapIt = localMap.begin(); mapIt != localMap.end(); ++mapIt)
	{
		combinedNeighborsVector.push_back(mapIt->first);
	}
	return combinedNeighborsVector;
}

vector<i64> getSizesOfPivotNeighborSetsForAllSubGraphlets(vector<vector<VertexIdx>> allMinNeighborSetsRj)
{
	vector<i64> allPivotSetSizesRj;	
	for(unsigned int i = 0; i < allMinNeighborSetsRj.size(); i++)
	{
		allPivotSetSizesRj.push_back(allMinNeighborSetsRj[i].size());
	}

	return allPivotSetSizesRj; 
}

vector<vector<int>> getCombinations(int N, int q)
{
	// cout << "Getting Combinations ....\n";
    std::string bitmask(q, 1); // q leading 1's
    bitmask.resize(N, 0); // N-q trailing 0's
    vector<vector<int>> nCkSets;

    // print integers and permute bitmask
    do {
        vector<int> aset;
        for (int i = 0; i < N; ++i) // [0..N-1] integers
        {
            if (bitmask[i])
            {
				// std::cout << " " << i;
              	aset.push_back(i);
            }
        }
        // std::cout << std::endl;
        nCkSets.push_back(aset);

    } while (std::prev_permutation(bitmask.begin(), bitmask.end()));

    return nCkSets;
}


/*
vector<VertexIdx> getIntersectionOfTwoVectors(vector<VertexIdx> set1, vector<VertexIdx> set2)
{
	vector<VertexIdx> intersectedSet;
	unordered_map<VertexIdx, int> intersectionMap;

	if(set1.size() > set2.size())
	{
		for(unsigned int i = 0; i < set1.size(); i++)
		{
			intersectionMap[set1[i]] = 1;
		}

		for(unsigned int i = 0; i < set2.size(); i++)
		{
			if(intersectionMap.find(set2[i]) != intersectionMap.end())
			{
				intersectedSet.push_back(set2[i]);
			}
		}
	}
	else
	{
		for(unsigned int i = 0; i < set2.size(); i++)
		{
			intersectionMap[set2[i]] = 1;
		}

		for(unsigned int i = 0; i < set1.size(); i++)
		{
			if(intersectionMap.find(set1[i]) != intersectionMap.end())
			{
				intersectedSet.push_back(set1[i]);
			}
		}	
	}

	return intersectedSet;
}
*/


/*
unordered_map<VertexIdx, int> getUnionOfNeighbors(vector<VertexIdx> set1, vector<VertexIdx> set2)
{
	unordered_map<VertexIdx, int> combinedSet;
	for(unsigned int  i = 0; i < set1.size(); i++)
	{
		combinedSet[set1[i]] = 1;
	}
	for(unsigned int  j = 0; j < set1.size(); j++)
	{
		combinedSet[set2[j]] = 1;
	}

	return combinedSet;
}
*/
