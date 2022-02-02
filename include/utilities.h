#ifndef UTILITIES_H 
#define UTILITIES_H

using namespace std;


uint64_t next();
VertexIdx getRandomStartPoint(Count numVertices);
// Read file and get number of vertices and all the edges as a vector of vectors...
pair<Count, vector<ePair>> getSetOfEdges(string graphFileName);
vector<VertexIdx> getUnionOfNeighbors(vector<VertexIdx> set1, vector<VertexIdx> set2);
vector<VertexIdx> getIntersectionOfTwoVectors(vector<VertexIdx> set1, vector<VertexIdx> set2);
int sampleFromDiscreteDistOfR_prev(vector<i64> weightVector);
int sampleFromUniformDist(int maxSize);
vector<VertexIdx> unorderedMapToVec(unordered_map<VertexIdx, int> localMap);
vector<i64> getSizesOfPivotNeighborSetsForAllSubGraphlets(vector<vector<VertexIdx>> allMinNeighborSetsRj);
vector<vector<int>> getCombinations(int N, int q);
int notInList(vector<VertexIdx> initialRandomWalkNodes, VertexIdx nextNode);

#endif 	// UTILITIES_H

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