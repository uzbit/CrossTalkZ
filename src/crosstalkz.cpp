/*
CrossTalkZ - Statistical tool to assess crosstalk enrichment between node groupings in a network.
Copyright (C) 2013  Ted McCormack

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contents:
This file contains function definitions for loading gene groups, 
randomizing Graphs, calculating statistics, and functions for the
various tests that are presented in the paper: 
"Statistical assessment of gene group crosstalk enrichment in networks", 
2013, PLOS ONE.

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/graph/iteration_macros.hpp>

#include "crosstalkz.h"
#include "defines.h"

using namespace std;
using namespace boost;

int numSimIter = 100;
float cutoffScore = 0.0;
bool useCutoff = false;
int modeFlag = MODE_0;
bool allVsall = true;
int methodFlag = METHOD_DEFAULT;
int minimumGenesForGroup = 10;
bool doClusteringCoeff = false;
bool doHyper = false;

//	the 0th element is the  from origNetwork with gene string
//	the 1st element is the  from randNetwork with gene string
map<string, vector<Graph::Node> > geneVertMap; 	
float minObsLinks = 3;		//not used
float minExpLinks = 0.3;	//not used

#define DEGREE_BIN(x)		((int)round(log(x)+1))
#define MIN(a,b) ((a<b)?a:b)
#define MAX(a,b) ((a>b)?a:b)

#if VERBOSE
vector<float> rVals;
vector<float> smetricRatio;
#endif

/****************UTILITIY FUNCTIONS****************/

long double nCk(int n, int k){
	long double ret = 1;
	if (k >= 0 && n >= k)
		for(int i = 1; i <= k; i++)
			ret *= (n-(k-i))/(i+0.0);
	return ret;
}

long double pHyper(int n, int m, int k, int N){
	return nCk(m, k)*nCk(N-m, n-k)/nCk(N,n);
}

void zeroLinks(Graph &network)
{
	for (Graph::link_range_t er = network.getLinks(); er.first != er.second; er.first++)
		network.properties(*er.first).weight = 0;
}

template<class T1, class T2 >
bool keyInMap(const T1 &k, const map<T1, T2 > &m)
{
	return m.find(k) != m.end();
}

long double calculatePvalueFromZscore(long double x)
{
	return erfcl(fabsl(x)/M_SQRT2);
}

template <class T>
float calculateReducedChiSquare(vector<T> &dataSet, const float mean, const float stdDev);

//float calculateClusteringCoeffForGroup(const Graph &graph, const vector<GeneGroup> &groups, int groupIndex, map<string, vector<Graph::Node> > &geneVertMap, int geneVertIndex);

bool groupSort(const GeneGroup& g1, const GeneGroup& g2)
{
  return (g1.groupId < g2.groupId);
}

bool pValueSort(const pair<string, float> &p1, const pair<string, float> &p2)
{
  return (p1.second < p2.second);
}

void writeLog(const Graph &origNet, const Graph &randNet){
	long long s1 = calculateSmetricNetwork(randNet);
	long long s2 = calculateSmetricNetwork(origNet);
	cout << "Random network s-metric = " << s1 << endl;
	cout << "Original network s-metric = " << s2 << endl;
	cout << "Ratio random/original s-metric = " << (double)(s1)/(double)(s2) << endl;
	float r1 = calculateRfromNetwork(randNet);
	float r2 = calculateRfromNetwork(origNet);
	rVals.push_back(r1);
	smetricRatio.push_back((double)s1/(double)s2);
	cout << "\nRandom network r-metric = " << r1 << endl;
	cout << "Original network r-metric = " << r2 << endl;
	cout << "Ratio random/original r-metric = " << (double)(r1)/(double)(r2) << endl;
}

//validateConnectivities:
//	Checks the node connectivities of the randNet network against the connectivities of the origNet network.
//	Generates a list of errors sorted by size of differences in increasing order.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph to compare with origNet
//	errors: a vector that afterwards contains pairs of vertices of origNet and randNet that 
//			have different connectivities between origNet and randNet sorted from smallest difference to largest
//	returns true if randNet has the same connectivities as origNet else false
bool validateConnectivities(const Graph &origNet, 
							const Graph &randNet, 
							vector<pair<Graph::Node, Graph::Node> > &errors);

//fixConnectivityErrors:
//	Fixes the node connectivities of the randNet network to match connectivities of the origNet network
//	using the errors found by validateConnectivities
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a randomized Graph to fix
//	errors: a vector containing pairs of vertices of origNet and randNet that 
//			have different connectivities between origNet and randNet sorted from smallest difference to largest
void fixConnectivityErrors(const Graph &origNet, 
						   Graph &randNet, 
						   vector<pair<Graph::Node, Graph::Node> > &errors);

 
/****************FUNCTION DEFINITIONS****************/


void readGeneGroups(const Graph &origNet, vector<GeneGroup> &groups, map<string, vector<string> > &geneGroupMap, string path, stringstream &ss)
{
	ifstream file;
	string line, currentGroupID, currentGeneID;
	vector<string> lineVals; 
	currentGroupID = "";
	bool geneInNetwork;
	map<string, int> groupIndexMap;
	map<string, bool> countNotInNetwork, totalGenesInput;
	
	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
	
	cout << endl <<"Reading groups from "<< path << " ..." << endl;

	groups.clear();
	geneGroupMap.clear();
	groupIndexMap.clear();
	
	while(getline(file, line))
	{
		GeneGroup thisGroup;
		thisGroup.inputFilePath = path;
		
		split(lineVals, line, is_any_of(", \t"));
		for (int c = 0; c < (int)lineVals.size(); c++) //clean out peoples trash.
			if (lineVals[c] == ""){
				lineVals.erase(lineVals.begin()+c);
				c = 0;
			}
		if (lineVals.size() < 2)
			continue;
		
		to_upper(lineVals[GROUP_GENE]);
		to_upper(lineVals[GROUP_ID]);
		if (lineVals.size() > GROUP_SPE)
			to_lower(lineVals[GROUP_SPE]);
		if (lineVals.size() > GROUP_SYS)
			to_upper(lineVals[GROUP_SYS]);
			
		currentGeneID = lineVals[GROUP_GENE];
		currentGeneID.erase(currentGeneID.find_last_not_of(" \n\r\t")+1);
		totalGenesInput[currentGeneID] = true;
		
		geneInNetwork = (getNodeById(origNet, currentGeneID) != NULL);
		
		currentGroupID = lineVals[GROUP_ID];
		currentGroupID.erase(currentGroupID.find_last_not_of(" \n\r\t")+1);
		
		if (lineVals.size() > GROUP_SPE){
			thisGroup.groupSpe = lineVals[GROUP_SPE];
			thisGroup.groupSpe.erase(thisGroup.groupSpe.find_last_not_of(" \n\r\t")+1);
		}
		if (lineVals.size() > GROUP_SYS){
			thisGroup.groupSys = lineVals[GROUP_SYS];
			thisGroup.groupSys.erase(thisGroup.groupSys.find_last_not_of(" \n\r\t")+1);
		}
		if (lineVals.size() > GROUP_DESC){
			thisGroup.groupDesc = lineVals[GROUP_DESC];
			thisGroup.groupDesc.erase(thisGroup.groupDesc.find_last_not_of("\n\r\t")+1);
		}
		
		thisGroup.groupGenes.clear();
		
		if (!keyInMap(currentGroupID, groupIndexMap))
		{
			thisGroup.groupId = currentGroupID;
			groups.push_back(thisGroup);
			groupIndexMap[currentGroupID] = groups.size()-1;
		}
		
		if(geneInNetwork)
		{
			groups[groupIndexMap[currentGroupID]].groupGenes.push_back(currentGeneID);
			geneGroupMap[currentGeneID].push_back(currentGroupID);
		}
		else
			countNotInNetwork[currentGeneID] = true;
	}
	file.close();
	
	int totalGroups =  (int)groups.size();
	
	//remove groups from the list that do not meet some requirements
	for (int i = 0; i < (int)groups.size();)
		if ((int)groups[i].groupGenes.size() < minimumGenesForGroup) //this is the requirement to remove
		{
			//also need to remove all occurances from the geneGroup Map
			for (map<string, vector<string> >::iterator it = geneGroupMap.begin(); it != geneGroupMap.end(); it++)
				for (int j = 0; j < (int)it->second.size();)
					if (it->second[j] == groups[i].groupId)
					{
						it->second.erase(it->second.begin()+j);
						j = 0;
					}
					else
						j++;

			groups.erase(groups.begin()+i);
			i = 0;			
		}
		else
			i++;
		
	if (!groups.size())
	{
		cout << "There were no valid groups loaded. Verify group file format and existence in network." << endl;
		exit(1);
	}

	sort(groups.begin(), groups.end(), groupSort);


	map<string, bool> testUnique;
	for (int i = 0; i < (int)groups.size(); i++)
	for (int j = 0; j < (int)groups[i].groupGenes.size(); j++)
		if(!keyInMap(groups[i].groupGenes[j], testUnique))
			testUnique[groups[i].groupGenes[j]] = true;

	ss << "Total number of groups input: " << totalGroups << endl;
	ss << "Total number of unique genes in the set of groups: " << totalGenesInput.size() << endl;
	ss << "Number of groups with at least " << minimumGenesForGroup << " gene members (final number of groups): " << groups.size() << endl;
	ss << "Number of unique group genes not found in the network: " << countNotInNetwork.size() << endl;
	ss << "Number of unique genes in the set of groups and in the network: "<< testUnique.size() << endl;

#if VERBOSE	
/*	For degree distribution:
	map <int, int> degHist;
	for (Graph::node_range_t vr = origNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		degHist[origNet.getNodeDegree(*vr.first)] += 1;
	}
	
	for (map<int,int>::iterator it=degHist.begin() ; it != degHist.end(); it++)
    	cout << (*it).first << " => " << (*it).second << endl;
*/
	
	/*for (map<string, vector<string>::iterator i = geneToGroupMap.begin() ; i != geneToGroupMap.end(); i++)
		for (int j = 0; j < i.second.size(); j++)
			cout << i.first << "\t" << i.second[j];*/
	/*for (int i = 0; i < (int)groups.size(); i++)
	{	
		cout << groups[i].groupId << endl;
		for (int j = 0; j < (int)groups[i].groupGenes.size(); j++)
			cout <<"\t"<< groups[i].groupGenes[j] << endl;
	}
	cout << groups[groups.size()-1].groupGenes.size() << "<< size\n";*/
#endif

}


int generateRandomNetworkLinkSwap(const Graph &origNet, Graph &randNet)
{
	vector<Graph::Link> links;
	LinkProperties link;
	link.weight = 1.0;
	
	for (Graph::link_range_t er = randNet.getLinks(); er.first != er.second; er.first++)
		links.push_back(*er.first);
	
	int randIndex1, randIndex2;
	bool test;
	Graph::Node v1, v2, v3, v4;
	map<string, bool> hasTested; 
	stringstream ss;
	int countSwaps = 0;
	
	while(hasTested.size()/2 <= links.size() && links.size() >= 2)
	{
		do
		{
			randIndex1 = rand()%links.size();
			randIndex2 = rand()%links.size();
			
			randNet.getNodesByLink(links[randIndex1], v1, v2);
			randNet.getNodesByLink(links[randIndex2], v3, v4);
	
			test = (randIndex1 != randIndex2)
				&& ((!randNet.hasLink(v1, v3) && !randNet.hasLink(v2, v4))
				|| (!randNet.hasLink(v1, v4) && !randNet.hasLink(v2, v3)))
				&& v1 != v3 && v2 != v3 && v4 != v1 && v4 != v2;
			
			if (!test && (randIndex1 != randIndex2))
			{
				ss << randIndex1 << " " << randIndex2;
				hasTested[ss.str()] = true;
				ss << randIndex2 << " " << randIndex1;
				hasTested[ss.str()] = true;
			}
			
			//cout << hasTested.size()/2 << " " << links.size() <<endl;
		} while (!test && hasTested.size()/2 <= links.size() && links.size()>=2);
		
		if (test)
		{
			if (!randNet.hasLink(v1, v3) && !randNet.hasLink(v2, v4))
			{
				randNet.RemoveLink(v1, v2);
				randNet.RemoveLink(v3, v4);
				randNet.AddLink(v1, v3, link);
				randNet.AddLink(v2, v4, link);	
				countSwaps++;
			}
			else if (!randNet.hasLink(v1, v4) && !randNet.hasLink(v2, v3))
			{
				randNet.RemoveLink(v1, v2);
				randNet.RemoveLink(v3, v4);
				randNet.AddLink(v1, v4, link);
				randNet.AddLink(v2, v3, link);
				countSwaps++;
			}
			if (randIndex1 > randIndex2)
			{
				links.erase(links.begin()+randIndex1);
				links.erase(links.begin()+randIndex2);
				hasTested.clear();
			}
			else
			{
				links.erase(links.begin()+randIndex2);
				links.erase(links.begin()+randIndex1);
				hasTested.clear();
			}
		}
	}
	
#if VERBOSE
	vector<pair<Graph::Node, Graph::Node> > errors;
	if(validateConnectivities(origNet, randNet, errors))
		printf("Conserved connectivity.\n");
	else
		printf("Failed to conserve connectivity.\n");
	
	writeLog(origNet, randNet);
#endif
	
	return 2*countSwaps;
}

bool generateRandomNetworkLabelSwap(const Graph &origNet, Graph &randNet, map<int, vector<Record> > &degRecordsMap)
{
	vector<Record> randRecords;
	
	//create a list of vertices paired with the connectivity 
	//as in the original network (refered to as records)
	for (Graph::node_range_t vr = randNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		Record record;
		record.node = *vr.first;
		record.degree = randNet.getNodeDegree(record.node); //at this point randNet is exact copy of origNet
		randRecords.push_back(record);
	}
	
	//randomize the list of records
	random_shuffle(randRecords.begin(), randRecords.end());
	
	int randIndex;
	string swaplabel;
	vector<Record> *vec1;
	for (int k = 0; k < (int)randRecords.size(); k++)
	{
		vec1 = &(degRecordsMap[DEGREE_BIN(randRecords[k].degree)]);
		randIndex = rand()%(vec1->size());
		
		swaplabel = randNet.properties(randRecords[k].node).geneId;
		
		randNet.properties(randRecords[k].node).geneId = randNet.properties((*vec1)[randIndex].node).geneId;
		randNet.properties((*vec1)[randIndex].node).geneId = swaplabel;
	}

#if VERBOSE
	writeLog(origNet, randNet);
#endif
	
	return true;	
}

bool generateRandomNetworkSecondOrder(const Graph &origNet, Graph &randNet, map<int, vector<Record> > &degRecordsMap)
{
	int randNum, randIndex;
	vector<Record> randRecordsAll, randRecordsAvail;
	vector<int> randRecordIndices;
	LinkProperties link;
	bool test = true;
	link.weight = 1.0;
	
	cout << "Generating random network... "<<endl;
	/* NOTE: Assuming that at this point randNet = origNet!!!!! */
		
	//create a list of vertices paired with the connectivity 
	//as in the original network (refered to as records)
	for (Graph::node_range_t vr = randNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		Record record;
		record.node = *vr.first;
		record.degree = randNet.getNodeDegree(record.node); //at this point randNet is exact copy of origNet
		//record.smetric = calculateSmetricNode(randNet, record.);
		randRecordsAll.push_back(record);
	}
   
	//randomize the list of records
	random_shuffle(randRecordsAll.begin(), randRecordsAll.end());
	
	//clear the random network of all links
	randNet.RemoveAllLinks();
	
	for (int k = 0; k < (int)randRecordsAll.size(); k++)
	{			
		for (int j = 0; j < (int)randNet.properties(randRecordsAll[k].node).connectedDegrees.size(); j++)
		{
			randRecordsAvail = degRecordsMap[randNet.properties(randRecordsAll[k].node).connectedDegrees[j]];
			
			//generate a list of indices of randRecordsAvail indices
			randRecordIndices.clear();
			for (int i = 0; i < (int)randRecordsAvail.size(); i++)
				randRecordIndices.push_back(i);
	
			//pick a random record index out of the list of indices
			randIndex = rand()%randRecordIndices.size();
			randNum = randRecordIndices[randIndex];
			
			//test if this  is ok to connect to (!test is ok)		
			test = (randRecordsAll[k].node == randRecordsAvail[randNum].node
					|| randNet.hasLink(randRecordsAll[k].node, randRecordsAvail[randNum].node) 
					|| randNet.getNodeDegree(randRecordsAll[k].node) == randRecordsAll[k].degree
					|| randNet.getNodeDegree(randRecordsAvail[randNum].node) == randRecordsAvail[randNum].degree);
			
			//if not find a different one and delete that index out of the list of record indices
			while(test)
			{
				randRecordIndices.erase(randRecordIndices.begin()+randIndex);
				
				if (!randRecordIndices.size())
					break;
						
				randIndex = rand()%randRecordIndices.size();
				randNum = randRecordIndices[randIndex];
					
				test = (randRecordsAll[k].node== randRecordsAvail[randNum].node
						|| randNet.hasLink(randRecordsAll[k].node, randRecordsAvail[randNum].node) 
						|| randNet.getNodeDegree(randRecordsAll[k].node) == randRecordsAll[k].degree
						|| randNet.getNodeDegree(randRecordsAvail[randNum].node) == randRecordsAvail[randNum].degree);
			}
					
			if (!test) 
			{
				bool doBreak = false;
				
				//add the link to the network
				randNet.AddLink(randRecordsAll[k].node, randRecordsAvail[randNum].node, link);	
													
				//instead of checking all the records, only check two (the extra case is if k < randNum and k is removed)
				if (randRecordsAll[k].degree == randNet.getNodeDegree(randRecordsAll[k].node))
				{
					randRecordsAll.erase(randRecordsAll.begin()+k);
					doBreak = true;
				}
				if (randRecordsAvail[randNum].degree == randNet.getNodeDegree(randRecordsAvail[randNum].node))
				{
					randRecordsAvail.erase(randRecordsAvail.begin()+randNum);
					doBreak = true;
				}
				if (doBreak)
				{
					if (randRecordIndices.size())
						k --;
					break;
				}
			}
			
			if (!randRecordsAll.size() || !randRecordsAvail.size() || !randRecordIndices.size())
				break;
					
		}//end for j
		//cout << "k: " << k << "\tsize: " << randRecordsAll.size() << endl;
	}//end for k	

	//validate and fix the connectivity errors
	vector<pair<Graph::Node, Graph::Node> > errors;
	if (!validateConnectivities(origNet, randNet, errors))
		fixConnectivityErrors(origNet, randNet, errors);

	if (validateConnectivities(origNet, randNet, errors))
	{

#if VERBOSE
		writeLog(origNet, randNet);
#endif
		return true;
	}
	
#if DEBUG
	for (int i = 0; i < (int)errors.size(); i++)
		printf("error %d has %d in orig and %d in rand\n", i, origNet.getNodeDegree(errors[i].first), randNet.getNodeDegree(errors[i].second));
#endif

	int sum = 0;
	for (int i = 0; i < (int)errors.size(); i++)
		sum += (int)abs(origNet.getNodeDegree(errors[i].first) - randNet.getNodeDegree(errors[i].second));	
	printf("***Warning*** Randomization failed to conserve connectivities.\n");
	printf("***Warning*** There was a difference of %d links between the original and randomized network\n", sum);

	return false;
}

bool generateRandomNetworkAssignment(const Graph &origNet, Graph &randNet)
{
	int randNum, randIndex, numToGo;
	vector<pair<Graph::Node, int> > randRecords;
	vector<int> randRecordIndices;
	LinkProperties link;
	bool test = true, doBreak = false;
	
	cout << "Generating random network... "<<endl;

	link.weight = 1.0;
			
	//create a list of vertices paired with the connectivity 
	//as in the original network (refered to as records)
	for (Graph::node_range_t vr = randNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		pair<Graph::Node, int> record;
		record.first = *vr.first;
		record.second = randNet.getNodeDegree(record.first); //at this point randNet is exact copy of origNet
		randRecords.push_back(record);
	}
   
	//randomize the list of records
	random_shuffle(randRecords.begin(), randRecords.end());
				
	//clear the random network of all links
	randNet.RemoveAllLinks();
	
	for (int k = 0; k < (int)randRecords.size();) //k++ below
	{		
		numToGo = randRecords[k].second - randNet.getNodeDegree(randRecords[k].first);
				
		//generate a list of indices of randRecords
		randRecordIndices.clear();
		for (int i = 0; i < (int)randRecords.size(); i++)
			randRecordIndices.push_back(i);
			
		for (int j = 0; j < numToGo; j++)
		{
			//pick a random record index out of the list of indices
			randIndex = rand()%randRecordIndices.size();
			randNum = randRecordIndices[randIndex];
			doBreak = false;
				
			//test if this  is ok to connect to (!test is ok)		
			test = (k == randNum
					|| randNet.hasLink(randRecords[k].first, randRecords[randNum].first) 
					|| randNet.getNodeDegree(randRecords[k].first) == randRecords[k].second
					|| randNet.getNodeDegree(randRecords[randNum].first) == randRecords[randNum].second);
				
			//if not find a different one and delete that index out of the list of record indices
			while(test)
			{
				randRecordIndices.erase(randRecordIndices.begin()+randIndex);
				
				if (!randRecordIndices.size())
					break;
						
				randIndex = rand()%randRecordIndices.size();
				randNum = randRecordIndices[randIndex];
					
				test = (k == randNum
						|| randNet.hasLink(randRecords[k].first, randRecords[randNum].first) 
						|| randNet.getNodeDegree(randRecords[k].first) == randRecords[k].second
						|| randNet.getNodeDegree(randRecords[randNum].first) == randRecords[randNum].second);
			}
					
			if (!test) 
			{
				//add the link to the network
				randNet.AddLink(randRecords[k].first, randRecords[randNum].first, link);	
													
				//instead of checking all the records, only check two (the extra case is if k < randNum and k is removed)
				if (randRecords[k].second == randNet.getNodeDegree(randRecords[k].first))
				{
					randRecords.erase(randRecords.begin()+k);
					doBreak = true;
				}
				if (randRecords[randNum].second == randNet.getNodeDegree(randRecords[randNum].first))
				{
					randRecords.erase(randRecords.begin()+randNum);
					doBreak = true;
				}
				else if (randNum > 0)
					if (randRecords[randNum-1].second == randNet.getNodeDegree(randRecords[randNum-1].first))
					{
						randRecords.erase(randRecords.begin()+(randNum-1));
						doBreak = true;
					}
						
				if (doBreak)
				{
					k = 0;
					break;
				}
			}
						
			if (!randRecords.size() || !randRecordIndices.size())
				break;
					
		}//end for j
		
		if (!doBreak)
			k++;
		
	}//end for k	

	//validate and fix the connectivity errors
	vector<pair<Graph::Node, Graph::Node> > errors;
	if (!validateConnectivities(origNet, randNet, errors))
		fixConnectivityErrors(origNet, randNet, errors);

	if (validateConnectivities(origNet, randNet, errors))
	{
#if VERBOSE
		writeLog(origNet, randNet);
#endif
		return true;
	}
	
#if DEBUG
	for (int i = 0; i < (int)errors.size(); i++)
		printf("error %d has %d in orig and %d in rand\n", i, 
			origNet.getNodeDegree(errors[i].first), randNet.getNodeDegree(errors[i].second));
#endif

	int sum = 0;
	for (int i = 0; i < (int)errors.size(); i++)
		sum += (int)abs(origNet.getNodeDegree(errors[i].first) - randNet.getNodeDegree(errors[i].second));	
	printf("***Warning*** Randomization failed to conserve connectivities.\n");
	printf("***Warning*** There was a difference of %d links between the original and randomized network\n", sum);
	
	return false;
}


void fixConnectivityErrors(const Graph &origNet, Graph &randNet, 
				vector<pair<Graph::Node, Graph::Node> > &errors)
{
	int numToGo;
	Graph::Node v1=NULL, v2=NULL;
	Graph::link_range_t er;
	Graph::link_iter lastLink;
	LinkProperties link;
	bool test = false;
			
	link.weight = 1.0;


	//This algorithm has two parts. 
	//First, look for nodes that have an odd number of connectivity errors. (There should be an even number of these)
	//In order to fix these errors, you have to break a link to a node with only one
	//connected link. Then you can form an link between the node in errors and the one with odd degree
	//and use the other node to connect to the next node in errors that has odd connectivity error.	 
	//Second, after the odd errors are fixed, it fixes the even connectivity errors by finding sufficient links
	//to break and adding both links to the node with even errors until it has the correct connectivity.

#if VERBOSE	
	cout << "error size: " << errors.size() << endl;
	clock_t start = clock();
#endif	
	//this case fixes all of the odd numbered links first
	for (int i = 0; i < (int)errors.size(); i++)
	{
		numToGo = origNet.getNodeDegree(errors[i].first) - randNet.getNodeDegree(errors[i].second);
		if (numToGo % 2 == 1)
		{
			int nextOdd = -1;
			for (int k = i+1; k < (int)errors.size(); k++)
				if ((origNet.getNodeDegree(errors[k].first) - randNet.getNodeDegree(errors[k].second)) % 2 == 1)
					nextOdd = k;
			
			if (nextOdd == -1) //should not happen, just for safety. 
				break;
			  
			er = randNet.getLinks();
			for (; er.first != er.second; )
			{
				randNet.getNodesByLink(*er.first, v1, v2);
			
				test = (randNet.getNodeDegree(v1)%2 == 1) || (randNet.getNodeDegree(v2)%2 == 1);
				test = test && ((v1 != errors[i].second && v2 != errors[nextOdd].second) && (v1 != errors[nextOdd].second && v2 != errors[i].second));
				test = test && (((!randNet.hasLink(v1, errors[i].second) && !randNet.hasLink(v2, errors[nextOdd].second)) || (!randNet.hasLink(v1, errors[nextOdd].second) && !randNet.hasLink(v2, errors[i].second))));
				
				er.first++;
				if (test)
					break;
			}
			//	cout << randNet.properties(v1).geneId << " " << randNet.getNodeDegree(v1) << endl;
			//	cout << randNet.properties(v2).geneId << " " << randNet.getNodeDegree(v2) << endl;

			// in this case there will always be an even number of errors with 
			// odd numbers of connections
			if (test)
			{
				bool didAdd = false;
				if (!randNet.hasLink(v1, errors[i].second) && !randNet.hasLink(v2, errors[nextOdd].second))
				{
					randNet.AddLink(v1, errors[i].second, link, link);
					randNet.AddLink(v2, errors[nextOdd].second, link, link);
					didAdd = true;
				}
				else if (!randNet.hasLink(v2, errors[i].second) && !randNet.hasLink(v1, errors[nextOdd].second))
				{
					randNet.AddLink(v2, errors[i].second, link, link);
					randNet.AddLink(v1, errors[nextOdd].second, link, link);
					didAdd = true;
				}
				if (didAdd)
					randNet.RemoveLink(v1, v2);
			}//end if test
		}//end if odd
	}
	
#if VERBOSE
	printf("odd loops took: %f sec\n", (clock()-start)/(CLOCKS_PER_SEC+0.0));
	start = clock();
#endif

	
	//this case fixes all of the even ones.
	for (int i = 0; i < (int)errors.size(); i++)
	{
		numToGo = origNet.getNodeDegree(errors[i].first) - randNet.getNodeDegree(errors[i].second);
		//at this point numToGo should be even
		for (int j = 0; j < (numToGo/2); j++)
		{
			if (j == 0) er = randNet.getLinks();
			for (; er.first != er.second;)
			{
				randNet.getNodesByLink((*er.first), v1, v2);
				test = (v1 != errors[i].second 
					 && v2 != errors[i].second
					 && !randNet.hasLink(v1, errors[i].second)
					 && !randNet.hasLink(v2, errors[i].second));
				
				er.first++;				
				if (test)
					break;
			}
			
			if (test)
			{
				randNet.RemoveLink(v1, v2);
				randNet.AddLink(v1, errors[i].second, link, link);
				randNet.AddLink(v2, errors[i].second, link, link);
			}
		}
	}
#if VERBOSE
	printf("even loops took: %f sec\n", (clock()-start)/(CLOCKS_PER_SEC+0.0));
#endif 

}

bool validateConnectivities(const Graph &origNet, const Graph &randNet, 
				vector<pair<Graph::Node, Graph::Node> > &errors)
{
	bool valid = true;
	Graph::Node v1, v2;

	errors.clear();
	for (Graph::node_range_t vp1 = origNet.getNodes(), vp2 = randNet.getNodes(); vp1.first != vp1.second && vp2.first != vp2.second; )
	{
		v1 = *vp1.first;
		v2 = *vp2.first;
		if (origNet.properties(v1).geneId != randNet.properties(v2).geneId)
		{
			valid = false;
			printf("Error: node mismatch %s in orig is %s in rand\n", origNet.properties(v1).geneId.c_str(), randNet.properties(v2).geneId.c_str());
			printf("This should never happen...\n");
			exit(1);
		}
		if (origNet.getNodeDegree(v1) != randNet.getNodeDegree(v2))
		{
			int i = 0;
			pair<Graph::Node, Graph::Node> thisPair;
			thisPair.first = v1;
			thisPair.second = v2;
			//sort the errors by increasing difference in connectivity from original network
			if (errors.size())
				for (i = 0; i < (int)errors.size(); i++)
					if (origNet.getNodeDegree(errors[i].first) - randNet.getNodeDegree(errors[i].second) >
						origNet.getNodeDegree(v1) - randNet.getNodeDegree(v2))
							break;
						
			errors.insert(errors.begin()+i, thisPair);
			valid = false;
		}
		vp1.first++; 
		vp2.first++;
	}
#if VERBOSE	
	for (int i = 0; i < (int)errors.size(); i++)
	{	
		printf("Error %d has %d in orig and %d in rand, delta = %d\n", i, 
		origNet.getNodeDegree(errors[i].first), randNet.getNodeDegree(errors[i].second), 
		origNet.getNodeDegree(errors[i].first)-randNet.getNodeDegree(errors[i].second));
	}
#endif		

	return (valid && (origNet.getNodeCount() == randNet.getNodeCount()) && (origNet.getLinkCount() == randNet.getLinkCount()));
}

bool getTest(vector<string> *ggmp1, const string &g1, const int &p1s, vector<string> *ggmp2, const string &g2, const int &p2s)
{
	bool test = false;
	int k;
	switch (modeFlag)
	{
		//this code could be simplified, but it might be confusing
		default:
		case MODE_0:
			//if either gene is in both groups, dont count
			for (k = 0; k < p1s; k++)
				if ((*ggmp1)[k] == g2)
				{
					test = true;
					break;
				}
			if (!test)
			for (k = 0; k < p2s; k++)
				if ((*ggmp2)[k] == g1)
				{
					test = true;
					break;
				}	
			break;
		case MODE_1:
			//if both genes are in both groups, dont count
			for (k = 0; k < p1s; k++)
				if ((*ggmp1)[k] == g2)
				{
					test = true;
					break;
				}
			for (k = 0; k < p2s; k++)
				if ((*ggmp2)[k] == g1)
				{
					test = (test && true);
					break;
				}	
			break;
	}
	return test;
}


void countLinksForGroupsAll(Graph &randNet,
				vector<GeneGroup> &groups,
				map<string, Stats > &groupStats,
				map<string, vector<string> > &geneGroupMap) 
{
	Graph::Node v1, v2;
	string groupsVsStr, p1, p2, g1, g2;
	Stats *thisGroupStats;
	vector<string> *ggmp1, *ggmp2; //using these pointras for speed optimization
	int p1s, p2s;
	
	clock_t start = clock();

	cout << "Counting links between groups...";flush(cout);
	
	for (int i = 0; i < (int)groups.size(); i++)
		for (int j = 0; j <= i; j++)
		{
			groupsVsStr = groups[i].groupId + "_vs_" + groups[j].groupId;
			groupStats[groupsVsStr].linkCount.push_back(0);
			/*if (doClusteringCoeff)
			{
				groupStats[groupsVsStr].clusteringCoeff.push_back(	
				calculateClusteringCoeffForTwoGroups(randNet, groups[i], groups[j]));
			}*/
		}
	
	for (Graph::link_range_t er = randNet.getLinks(); er.first != er.second; er.first++)
	{
		randNet.getNodesByLink((*er.first), v1, v2);
		
		p1 = randNet.properties(v1).geneId;
		p2 = randNet.properties(v2).geneId;
		
		ggmp1 = &(geneGroupMap[p1]);
		ggmp2 = &(geneGroupMap[p2]);
				
		p1s = ggmp1->size();
		p2s = ggmp2->size();
		
		for (int i = 0; i < p1s; i++)
		{
			g1 = (*ggmp1)[i];
			
			for (int j = 0; j < p2s; j++)
			{
				g2 = (*ggmp2)[j];
				
				if (g1 >= g2)
					groupsVsStr = g1 + "_vs_" + g2;
				else
					groupsVsStr = g2 + "_vs_" + g1;

				thisGroupStats = &(groupStats[groupsVsStr]);
	
				//if (keyInMap(groupsVsStr, groupStats))
				{
					if (g2 == g1) //same group
						thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;
					else //different groups
					{
						if (!getTest(ggmp1, g1, p1s, ggmp2, g2, p2s))
							thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;		
					}
				}
			}
		}
	}
		
	cout << "done in " << (clock()-start)/(CLOCKS_PER_SEC+0.0) <<" seconds." << endl;
}

void countLinksForGroups12(Graph &randNet,
								vector<GeneGroup> &groups1,
								vector<GeneGroup> &groups2,
								map<string, Stats > &groupStats,
								map<string, vector<string> > &geneGroupMap1,
								map<string, vector<string> > &geneGroupMap2) 
{
	Graph::Node v1, v2;
	string groupsVsStr, p1, p2, g1, g2;
	Stats *thisGroupStats;
	vector<string> *ggmp1, *ggmp2;
	int p1s, p2s;
	
	clock_t start = clock();

	cout << "Counting links between groups...";flush(cout);
	
	for (int i = 0; i < (int)groups1.size(); i++)
		for (int j = 0; j < (int)groups2.size(); j++)
		{
			groupsVsStr = groups1[i].groupId + "_vs_" + groups2[j].groupId;
			groupStats[groupsVsStr].linkCount.push_back(0);
			/*if (doClusteringCoeff)
			{
				groupStats[groupsVsStr].clusteringCoeff.push_back(	
				calculateClusteringCoeffForTwoGroups(randNet, groups1[i], groups2[j]));
			}*/

		}
	
	for (Graph::link_range_t er = randNet.getLinks(); er.first != er.second; er.first++)
	{
		randNet.getNodesByLink((*er.first), v1, v2);
		
		p1 = randNet.properties(v1).geneId;
		p2 = randNet.properties(v2).geneId;
		
		bool p11 = keyInMap(p1, geneGroupMap1);
		bool p12 = keyInMap(p1, geneGroupMap2);
		bool p21 = keyInMap(p2, geneGroupMap1);
		bool p22 = keyInMap(p2, geneGroupMap2);
		
		if (p11 && p22)
		{
			ggmp1 = &(geneGroupMap1[p1]);
			ggmp2 = &(geneGroupMap2[p2]);
			p1s = ggmp1->size();
			p2s = ggmp2->size();
			
			for (int i = 0; i < p1s; i++)
				for (int j = 0; j < p2s; j++)
				{
					g1 = (*ggmp1)[i];
					g2 = (*ggmp2)[j];
					
					groupsVsStr = g1 + "_vs_" + g2;
					thisGroupStats = &(groupStats[groupsVsStr]);
			
					if (g2 == g1) //same group
						thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;
					else //different groups
					{
						if (!getTest(ggmp1, g1, p1s, ggmp2, g2, p2s))
							thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;
					}
				}		
		}
		
		if (p12 && p21)
		{
			ggmp1 = &(geneGroupMap1[p2]);
			ggmp2 = &(geneGroupMap2[p1]);
			p1s = ggmp1->size();
			p2s = ggmp2->size();
			
			for (int i = 0; i < p1s; i++)
				for (int j = 0; j < p2s; j++)
				{
					g1 = (*ggmp1)[i];
					g2 = (*ggmp2)[j];
					
					groupsVsStr = g1 + "_vs_" + g2;
					thisGroupStats = &(groupStats[groupsVsStr]);
			
					if (g2 == g1) //same group
						thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;
					else //different groups
					{
						if (!getTest(ggmp1, g1, p1s, ggmp2, g2, p2s))
							thisGroupStats->linkCount[thisGroupStats->linkCount.size()-1] += 1;				
					}
				}		
		}
	} //end for each link in random network
	
	cout << "done in " << (clock()-start)/(CLOCKS_PER_SEC+0.0) <<" seconds." << endl;
}

template <class T>
void calcStatFromVec(vector<T> &vec, int s, float &mean, float &std)
{
	mean = 0;
	for (int c = 0; c < s; c++)
		mean += vec[c];
	mean /= (s+0.0);	
		
	std = 0;	
	for (int c = 0; c < s; c++)
		std += powf(vec[c] - mean, 2.0);
			
	std /= (s+0.0);
	std = sqrt(std);			
}


void calculateAndWriteResultsAll(Graph &origNet,
								 vector<GeneGroup> &groups,
								 map<string, Stats > &groupStats,
								 map<string, vector<string> > &geneGroupMap,
								 string path) 
{
	ofstream file;
	float NobservedLinks, NexpectedLinks, stdDev;// CCobserved, CCexpected, , stdDevCC;
	string groupsVsStr, g1, g2;
	Stats *thisGroupStats;
	vector<pair<string, double> > sortedPValuesIntra;
	vector<pair<string, double> > sortedPValuesInter;
	vector<int> *gsm;
	map<string, int > kSuccess, nDraws, mSuccesses;
	map<string, Stats > observedGroupStats;
	int gss = 0, N = 0; 
	
	cout << "Calculating results... " << endl;
	
	countLinksForGroupsAll(origNet, groups, observedGroupStats, geneGroupMap);
	
	for (int i = 0; i < (int)groups.size(); i++)
	{
		for (int j = 0; j <= i; j++)
		{
			g1 = groups[i].groupId;
			g2 = groups[j].groupId;
			
			if (g1 >= g2)
				groupsVsStr = g1 + "_vs_" + g2;
			else
				groupsVsStr = g2 + "_vs_" + g1;
			
			thisGroupStats = &(groupStats[groupsVsStr]);
					
			gsm = &(thisGroupStats->linkCount);
			gss = (int)gsm->size();
									
			calcStatFromVec((*gsm), gss, NexpectedLinks, stdDev);
			
			NobservedLinks = observedGroupStats[groupsVsStr].linkCount[0];
				
			thisGroupStats->expectedLinks = NexpectedLinks;
			thisGroupStats->observedLinks = NobservedLinks;
			
			bool testValid = (stdDev != 0.0); //add other tests here
			
			if (testValid)
			{
				thisGroupStats->zScore = ((NobservedLinks - NexpectedLinks)/stdDev);
				thisGroupStats->pValue = calculatePvalueFromZscore(thisGroupStats->zScore);
				thisGroupStats->stdDev = stdDev;
				thisGroupStats->chiSqr = calculateReducedChiSquare((*gsm), NexpectedLinks, stdDev);

				if (i == j)
				  sortedPValuesIntra.push_back(pair<string, float>(groupsVsStr, thisGroupStats->pValue));
				else
				  sortedPValuesInter.push_back(pair<string, float>(groupsVsStr, thisGroupStats->pValue)); 
			}
			
			/*
			if (doClusteringCoeff)
			{
				calcStatFromVec(groupStats[groupsVsStr].clusteringCoeff, (int)groupStats[groupsVsStr].clusteringCoeff.size(), CCexpected, stdDevCC);
			
				CCobserved = observedGroupStats[groupsVsStr].clusteringCoeff[0];
				
				expectedCC[groupsVsStr] = CCexpected;
				observedCC[groupsVsStr] = CCobserved;
				stdsCC[groupsVsStr] = stdDevCC;
				
				testValid = (stdDevCC != 0.0); //add other tests here
			
				if (testValid)
				{
					zScoresCC[groupsVsStr] = ((CCobserved - CCexpected)/stdDevCC);
					pValuesCC[groupsVsStr] = calculatePvalueFromZscore(zScoresCC[groupsVsStr]);
					stdsCC[groupsVsStr] = stdDevCC;
					
					if (i == j)
					  sortedPValuesIntra.push_back(pair<string, float>(groupsVsStr, pValues[groupsVsStr]));
					else
					  sortedPValuesInter.push_back(pair<string, float>(groupsVsStr, pValues[groupsVsStr])); 
				}
			}*/
			
		}
	}
	
	sort(sortedPValuesIntra.begin(), sortedPValuesIntra.end(), pValueSort);
	sort(sortedPValuesInter.begin(), sortedPValuesInter.end(), pValueSort);
	
	for (int c = 1; c < (int)sortedPValuesIntra.size(); c++)
	{
		sortedPValuesIntra[c].second *= sortedPValuesIntra.size()/(sortedPValuesIntra.size() - c+0.0);
		if (sortedPValuesIntra[c].second > 1.0)
			sortedPValuesIntra[c].second = 1.0;
	}	 
	
	for (int c = 1; c < (int)sortedPValuesInter.size(); c++)
	{
		sortedPValuesInter[c].second *= sortedPValuesInter.size()/(sortedPValuesInter.size() - c+0.0);
		if (sortedPValuesInter[c].second > 1.0)
			sortedPValuesInter[c].second = 1.0;
	}
	
	if (doHyper){
		N = getTotalInputUniqueGeneCount(groups[0].inputFilePath, groups[0].inputFilePath);
		
		cout << "P-hyper using N (total unique genes in the two grops) = " << N << endl;
	
		for (int i = 0; i < (int)groups.size(); i++)
		{
			for (int j = 0; j <= i; j++)
			{
				groupsVsStr = groups[i].groupId + "_vs_" + groups[j].groupId;
				int countShared = 0;
				for (int k = 0; k < (int)groups[i].groupGenes.size(); k++)
					for (int l = 0; l < (int)groups[j].groupGenes.size(); l++)
						if (groups[i].groupGenes[k] == groups[j].groupGenes[l])
							countShared++;
				kSuccess[groupsVsStr] = countShared;
				nDraws[groupsVsStr] = MIN(groups[i].groupGenes.size(), groups[j].groupGenes.size());
				mSuccesses[groupsVsStr] = MAX(groups[i].groupGenes.size(), groups[j].groupGenes.size());
			}
		}
	}

	cout << "Writing results to "<< path << " ..." << endl;
	
	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
#if !DEBUG	
	file << "Scores for intra-group comparisons:"<<endl;
#endif	
	file << "PAIR\ttype1 type2\tintra/inter\tObserved links\tExpected Links\tZscore\tp-value\tpFDR\tstdDev\tReduced ChiSqr"<<endl;
	for (int i = 0; i < (int)groups.size(); i++)
	{
		groupsVsStr = groups[i].groupId + "_vs_" + groups[i].groupId;
		thisGroupStats = &(groupStats[groupsVsStr]);
		
		bool testValid = (thisGroupStats->stdDev != 0.0);
		long double fdr = 0;
		
		for (int c = 0; c < (int)sortedPValuesIntra.size() ; c++)
			if (groupsVsStr == sortedPValuesIntra[c].first)
			{
				fdr = sortedPValuesIntra[c].second;
				break;
			}	 

		file.precision(6);
		file	<< groupsVsStr << "\t" 
			<< groups[i].groupSys <<" "<< groups[i].groupSys << "\t" 
			<< "intra" << "\t"
			<< (int)thisGroupStats->observedLinks << "\t"
			<< thisGroupStats->expectedLinks << "\t" 
			<< ((testValid)?lexical_cast<string>((float)thisGroupStats->zScore):"NA") << "\t"
			<< ((testValid)?lexical_cast<string>(thisGroupStats->pValue):"NA") << "\t"
			<< ((testValid)?lexical_cast<string>((float)fdr):"NA") << "\t"
			<< ((testValid)?lexical_cast<string>((float)thisGroupStats->stdDev):"NA") << "\t"
			<< thisGroupStats->chiSqr << "\t"
			<< ((testValid)?(doHyper?lexical_cast<string>(pHyper(nDraws[groupsVsStr], mSuccesses[groupsVsStr], kSuccess[groupsVsStr], N)):""):"NA")<< "\t"	
			<< endl;
					
	}

#if !DEBUG	
	file << endl << "Scores for inter-group comparisons:"<<endl;
	file << "PAIR\ttype1 type2\tintra/inter\tObserved links\tExpected Links\tZscore\tp-value\tpFDR\tstdDev\tReduced ChiSqr\tp-hyper"<<endl;
#endif
		
	for (int i = 0; i < (int)groups.size(); i++)
	{
		for (int j = 0; j <= i; j++)
		{
			if (i==j)continue;
			
			groupsVsStr = groups[i].groupId + "_vs_" + groups[j].groupId;
			thisGroupStats = &(groupStats[groupsVsStr]);
		
			bool testValid = (thisGroupStats->stdDev != 0.0);
			long double fdr = 0;
			
			for (int c = 0; c < (int)sortedPValuesInter.size() ; c++)
				if (groupsVsStr == sortedPValuesInter[c].first)
				{
					fdr = sortedPValuesInter[c].second;
					break;
				}	 
			
			file.precision(6);
			file	<< groupsVsStr << "\t" 
				<< groups[i].groupSys <<" "<< groups[i].groupSys << "\t" 
				<< "inter" << "\t"
				<< (int)thisGroupStats->observedLinks << "\t"
				<< thisGroupStats->expectedLinks << "\t" 
				<< ((testValid)?lexical_cast<string>((float)thisGroupStats->zScore):"NA") << "\t"
				<< ((testValid)?lexical_cast<string>(thisGroupStats->pValue):"NA") << "\t"
				<< ((testValid)?lexical_cast<string>((float)fdr):"NA") << "\t"	
				<< ((testValid)?lexical_cast<string>((float)thisGroupStats->stdDev):"NA") << "\t"
				<< thisGroupStats->chiSqr << "\t"
				<< ((testValid)?(doHyper?lexical_cast<string>(pHyper(nDraws[groupsVsStr], mSuccesses[groupsVsStr], kSuccess[groupsVsStr], N)):""):"NA")<< "\t"
				<< endl;
		}
	}
	
	file.close();
	
/*#if DEBUG
	file.open("groupStats.out");
	for (int i = 0; i < (int)groups.size(); i++)
		for (int j = 0; j <= i; j++)
		{
			groupsVsStr = groups[i].groupId;
			groupsVsStr += "_vs_";
			groupsVsStr += groups[j].groupId;
			file << groupsVsStr<<"\t";
			for (int c = 0; c < (int)groupStats[groupsVsStr].linkCount.size(); c++)
				file << groupStats[groupsVsStr].linkCount[c] << "\t";
			file << endl;
		}
	file.close();
#endif*/
			
	
}

void calculateAndWriteResults12(Graph &origNet,
								 vector<GeneGroup> &groups1,
								 vector<GeneGroup> &groups2,
								 map<string, Stats > &groupStats,
								 map<string, vector<string> > &geneGroupMap1,
								 map<string, vector<string> > &geneGroupMap2,
								 string path) 
{
	ofstream file;
	float NobservedLinks, NexpectedLinks, stdDev;
	string groupsVsStr, g1, g2, p1, p2;
	Stats *thisGroupStats;
	vector<pair<string, long double> > sortedPValues;
 	vector<int> *gsm;
	map<string, int > kSuccess, nDraws, mSuccesses;
	map<string, Stats > observedGroupStats;
	int gss = 0, N = 0; 
						
	cout << "Calculating results... " << endl;

	countLinksForGroups12(origNet, groups1, groups2, observedGroupStats, geneGroupMap1, geneGroupMap2);
	
	for (int i = 0; i < (int)groups1.size(); i++)
	{
		for (int j = 0; j < (int)groups2.size(); j++)
		{
			g1 = groups1[i].groupId;
			g2 = groups2[j].groupId;
			
			groupsVsStr = g1 + "_vs_" + g2;
			thisGroupStats = &(groupStats[groupsVsStr]);
		
			gsm = &(thisGroupStats->linkCount);
			gss = (int)gsm->size();
			
			for (int c = 0; c < gss; c++)
				(*gsm)[c] *= ((g1==g2)?0.5:1.0); //links between same groups counted twice 
			
			calcStatFromVec((*gsm), gss, NexpectedLinks, stdDev);
			
			NobservedLinks = observedGroupStats[groupsVsStr].linkCount[0]*((g1==g2)?0.5:1.0);
						
			thisGroupStats->observedLinks = NobservedLinks; 
			thisGroupStats->expectedLinks = NexpectedLinks;
			
			bool testValid = (stdDev != 0.0); //add other tests here
			
			if (testValid)
			{
				thisGroupStats->zScore = ((NobservedLinks - NexpectedLinks)/stdDev);
				thisGroupStats->pValue = calculatePvalueFromZscore(thisGroupStats->zScore);
				thisGroupStats->stdDev = stdDev;
			    thisGroupStats->chiSqr = calculateReducedChiSquare((*gsm), NexpectedLinks, stdDev);

				int c = 0;
			    for (c = 0; c < (int)sortedPValues.size(); c++)
			  	  if (sortedPValues[c].second > thisGroupStats->pValue)
			  		  break;
				sortedPValues.insert(sortedPValues.begin()+c, pair<string, float>(groupsVsStr, thisGroupStats->pValue));
			}			
		}
	}
	
	for (int c = 1; c < (int)sortedPValues.size(); c++)
	{
		sortedPValues[c].second *= sortedPValues.size()/(sortedPValues.size() - c+0.0);
		if (sortedPValues[c].second > 1.0)
			sortedPValues[c].second = 1.0;
	}

	if (doHyper){
		N = getTotalInputUniqueGeneCount(groups1[0].inputFilePath, groups2[0].inputFilePath);
		
		cout << "P-hyper using N (total unique genes in the two grops) = " << N << endl;
	
		for (int i = 0; i < (int)groups1.size(); i++)
		{
			for (int j = 0; j < (int)groups2.size(); j++)
			{
				groupsVsStr = groups1[i].groupId + "_vs_" + groups2[j].groupId;
				int countShared = 0;
				for (int k = 0; k < (int)groups1[i].groupGenes.size(); k++)
					for (int l = 0; l < (int)groups2[j].groupGenes.size(); l++)
						if (groups1[i].groupGenes[k] == groups2[j].groupGenes[l])
							countShared++;
				kSuccess[groupsVsStr] = countShared;
				nDraws[groupsVsStr] = MIN(groups1[i].groupGenes.size(), groups2[j].groupGenes.size());
				mSuccesses[groupsVsStr] = MAX(groups1[i].groupGenes.size(), groups2[j].groupGenes.size());
			}
		}
	}

	cout << "Writing results to "<< path << " ..." << endl;
	
	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
		
	file << "PAIR\ttype1 type2\tintra/inter\tObserved links\tExpected Links\tZscore\tp-value\tpFDR\tstdDev\tReduced ChiSqr\tp-hyper"<<endl;	
	for (int i = 0; i < (int)groups1.size(); i++)
	{
		for (int j = 0; j < (int)groups2.size(); j++)
		{	
			groupsVsStr = groups1[i].groupId + "_vs_" + groups2[j].groupId;
			thisGroupStats = &(groupStats[groupsVsStr]);
		
			bool testValid = (thisGroupStats->stdDev != 0.0);

			long double fdr = 0;
			
			if (testValid)
			{
				for (int c = 0; c < (int)sortedPValues.size() ; c++)
					if (groupsVsStr == sortedPValues[c].first)
					{
						fdr = sortedPValues[c].second;
						break;
					}	 
				
				file.precision(6);
				file	<< groupsVsStr << "\t" 
					<< groups1[i].groupSys <<" "<< groups2[j].groupSys << "\t" 
					<< ((groups1[i].groupId == groups2[j].groupId)?"intra":"inter")<< "\t"					
					<< (int)thisGroupStats->observedLinks << "\t"
					<< thisGroupStats->expectedLinks << "\t" 
					<< thisGroupStats->zScore << "\t"
					<< thisGroupStats->pValue << "\t"
					<< fdr << "\t"
					<< thisGroupStats->stdDev << "\t"
					<< thisGroupStats->chiSqr << "\t"
					<< ((doHyper)?lexical_cast<string>(pHyper(nDraws[groupsVsStr], mSuccesses[groupsVsStr], kSuccess[groupsVsStr], N)):"")<< "\t"
					<< endl;
	
			}
			else
			{
				file	<< groupsVsStr << "\t" 
					<< groups1[i].groupSys <<" "<< groups2[j].groupSys << "\t" 
					<< ((groups1[i].groupId == groups2[j].groupId)?"intra":"inter")<< "\t"
					<< (int)thisGroupStats->observedLinks << "\t"
					<< thisGroupStats->expectedLinks << "\t" 
					<< "NA" << "\t"
					<< "NA" << "\t"
					<< "NA" << "\t"
					<< "NA" << "\t"
					<< "NA" << "\t"
					<< ((doHyper)?"NA":"") << "\t"
					<< endl;
			}
		}
	}
	
	file.close();
	
/*#if DEBUG
	file.open("groupStats.out");
	for (int i = 0; i < (int)groups.size(); i++)
		for (int j = 0; j <= i; j++)
		{
			groupsVsStr = groups[i].groupId;
			groupsVsStr += "_vs_";
			groupsVsStr += groups[j].groupId;
			file << groupsVsStr<<"\t";
			for (int c = 0; c < (int)groupStats[groupsVsStr].linkCount.size(); c++)
				file << groupStats[groupsVsStr].linkCount[c] << "\t";
			file << endl;
		}
	file.close();
#endif*/
			
	
}


void generateMaps(Graph &origNet,
				  Graph &randNet,
				  map<int, vector<Record> > &degRecordsMap)
{
	cout << "Generating maps...";flush(cout);
		
	Graph::Node v1;
	map<string, Graph::Node> g2vMap;
	string gene;
	
	geneVertMap.clear();
	degRecordsMap.clear();
	
	for (Graph::node_range_t vr = randNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		v1 = *vr.first;
		g2vMap[randNet.properties(v1).geneId] = v1;
		
		Record record;
		record.node = v1;
		record.degree = randNet.getNodeDegree(record.node); //at this point randNet is exact copy of origNet
		//record.smetric = calculateSmetricNode(randNet, record.);
		degRecordsMap[DEGREE_BIN(record.degree)].push_back(record);
	}
	for (Graph::node_range_t vr = origNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		v1 = *vr.first;
		gene = origNet.properties(v1).geneId;
		geneVertMap[gene].push_back(v1);
		geneVertMap[gene].push_back(g2vMap[gene]);		
	}
	
	for (Graph::node_range_t vr = randNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		v1 = *vr.first;
		for (Graph::adjacency_node_range_t avr = randNet.getAdjacentNodes(v1); avr.first != avr.second; avr.first++)
			origNet.properties(geneVertMap[randNet.properties(v1).geneId][0]).connectedDegrees.push_back(DEGREE_BIN(randNet.getNodeDegree(*avr.first)));
	}
	
	for (Graph::node_range_t vr = origNet.getNodes(); vr.first != vr.second; vr.first++)
		randNet.properties(geneVertMap[origNet.properties(*vr.first).geneId][1]).connectedDegrees = origNet.properties(*vr.first).connectedDegrees;	
	
	g2vMap.clear();
	cout << "done." << endl;
}


Graph::Node getNodeById(const Graph &g, const string &Id)
{
	if (!keyInMap(Id, geneVertMap))
		return NULL;
	return (geneVertMap[Id][g.id]);
}

void copyOrigToRand(Graph &origNet, Graph &randNet)
{
	Graph::Node v1, v2;
	LinkProperties link;
	link.weight = 1.0;

	//cout << "num links before: " << randNet.getLinkCount() << endl;
	randNet.RemoveAllLinks();
	//cout << "num links after: " << randNet.getLinkCount() << endl;
	
	for (Graph::link_range_t er = origNet.getLinks(); er.first != er.second; er.first++)
	{
		origNet.getNodesByLink(*er.first, v1, v2);
		link.weight = origNet.properties(*er.first).weight;
		randNet.AddLink(geneVertMap[origNet.properties(v1).geneId][1], geneVertMap[origNet.properties(v2).geneId][1], link); 	
	}
}


void printNetwork(const Graph &network)
{
	cout << "Nodes:" << endl;
	for (Graph::node_range_t vr = network.getNodes(); vr.first != vr.second; vr.first++)
	{
		Graph::Node v1 = *vr.first;
		cout << network.properties(v1).geneId << " " << network.getNodeDegree(v1) << endl;
		
		for (Graph::adjacency_node_range_t avr = network.getAdjacentNodes(v1); avr.first != avr.second; avr.first++)
		{	
			Graph::Node v2 = *avr.first;
			Graph::LinkPair linkPair;
			if (network.getLinkPair(v1, v2, linkPair) == LINK_BOTH)
				cout << "\t" << network.properties(v2).geneId <<"\t"<< network.properties(linkPair.first).weight << "\t"<< network.properties(linkPair.second).weight << endl;
		}
   }
}

string getMethodString(int m)
{
	switch (m)
	{
		case METHOD_LINKSWAP:
			return string("Link Permutation");
		case METHOD_ASSIGN:
			return string("Link Assignment");
		case METHOD_ASSIGN_SECOND:
			return string("Link Assignment + Second-order");
		case METHOD_LABELSWAP:
			return string("Node Label Permutation");	
	}
	return string("");
}

int getTotalInputUniqueGeneCount(string path1, string path2)
{
	string line;
	vector<string> lineVals; 
	map<string, int > uniqueGenesInPool;
	ifstream file;
	int ret = 0;
	
	file.open(path1.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path1 << endl;
		exit(1);
	}
	
	while(getline(file, line))
	{
		split(lineVals, line, is_any_of(", \t"));
		for (int c = 0; c < (int)lineVals.size(); c++) //clean out peoples trash.
			if (lineVals[c] == ""){
				lineVals.erase(lineVals.begin()+c);
				c = 0;
			}
		if (lineVals.size() < 2)
			continue;
		
		to_upper(lineVals[GROUP_GENE]);
		if (!keyInMap(lineVals[GROUP_GENE], uniqueGenesInPool)){
			uniqueGenesInPool[lineVals[GROUP_GENE]] = 1;
			ret++;
		}
	}
	file.close();
	
	file.open(path2.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path2 << endl;
		exit(1);
	}
	
	while(getline(file, line))
	{
		split(lineVals, line, is_any_of(", \t"));
		for (int c = 0; c < (int)lineVals.size(); c++) //clean out peoples trash.
			if (lineVals[c] == ""){
				lineVals.erase(lineVals.begin()+c);
				c = 0;
			}
		if (lineVals.size() < 2)
			continue;
		
		to_upper(lineVals[GROUP_GENE]);
		if (!keyInMap(lineVals[GROUP_GENE], uniqueGenesInPool)){
			uniqueGenesInPool[lineVals[GROUP_GENE]] = 1;
			ret++;
		}
	}
	file.close();

	return ret;
}

long long calculateSmetricNetwork(const Graph &g)
{
	long long s_metric = 0;
	Graph::Node v1, v2;
	for (Graph::link_range_t er = g.getLinks(); er.first != er.second; er.first++)
	{
		g.getNodesByLink(*er.first, v1, v2);
		s_metric += g.getNodeDegree(v1)*g.getNodeDegree(v2);
	}
	return s_metric;	
}

float calculateRfromNetwork(const Graph &g)
{
	long double invLinkCount = 1.0/g.getLinkCount();
	long double term1=0, term2=0, term3=0;
	long double r = 0;
	Graph::Node v1, v2;
	int d1, d2;
	for (Graph::link_range_t er = g.getLinks(); er.first != er.second; er.first++)
	{
		g.getNodesByLink(*er.first, v1, v2);
		d1 = g.getNodeDegree(v1);
		d2 = g.getNodeDegree(v2);
		term1 += d1*d2;
		term2 += 0.5*(d1+d2);
		term3 += 0.5*(d1*d1+d2*d2);
	}
	r = invLinkCount*term1 - powl(invLinkCount*term2, 2.0);
	r /= (invLinkCount*term3 - powl(invLinkCount*term2, 2.0));
	return r; 	
}

int calculateSmetricNode(const Graph &g, Graph::Node v)
{
	int sum = 0;
	for (Graph::adjacency_node_range_t avr = g.getAdjacentNodes(v); avr.first != avr.second; avr.first++)
		sum += g.getNodeDegree(*avr.first);
	return sum;	
}


#define LOWER_BIN -1.6
#define HIGHER_BIN 1.6
#define BIN_WIDTH 0.4
#define NUM_BINS (((HIGHER_BIN-LOWER_BIN)/BIN_WIDTH)+1)
template <class T>
float calculateReducedChiSquare(vector<T> &dataSet, const float mean, const float stdDev)
{
	map<int, int> bins;
	int N = dataSet.size();
	vector<float> normData;
	vector<float> bins_vals;
	vector<int> bins_counts;
	
	//normalize the dataset
	for (int i=0; i < N; i++)
	{	 
		normData.push_back((dataSet[i]-mean)/stdDev);
		//cout << "data: " << dataSet[i] << " normedData: " << normData[i] << endl;
		
	}
	//create the bins [-inf~-1e9 to LOWER_BIN, LOWER_BIN+NUM_BINS*BIN_WIDTH, HIGHER_BIN to inf~1e9) 
	bins_vals.push_back(-1e9);
	for (int i=0; i < NUM_BINS; i++)
		bins_vals.push_back(LOWER_BIN+i*BIN_WIDTH);	
	bins_vals.push_back(1e9);
	
	//count the frequency in each bin
	bins_counts.resize(bins_vals.size(), 0);
	for (int i=0; i < N; i++)
		for(int j=1; j < NUM_BINS+2; j++)
		{
			if (normData[i] < bins_vals[j])
			{
				bins_counts[j-1]++;
				break;
			}
		}
	
	//calculate chi-square
	float chiSqr = 0, exp = 0;
	for(int j=0; j < NUM_BINS+1; j++)
	{	
	//	cout << bins_vals[j] << " has: " << bins_counts[j]<<endl;
		exp = N*0.5*(erfc(-bins_vals[j+1]/M_SQRT2)-erfc(-bins_vals[j]/M_SQRT2));
		chiSqr += powf((bins_counts[j]-exp)/stdDev, 2.0);
	//	cout<<"Expected: " << exp <<endl;
	}	
	//return reduced chi-squared  = chiSqr/(N - #constraints=3 (calculated mean, stdDev and N) - 1)
	//cout << "chiSqr: " << chiSqr << ", red chiSqr: " << chiSqr/(N-3) << endl;
	return chiSqr/(N-3);
}

float calculateClusteringCoeffForTwoGroups(const Graph &graph, const GeneGroup &group1, const GeneGroup &group2)
{
	
	if (group1.groupId == group2.groupId)
		return calculateClusteringCoeffForGroupOnly(graph, group1);
	else
		return 0;
		
	{	//this next block will do all vs all instead of just intra (take out else return 0 above)
		GeneGroup combined = group1;
		bool hasGene;
		for (int i = 0; i < (int)group2.groupGenes.size(); i++)
		{
			hasGene = false;
			for (int j = 0; j < (int)group1.groupGenes.size(); j++)
				if (group2.groupGenes[i] == group1.groupGenes[j])
				{
					hasGene = true;
					break;
				}
			if (!hasGene)
				combined.groupGenes.push_back(group2.groupGenes[i]);
		}
		return calculateClusteringCoeffForGroupOnly(graph, combined);
	}
}

float calculateClusteringCoeffForGroup(const Graph &graph, const GeneGroup &group)
{
	Graph::Node v1;
	float numLinksInSubGraph, v1Deg, ret;
	
	ret = 0;
	for (int i = 0; i < (int)group.groupGenes.size(); i++)
	{
		v1 = getNodeById(graph, group.groupGenes[i]);
		v1Deg = graph.getNodeDegree(v1);
		numLinksInSubGraph = 0;
		
		if (v1Deg <= 1)
			continue;
		
		for (Graph::adjacency_node_range_t avr1 = graph.getAdjacentNodes(v1); avr1.first != avr1.second; avr1.first++)
			for (Graph::adjacency_node_range_t avr2 = graph.getAdjacentNodes(v1); avr2.first != avr2.second; avr2.first++)
			{
				if (avr1.first == avr2.first)
					continue;
				if(graph.hasLink(*avr1.first, *avr2.first))
					numLinksInSubGraph += 1.0;		
			}
			
		if (numLinksInSubGraph < 1)
			continue; 
		else		
			ret += (numLinksInSubGraph)/(v1Deg * (v1Deg-1.0));
	}

	if (ret == 0)
		return 0;
		
	return ret/group.groupGenes.size();
}


float calculateClusteringCoeffForGroupOnly(const Graph &graph, const GeneGroup &group)
{
	Graph::Node v1;
	float numLinksInSubGraph, v1Deg, ret;
	ret = 0;
	map<Graph::Node, bool> adjGenesInGroup;
	map<Graph::Node, bool>::iterator it;
	
	for (int i = 0; i < (int)group.groupGenes.size(); i++)
	{
		v1 = getNodeById(graph, group.groupGenes[i]);
		v1Deg = graph.getNodeDegree(v1);
		numLinksInSubGraph = 0;
		
		if (v1Deg <= 1)
			continue;
		
		adjGenesInGroup.clear();
		for (Graph::adjacency_node_range_t avr1 = graph.getAdjacentNodes(v1); avr1.first != avr1.second; avr1.first++)
		{
			for (int j = 0; j < (int)group.groupGenes.size(); j++)
				if (graph.properties(*avr1.first).geneId == group.groupGenes[j])
				{
					adjGenesInGroup[*avr1.first] = true;
					break;
				}
		}
		
		for (Graph::adjacency_node_range_t avr1 = graph.getAdjacentNodes(v1); avr1.first != avr1.second; avr1.first++)
			for (it = adjGenesInGroup.begin(); it != adjGenesInGroup.end(); it++)
				if (graph.hasLink(*avr1.first, (*it).first))
					numLinksInSubGraph += 1.0;	
		
		if (numLinksInSubGraph < 1)
			continue; 
		else		
			ret += (numLinksInSubGraph)/(v1Deg * (v1Deg-1.0));			
	}
	if (ret == 0)
		return 0;
		
	return ret/group.groupGenes.size();
}


#if SPLIT_GROUPS
void splitGroups(vector<GeneGroup> &groups, map<string, vector<string> > &geneGroupMap)
{
	cout << "\n***Warning*** splitting groups!\n";
	geneGroupMap.clear();
	
	vector<GeneGroup> groupsNew;
	GeneGroup group1, group2;
	for (int i = 0; i < (int)groups.size(); i++)
	{
		group1.groupGenes.clear();
		group2.groupGenes.clear();
		
		random_shuffle(groups[i].groupGenes.begin(), groups[i].groupGenes.end());
		
		group1.groupId = groups[i].groupId;
		group2.groupId = groups[i].groupId;
		
		group1.groupId += "_1";
		group2.groupId += "_2";
		
		group1.groupSpe = groups[i].groupSpe;
		group2.groupSpe = groups[i].groupSpe;
		
		group1.groupSys = groups[i].groupSys;
		group2.groupSys = groups[i].groupSys;
		
		group1.groupDesc = groups[i].groupDesc;
		group2.groupDesc = groups[i].groupDesc;
		
		//printf("\n");
		for(int j = 0; j < (int)(groups[i].groupGenes.size()/2); j++)
		{
			geneGroupMap[groups[i].groupGenes[j]].push_back(group1.groupId);
			group1.groupGenes.push_back(groups[i].groupGenes[j]);
		//	printf("%s\t%s\n", group1.groupId.c_str(), groups[i].groupGenes[j].c_str());
		}
		for(int j = (int)(groups[i].groupGenes.size()/2); j < (int)groups[i].groupGenes.size(); j++)
		{
			geneGroupMap[groups[i].groupGenes[j]].push_back(group2.groupId);
			group2.groupGenes.push_back(groups[i].groupGenes[j]);
		//	printf("%s\t%s\n", group2.groupId.c_str(), groups[i].groupGenes[j].c_str());
		}
		/*printf("\n");
		for(int j = 0; j < groups[i].groupGenes.size(); j++)
			printf("%s\t%s\n", groups[i].groupId.c_str(), groups[i].groupGenes[j].c_str());
		*/	
		groupsNew.push_back(group1);
		groupsNew.push_back(group2);
	}

	groups.clear();
	for (int i = 0; i < (int)groupsNew.size(); i++)
		groups.push_back(groupsNew[i]);

}

#endif //if SPLIT_GROUPS

#if RAND_GROUPS
#define PI 3.14159265358979323846
//borrowed this code for randn_trig. it returns a random number for a
//normal distribution with mean mu and std sigma
double randn_trig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double dist, angle;
	
	//	If no deviate has been stored, the standard Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose a pair of uniformly distributed deviates, one for the
		//	distance and one for the angle, and perform transformations
		dist=sqrt( -2.0 * log(double(rand()) / double(RAND_MAX)) );
		angle=2.0 * PI * (double(rand()) / double(RAND_MAX));
		
		//	calculate and store first deviate and set flag
		storedDeviate=dist*cos(angle);
		deviateAvailable=true;
		
		//	calcaulate return second deviate
		return dist * sin(angle) * sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}

void randomizeGroups(Graph &origNet, vector<GeneGroup> &groups)
{
	cout << "*** WARNING: Randomizing groups ***" << endl;

	/*NORMAL DIST OF GROUPS CONNECTIVITY
	float mean = 500, var = 100;
	vector<string> allGenes;
	int maxDeg = 0;
	for (Graph::node_range_t vr = origNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		allGenes.push_back(origNet.properties(*vr.first).geneId);
		int deg = origNet.getNodeDegree(*vr.first);
		if (deg > maxDeg)
			maxDeg = deg;
	}
	
	int thisDeg;
	for (int i=0; i < groups.size(); i++)
	{
		thisDeg = randn_trig(mean, var);
		while(thisDeg <= 0)
			thisDeg = randn_trig(mean, var);
			
		groups[i].groupGenes.clear();
		for (int degCount = 0; degCount <= thisDeg;)
		{
			string which = allGenes[rand()%allGenes.size()];
			bool hasAlready = false;
			do {
				hasAlready = false;
				for (int k = 0; k < groups[i].groupGenes.size(); k++)
					if (groups[i].groupGenes[k] == which)
						{hasAlready = true;break;}
				if (hasAlready)
					which = allGenes[rand()%allGenes.size()];
			}while(hasAlready);
			
			int degWhich = origNet.getNodeDegree(getNodeById(origNet, which));
			while(degCount+degWhich > thisDeg+rand()%100)
			{
				which = allGenes[rand()%allGenes.size()];
				degWhich = origNet.getNodeDegree(getNodeById(origNet, which));
			}
		
			groups[i].groupGenes.push_back(which);
			degCount+=degWhich;
		}
		
	}*/
	
	
	/*CONSERVATIVE*/
	int conn1, randNum, highest = 0;
	map<int, vector<Graph::Node> > degMap;
	Graph::Node v1;
	
	for (Graph::node_range_t vr = origNet.getNodes(); vr.first != vr.second; vr.first++)
	{
		v1 = *vr.first;
		conn1 = DEGREE_BIN((float)origNet.getNodeDegree(v1));
		degMap[conn1].push_back(v1);
		if (conn1 > highest)
			highest = conn1;
	}
	
	for (int i = 0; i < (int)groups.size(); i++)
	{
		vector<string> origGenes = groups[i].groupGenes;
		
		for (int j = 0; j < (int)groups[i].groupGenes.size(); j++)
		{
			conn1 = DEGREE_BIN((float)origNet.getNodeDegree(getNodeById(origNet, groups[i].groupGenes[j])));
			int k, l;
			bool hasAlready;
			do {
				k = 0;
				hasAlready = false;
				randNum = rand()%(degMap[conn1].size());
				for (k = 0; k <= j; k++)
				{
					if (groups[i].groupGenes[k] == origNet.properties(degMap[conn1][randNum]).geneId)
					{
						hasAlready = true;
						break;
					}
				}
				
				for (l = 0; l < (int)origGenes.size(); l++)
				{
					if (origGenes[l] == origNet.properties(degMap[conn1][randNum]).geneId)
					{
						hasAlready = true;
						break;
					}
				}
				
				if (hasAlready && (k > j || l == origGenes.size()))
				{
					if (rand()%2 == 0)
					{
						while(conn1 > 1){
							conn1--;
							if (keyInMap(conn1, degMap))
								break;
						}
						//cout << "minus" << conn1 <<endl;
					}
					else
					{
						while(conn1 < highest){
							conn1++;
							if (keyInMap(conn1, degMap))
								break;
						}
						//cout << "plus" << conn1 <<endl;
					}
				}
				
			}while(hasAlready  && conn1 > 0 && conn1 <= highest);
			//cout << "Done" << endl;
			if (hasAlready && k == (int)groups[i].groupGenes.size())
			{
				//errGrps[groups[i].groupId] = groups[i].groupId;
				cout << groups[i].groupId << endl;
			}
			groups[i].groupGenes[j] = origNet.properties(degMap[conn1][randNum]).geneId;
		}
	}
}
#endif //if RAND_GROUPS




