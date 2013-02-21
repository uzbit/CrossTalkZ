/*

This file contains function declarations for loading gene groups, 
randomizing Graphs, calculating statistics, and functions for the
various tests that are presented in the paper: 
"Statistical assessment of gene group crosstalk enrichment in networks", 
2013, PLoS ONE.

Written by Theodore McCormack.

*/

#ifndef __CROSSTALKZ_H__
#define __CROSSTALKZ_H__

#include "types.h"
#include "defines.h"

using namespace std;
using namespace boost;

/****************GLOBAL VARIABLES AND TYPES****************/

extern int numSimIter;
extern float cutoffScore;
extern bool useCutoff;
extern int modeFlag;
extern bool allVsall;
extern int methodFlag;
extern int minimumGenesForGroup;
extern bool doClusteringCoeff;
extern bool doHyper;

extern map<string, vector<Graph::Node> > geneVertMap; 
		

#if VERBOSE
extern vector<float> rVals;
extern vector<float> smetricRatio;
#endif

/****************FUNCTION DECLARATIONS****************/

//readGeneGroups:
//	Loads the gene group table from GeneGroupTableFile and makes a mapping of genes to 
//	a list of group id's that the each gene belongs to.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	groups: a vector of GeneGroups in any state to store group information loaded from GeneGroupTableFile
//	geneGroupMap: a map in any state used for mapping a gene to a vector of group id's that the gene belongs to.
//  ss: a stringstream to print group statistics into, used to print to stdout, info file in main.cpp
//	path: a string containing the file to load group information from
void readGeneGroups(const Graph &origNet, 
					vector<GeneGroup> &groups, 
					map<string, vector<string> > &geneGroupMap, 
					string path,
					stringstream &ss);

//generateRandomNetworkLinkSwap:
//	Swap links as suggested by maslov and sneppen: link pair (a, b) and (c, d) become
//	(a, c) and (b, d) or (a, d) and (c, b)  
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph that is a copy of origNet and is used to store the new randomized version of origNet
// 	returns the number of swaps performed
int generateRandomNetworkLinkSwap(const Graph &origNet, Graph &randNet);

//generateRandomNetworkLabelSwap:
//	Permutates node labels that fall into the same ln(deg) bin.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph in any state that is used to store the new randomized version of origNet
//  degRecordsMap: a map from degree bin to records with that degree bin (from generateMaps below)
//	returns true if randomization conserved connectivities of the original else false
bool generateRandomNetworkLabelSwap(const Graph &origNet,
								Graph &randNet,
								map<int, vector<Record> > &degRecordsMap);

//generateRandomNetworkSecondOrder:
//	Best effort randomization of the original network. Attempts to conserve second-order assortativity.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph in any state that is used to store the new randomized version of origNet
//  degRecordsMap: a map from degree bin to records with that degree bin (from generateMaps below)
//	returns true if randomization conserved connectivities of the original else false
bool generateRandomNetworkSecondOrder(const Graph &origNet, 
						   Graph &randNet,
						   map<int, vector<Record> > &degRecordsMap);

//generateRandomNetworkAssignment:
//	Best effort randomization of the original network.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph in any state that is used to store the new randomized version of origNet
//	returns true if randomization conserved connectivities of the original else false
bool generateRandomNetworkAssignment(const Graph &origNet, 
						   Graph &randNet);

//countLinksForGroupsAll:
//	Calculates links in Graph randNet between groups.
//	randNet: a randomized Graph with validated connectivities
//	groups: a vector containing the list of groups from readGeneGroups
//	groupStats: a map from group1_vs_group2 string to a vector containing a float for each iteration 
//				corresponding to the number of links between group1 and group2
//	geneGroupMap: a map from a gene string to a vector containing group strings that the gene belongs to			
void countLinksForGroupsAll(Graph &randNet, 
							 vector<GeneGroup> &groups, 
							 map<string, Stats > &groupStats,
							 map<string, vector<string> > &geneGroupMap);

//countLinksForGroups12:
//	Counts Links in Graph randNet between groups 1 and 2.
//	randNet: a randomized Graph with validated connectivities
//	groups: a vector containing the list of groups from readGeneGroups
//	groupStats: a map from group1_vs_group2 string to a vector containing a float for each iteration 
//				corresponding to the number of links between group1 and group2
//	geneGroupMap: a map from a gene string to a vector containing group strings that the gene belongs to			
void countLinksForGroups12(Graph &randNet, 
							 vector<GeneGroup> &groups1, 
							 vector<GeneGroup> &groups2, 
							 map<string, Stats > &groupStats,
							 map<string, vector<string> > &geneGroupMap1,
							 map<string, vector<string> > &geneGroupMap2);

//void countLinks(Graph &origNet, Graph &randNet, Graph &resultsNet, map<string, vector<Graph::Node> > &geneVertMap);
//void writeConnectivityMatrix(Graph &network);

//calculateAndWriteResultsAll:
//	Calculates the observed and expected links between groups, z score, p value, etc 
//	and writes it to the file pointed to in path.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	groups: a vector containing the list of groups from readGeneGroups
//	groupStats: a map from group1_vs_group2 string to a vector containing a float for each iteration 
//				corresponding to the number of links between group1 and group2
//	geneGroupMap: a map from a gene string to a vector containing group strings that the gene belongs to			
//	path: a string containing the file to write group information to
void calculateAndWriteResultsAll(Graph &origNet,  
				  vector<GeneGroup> &groups, 
				  map<string, Stats > &groupStats,	
				  map<string, vector<string> > &geneGroupMap,
				  string path); 

//calculateAndWriteResults12:
//	Calculates the observed and expected links between groups 1 and 2, z score, p value, etc 
//	and writes it to the file pointed to in path.
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	groups1: a vector containing the list of groups from readGeneGroups
//	groups2: a vector containing the list of groups from readGeneGroups
//	groupStats: a map from group1_vs_group2 string to a vector containing a float for each iteration 
//				corresponding to the number of links between group1 and group2
//	geneGroupMap1: a map from a gene string to a vector containing group strings that the gene belongs to			
//	geneGroupMap2: a map from a gene string to a vector containing group strings that the gene belongs to			
//	path: a string containing the file to write group information to
void calculateAndWriteResults12(Graph &origNet,  
				  vector<GeneGroup> &groups1, 
				  vector<GeneGroup> &groups2, 
				  map<string, Stats > &groupStats,	
				  map<string, vector<string> > &geneGroupMap1,
				  map<string, vector<string> > &geneGroupMap2,
				  string path); 
					
//generateMaps:
//	Create maps that can be used for quickly accessing  elements or records given 
//	a  id string or degree bin respectively. (this function creates geneVertMap and degRecordsMap)
//
//	origNet: a Graph that contains the network loaded from boostgraphio
//	randNet: a Graph in any state that is used to store the new randomized version of origNet
//  degRecordsMap: a map from degree bin to records with that degree bin (from generateMaps below)
void generateMaps(Graph &origNet, 
				  Graph &randNet, 
				  map<int, vector<Record> > &degRecordsMap); 	



//utility functions
Graph::Node getNodeById(const Graph &g, const string &Id);
void printNetwork(const Graph &network);
string getMethodString(int m);
int getTotalInputUniqueGeneCount(string path1, string path2);
void copyOrigToRand( Graph &origNet, Graph &randNet);//Need special copy function to preserve  and link data.


//statistics functions
float calculateRfromNetwork(const Graph &g);
long long calculateSmetricNetwork(const Graph &g);
int calculateSmetricNode(const Graph &g, Graph::Node v);
float calculateClusteringCoeffForTwoGroups(const Graph &graph, const GeneGroup &group1, const GeneGroup &group2);
float calculateClusteringCoeffForGroup(const Graph &graph, const GeneGroup &group);
float calculateClusteringCoeffForGroupOnly(const Graph &graph, const GeneGroup &group);




#if RAND_GROUPS
void randomizeGroups(Graph &origNet,
					vector<GeneGroup> &groups); 
#endif

#if SPLIT_GROUPS
void splitGroups(vector<GeneGroup> &groups, map<string, vector<string> > &geneGroupMap);
#endif


#endif


