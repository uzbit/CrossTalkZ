/*

This file contains the entry point of CrossTalkZ.
It is called main-test because it contains many debug 
and test modes that can be enabled via #defines in defines.h
One should be able to rename this to main.cpp and compile without
issues. 

Written by Ted McCormack.

*/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>

#include <boost/program_options.hpp>
#include <algorithm>
#include <iterator>
#include <exception>

#include "crosstalkz.h"
#include "boostgraphio.h"

using namespace std;
using namespace boost;
using namespace program_options;
	
string NetworkFile;
string GroupsFile;
string GroupsFile1, GroupsFile2;
string AnalysisResultsFile;
string RandomGraphFile;
string ResultFileFormat("crosstalkz_%s.csv");
string InfoFileFormat("crosstalkz_%s.info");

bool userSpecifiedOutFile = false;
bool writeRandomGraphOnly = false;

void parseArgs(int argc, char *argv[]);
void printInfos(ostream &os);

#define BUF_SIZE 5000

#if VERBOSE
	#define COUNT_EDGES		1
#endif

/****************PROGRAM START****************/
int main(int argc, char *argv[])
{
	srand(time(NULL));

	try {
	
		Graph origNetwork;
		Graph randNetwork;
		BoostGraphIO bgio;
		vector<GeneGroup> groups;
		vector<GeneGroup> groups1, groups2;
		map<string, vector<string> > geneToGroupMap;
		map<string, vector<string> > geneToGroupMap1, geneToGroupMap2;
		map<string, Stats > groupStatistics;
		map<int, vector<Record> > degToRecordsMap;
		stringstream infoString;
		
		parseArgs(argc, argv);
	
		bgio.setGraph(&origNetwork);
		bgio.readGraph(NetworkFile);
		
		infoString << endl << "----NETWORK STATISTICS----" << endl;
		infoString << "Final number of unique nodes in the network: " << origNetwork.getNodeCount() << endl;
		infoString << "Final number of links in the network: " << origNetwork.getLinkCount() << endl;
		
			
		//cout << "orig r: " << calculateRfromNetwork(origNetwork) << endl;
		//bgio.writeTSVGraph(string("origNetwork.csv"));
		
		randNetwork = origNetwork; //dont need copy function yet, but must copy like this once and only once.	
		generateMaps(origNetwork, randNetwork, degToRecordsMap);
	
		origNetwork.id = 0;
		randNetwork.id = 1;
		
	#if VERBOSE	
		cout << endl << "Printing binning counts:" << endl;
		map<int, vector<Record> >::iterator it;
		for ( it=degToRecordsMap.begin() ; it != degToRecordsMap.end(); it++ )
			cout << "bin #: " << it->first << " count: " << it->second.size() << endl;
	#endif
		
		
		if (writeRandomGraphOnly)
		{
			switch(methodFlag)
			{
				case METHOD_ASSIGN:
				{
					copyOrigToRand(origNetwork, randNetwork);
					generateRandomNetworkAssignment(origNetwork, randNetwork);
					break;
				}
				case METHOD_ASSIGN_SECOND:
				{
					copyOrigToRand(origNetwork, randNetwork);
					generateRandomNetworkSecondOrder(origNetwork, randNetwork, degToRecordsMap);
					break;
				}
				case METHOD_LABELSWAP:
				{
					generateRandomNetworkLabelSwap(origNetwork, randNetwork, degToRecordsMap);
					break;
				}
				case METHOD_LINKSWAP:
				{
					copyOrigToRand(origNetwork, randNetwork);
					generateRandomNetworkLinkSwap(origNetwork, randNetwork);
					break;
				}	
				default:
					break;
			}
			cout << "\nWriting random graph to " << RandomGraphFile << endl;
			bgio.setGraph(&randNetwork);
			bgio.writeTSVGraph(RandomGraphFile);
			exit(0);
		}
	
		if (allVsall)
		{
			infoString << endl << "----GROUP STATISTICS----" << endl;
			readGeneGroups(origNetwork, groups, geneToGroupMap, GroupsFile, infoString);
			cout << infoString.str();
		}
		else
		{
			infoString << "----GROUP A STATISTICS----" << endl;
			readGeneGroups(origNetwork, groups1, geneToGroupMap1, GroupsFile1, infoString);
			infoString << "----GROUP B STATISTICS----" << endl;
			readGeneGroups(origNetwork, groups2, geneToGroupMap2, GroupsFile2, infoString);
			cout << endl <<  infoString.str();
		}
	
	#if SPLIT_GROUPS	
		splitGroups(groups, geneToGroupMap);
	#endif	
		
		/*
		for (int i = 0; i < (int)groups.size(); i++)
		for (int j = 0; j < (int)groups.size(); j++)
			cout << "Groups " << i << " " << j << " have cluster coeff: " <<  calculateClusteringCoeffForTwoGroups(origNetwork, groups[i], groups[j]) << endl;
	
		exit(0);
		*/
		
		if (origNetwork.getNodeCount() && (groups.size() || (groups1.size() && groups2.size())))
		{		
			string groupsVsStr;
			clock_t iterStart, start = clock();
			
			groupStatistics.clear();
			if (allVsall)
			{		
				for (int i = 0; i < (int)groups.size(); i++)
					for (int j = 0; j <= i; j++)
					{
						groupsVsStr = groups[i].groupId + "_vs_" + groups[j].groupId;
						groupStatistics[groupsVsStr].linkCount.clear();
						groupStatistics[groupsVsStr].clusteringCoeff.clear();
					}
			}
			else
			{
				for (int i = 0; i < (int)groups1.size(); i++)
					for (int j = 0; j < (int)groups2.size(); j++)
					{
						groupsVsStr = groups1[i].groupId + "_vs_" +  groups2[j].groupId;
						groupStatistics[groupsVsStr].linkCount.clear();
						groupStatistics[groupsVsStr].clusteringCoeff.clear();
					}
			}
				
			
/////////////////////////////////////////////////////////////////////////////////////////
#if RAND_GROUPS
map<string, float> zscores;
vector<vector<string> > backupGroups;
for (int i = 0; i < (int)groups.size(); i++)
	backupGroups.push_back(groups[i].groupGenes);

for(int j = 90; j < 100; j++)
{
	for (int i = 0; i < (int)groups.size(); i++)
	{
		groups[i].groupGenes.clear();
		groups[i].groupGenes = backupGroups[i];
	}
			
	cout << "Group randomization "<< j << endl;
	randomizeGroups(origNetwork, groups);
		
	float avg=0;
	int degCount = 0;	

	ofstream file;
	char str1[BUF_SIZE];
	sprintf(str1, "/Users/uzbit/Documents/Biophysics/project/CrossTalkZ/branches/trunk-smetric/analysis/negtest/groups/group_set_%d.hsa", j);
	file.open(str1);
	for (int i = 0; i < (int)groups.size(); i++)
	for (int k = 0; k < (int)groups[i].groupGenes.size(); k++)
		file << groups[i].groupGenes[k] << "\t" << groups[i].groupId << "\tNONE\tSPE\tDESC" <<endl;
		
	continue; //skip network randomization when generating random groups
#endif		
//////////////////////////////////////////////////////////////////////////////////////////
	
	#if COUNT_EDGES
		vector<int> edgeCounts;
	#endif			
		for (int i = 0; i < numSimIter; i++)
		{
			cout << "\nIteration "<< i+1 << " out of " << numSimIter << " ..." << endl;
			
			iterStart = clock();
			switch(methodFlag)
			{
				case METHOD_ASSIGN:
				{
					if (!generateRandomNetworkAssignment(origNetwork, randNetwork))
						copyOrigToRand(origNetwork, randNetwork);
					else
						printf("Randomized %d links between %d nodes in %f seconds.\n", randNetwork.getLinkCount(), randNetwork.getNodeCount(), (clock()-iterStart)/(CLOCKS_PER_SEC+0.0)); 
					break;
				}
				case METHOD_ASSIGN_SECOND:
				{
					if (!generateRandomNetworkSecondOrder(origNetwork, randNetwork, degToRecordsMap))
						copyOrigToRand(origNetwork, randNetwork);
					else
						printf("Randomized %d links between %d nodes in %f seconds.\n", randNetwork.getLinkCount(), randNetwork.getNodeCount(), (clock()-iterStart)/(CLOCKS_PER_SEC+0.0)); 
					break;
				}
				case METHOD_LABELSWAP:
				{
					generateRandomNetworkLabelSwap(origNetwork, randNetwork, degToRecordsMap);
					printf("Randomized labels for %d nodes in %f seconds.\n", randNetwork.getNodeCount(), (clock()-iterStart)/(CLOCKS_PER_SEC+0.0)); 
					break;
				}
				case METHOD_LINKSWAP:
				{
					copyOrigToRand(origNetwork, randNetwork);
					int swapped = generateRandomNetworkLinkSwap(origNetwork, randNetwork);
					printf("Swapped %d of %d links between %d nodes in %f seconds.\n", swapped, randNetwork.getLinkCount(), randNetwork.getNodeCount(), (clock()-iterStart)/(CLOCKS_PER_SEC+0.0)); 
					break;
				}
				
				default:
					break;
			}
			
	#if COUNT_EDGES
			Graph::Node v1, v2;
			int countLinksCommon = 0;
			for (Graph::link_range_t er = randNetwork.getLinks(); er.first != er.second; er.first++)
			{	
				randNetwork.getNodesByLink(*er.first, v1, v2);
				if (origNetwork.hasLink(geneVertMap[randNetwork.properties(v1).geneId][0], geneVertMap[randNetwork.properties(v2).geneId][0]))
					countLinksCommon ++; 	
			}
			edgeCounts.push_back(countLinksCommon);
			cout << "Networks have " << countLinksCommon << " common edges" << endl;
	#endif			
			if (allVsall)
				countLinksForGroupsAll(randNetwork, groups, groupStatistics, geneToGroupMap);
			else
				countLinksForGroups12(randNetwork, groups1, groups2, groupStatistics, geneToGroupMap1, geneToGroupMap2);	
		}
		cout << "\nFinished in " << (clock()-start+0.0)/CLOCKS_PER_SEC << " seconds." <<endl;
	
//////////////////////////////////////////////////////////////////////////////////////////
#if RAND_GROUPS
		char str[BUF_SIZE];
		sprintf(str, ResultFileFormat.c_str(), j);
		AnalysisResultsFile = str;
		calculateAndWriteResultsAll(origNetwork, groups, groupStatistics, geneToGroupMap, AnalysisResultsFile);
}
exit(0);
#endif
//////////////////////////////////////////////////////////////////////////////////////////
		//Write out the result and info files.
		char *str;
		char timebuf[100];
		time_t rawtime;	
		struct tm *timeinfo;
		time(&rawtime);
		timeinfo = localtime(&rawtime);
		strftime(timebuf,100,"%Y%m%d%H%M",timeinfo);
		
		if (userSpecifiedOutFile)
			AnalysisResultsFile = ResultFileFormat;
		else
		{
			str = new char[ResultFileFormat.size()+sizeof(timebuf)];
			sprintf(str, ResultFileFormat.c_str(), timebuf);
			AnalysisResultsFile = str;
			delete [] str;
		}
			
		if (allVsall)
			calculateAndWriteResultsAll(origNetwork, groups, groupStatistics, geneToGroupMap, AnalysisResultsFile);
		else
			calculateAndWriteResults12(origNetwork, groups1, groups2, groupStatistics, geneToGroupMap1, geneToGroupMap2, AnalysisResultsFile);
			
		if (userSpecifiedOutFile)
		{
			InfoFileFormat = ResultFileFormat;
			InfoFileFormat += ".info";
			str = new char[InfoFileFormat.size()];
			sprintf(str, "%s", InfoFileFormat.c_str());	
		}
		else
		{	
			str = new char[InfoFileFormat.size()+sizeof(timebuf)];
			sprintf(str, InfoFileFormat.c_str(), timebuf);	
		}
		ofstream file;
		file.open(str);
		printInfos(file);
		file << endl;
		file << infoString.str();
		file.close();
		delete [] str;
	
	#if DEBUG && VERBOSE
		float avg=0;
		for (int i=0;i<(int)rVals.size();i++)
			avg+=rVals[i];
		avg/=rVals.size();
		float std=0;
		for (int i=0;i<(int)rVals.size();i++)
			std+=powf(rVals[i]-avg,2.0);
		std/=rVals.size();
		std=sqrt(std);	
		printf("r avg = %f, std = %f. over %d iterations\n", avg, std, (int)rVals.size()); 
		
		avg = 0;
		for (int i=0;i<(int)smetricRatio.size();i++)
			avg+=smetricRatio[i];
		avg/=smetricRatio.size();
		std = 0;
		for (int i=0;i<(int)smetricRatio.size();i++)
			std+=powf(smetricRatio[i]-avg,2.0);
		std/=smetricRatio.size();
		std=sqrt(std);	
		printf("smetric ratio avg = %f, std = %f. over %d iterations\n", avg, std, (int)smetricRatio.size()); 
	
	#if COUNT_EDGES	
		avg = 0;
		for (int i=0;i<(int)edgeCounts.size();i++)
			avg+=edgeCounts[i];
		avg/=edgeCounts.size();
		std = 0;
		for (int i=0;i<(int)edgeCounts.size();i++)
			std+=powf(edgeCounts[i]-avg,2.0);
		std/=edgeCounts.size();
		std=sqrt(std);	
		printf("edge count avg = %f, std = %f. over %d iterations\n", avg, std, (int)smetricRatio.size()); 
	#endif
		
	#endif
			
		} //end if network has valid nodes and groups
		else
		{
			if (!origNetwork.getNodeCount())
				cout << "Network did not contain any vertices." << endl;
			else
				cout << "No valid groups." << endl; 
		}
	
	}
	catch (std::bad_alloc& ba){
		cout << endl;
		cout << "Caught bad_alloc exception with message: " << ba.what() << endl << endl;
		cout << "Error allocating memory. It is likely that the process ran out of memory." << endl;
		cout << "Reduce the size of input network or group sets, or get more memory." << endl;
	}
	catch (std::exception& e){
		cout << endl;
		cout << "Standard exception: " << e.what() << endl;
	}
   
	return 0;
}

void parseArgs(int argc, char *argv[])
{
	try{
		NetworkFile = "";
		GroupsFile = "";
		options_description desc("Allowed options", 300);
		desc.add_options()
			("help,h", "Produces this message.")
			("network,n", value< string >(&NetworkFile), 
				  "Path to a network file. Required.")
			("group,g", value< string >(&GroupsFile), "Path to a group file. Results are comparisons between all possible\ngroup pair combinations within this file.")
			("groupA,a", value< string >(&GroupsFile1), "Path to a group file. Results are comparisons between all possible\ngroup pair combinations between groupA and groupB. Requires groupB file.")
			("groupB,b", value< string >(&GroupsFile2), "Path to a group file. Results are comparisons between all possible\ngroup pair combinations between groupA and groupB. Requires groupA file.")
			("cutoff,c", value<float>(&cutoffScore), 
				  "Lowest link weight to include in network. If not specified, all links are included.")
			//("clusteringCoeff,C", value<bool>(&doClusteringCoeff)->default_value(doClusteringCoeff),
			//	  "Also calculate and output the statistics using clustering coefficient for groups. (EXPERIMENTAL)")
			("method,d", value<int>(&methodFlag)->default_value(methodFlag),
				  "Method 0: Link Permutation, swap links between nodes.\nMethod 1: Link Assignment, assign links uniformly randomly, conserve degree.\nMethod 2: Link Assignment + Second-order, same as 1 but attempt to conserve second-order properties also.\nMethod 3: Node Permutation, swap node labels only.")
			("iter,i", value<int>(&numSimIter)->default_value(numSimIter),
				  "Number of network randomizations.")
			("mode,m", value<int>(&modeFlag)->default_value(modeFlag),
				  "Mode 0: Link isn't counted if either gene belongs to both groups.\nMode 1: Link isn't counted if both genes belong to both groups.")
			("outputFile,o", value< string > (&ResultFileFormat), "User specified results file.")
			("phyper,p", value<bool>(&doHyper)->default_value(doHyper),
				  "Also calculate and write out the hypergeometric probaility of the overlap\nbetween each group pair for gene set enrichment analysis.")
			("writeGraph,w", value< string > (&RandomGraphFile), "Randomize original graph once and output graph to specified file.")
			("minGenes,x", value<int>(&minimumGenesForGroup)->default_value(minimumGenesForGroup),
				  "Set the lower bound on the minimum number of genes a group should have to be included in the analysis.")
			
			;

		positional_options_description p;
		
		variables_map vm;
		store(command_line_parser(argc, argv).
				  options(desc).positional(p).run(), vm);
		notify(vm);
	
		if (vm.count("help") || argc < 5
			|| !vm.count("network")
			|| (vm.count("group") && (vm.count("groupA") || vm.count("groupB")))
			|| (!vm.count("group") && !vm.count("groupA") && !vm.count("groupB"))
			|| (!vm.count("group") && (vm.count("groupA") && !vm.count("groupB")))
			|| (!vm.count("group") && (!vm.count("groupA") && vm.count("groupB"))) ) 
		{
			cout << "CrossTalkZ " << VERSION << endl;;
			cout << "Usage: "<<argv[0]<<" [options] -n NETWORK_FILE [-g GROUP_FILE] or [-a GROUP_A_FILE -b GROUP_B_FILE]" << endl;
			cout << desc;
			exit(1);
		}
		
		/*
		if (vm.count("clusteringCoeff"))
			doClusteringCoeff = true;
		else
			doClusteringCoeff = false;
		*/
		if (vm.count("outputFile"))
			userSpecifiedOutFile = true;
		else
			userSpecifiedOutFile = false;
			
		if (vm.count("writeGraph"))
			writeRandomGraphOnly = true;
		else
			writeRandomGraphOnly = false;
		
		if (vm.count("group"))
			allVsall = true;
		else
			allVsall = false;
		
		if (vm.count("cutoff"))
			useCutoff = true;
		else
			useCutoff = false;	
		
		if (modeFlag != MODE_0 && modeFlag != MODE_1)
		{
			cout << "Invalid Mode: "<< modeFlag << endl;
			exit(1);		
		}
		
		if (methodFlag != METHOD_DEFAULT && methodFlag != METHOD_ASSIGN_SECOND
			&& methodFlag != METHOD_ASSIGN && methodFlag != METHOD_LABELSWAP
			&& methodFlag != METHOD_LINKSWAP)
		{
			cout << "Invalid Method: "<< methodFlag << endl;
			exit(1);		
		}
		
		//cheap way to test if valid files
		ifstream file(NetworkFile.c_str());
		if (!file)
		{
			cout << "Error opening "<< NetworkFile << endl;
			exit(1);
		}
		file.close();
		if (allVsall)
		{
			file.open(GroupsFile.c_str());
			if (!file)
			{
				cout << "Error opening "<< GroupsFile << endl;
				exit(1);
			}
			file.close();
		}
		else
		{
			file.open(GroupsFile1.c_str());
			if (!file)
			{
				cout << "Error opening "<< GroupsFile1 << endl;
				exit(1);
			}
			file.close();
			file.open(GroupsFile2.c_str());
			if (!file)
			{
				cout << "Error opening "<< GroupsFile2 << endl;
				exit(1);
			}
			file.close();
		}
	}
	catch(std::exception& e)
	{
		cout << e.what() << endl;
		exit(1);
	} 
	
	printInfos(cout);
}

void printInfos(ostream &os)
{
	os << endl << "CrossTalkZ version: " << VERSION << "\nUsing the following parameters:"<< endl <<endl;
	os << "Network file:\t\t\t" << NetworkFile << endl;
	if (allVsall)
		os << "Group file:\t\t\t" << GroupsFile << endl;
	else
	{
		os << "Group A file:\t\t\t" << GroupsFile1 << endl;
		os << "Group B file:\t\t\t" << GroupsFile2 << endl;
	}
	if (userSpecifiedOutFile)
		os << "Result file:\t\t\t" << ResultFileFormat << endl;
	
	if (writeRandomGraphOnly)
		os << "Random network file:\t\t" << RandomGraphFile << endl;
	
	if (useCutoff)
		os << "Link cutoff:\t\t\t" << cutoffScore << endl;
	else
		os << "Link cutoff:\t\t\tnone" << endl;
		
	os << "Iterations:\t\t\t" << numSimIter << endl;
	os << "Link counting mode:\t\t" << modeFlag << endl;
	os << "Randomization method:\t\t" << getMethodString(methodFlag) << endl;
	//os << "Also use clustering coeff:\t" << doClusteringCoeff << endl;
	os << "Minimum genes for group:\t" << minimumGenesForGroup<< endl;
	
	os << endl;
}


