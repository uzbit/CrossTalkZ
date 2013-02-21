/*

This file contains BoostGraphIO class definition that is used to load
and write networks. It allows the user to easily load XGMML, 
simple TSV, and funcoup TSV graphs. It currently only writes 
simple TSV graphs. 

Written by Ted McCormack.

*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <libxml/tree.h>
#include <boost/algorithm/string.hpp>

#include "boostgraphio.h"
#include "defines.h"
#include "types.h"

using namespace std;

#define SIMPLE_TSV		3
#define FUNCOUP_TSV		7

void err(void *ctx, const char *msg, ...) 
{
}
 
Graph& BoostGraphIO::readGraph(string path)
{
	ifstream file;
	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
	file.close();
	
	xmlGenericErrorFunc handler = (xmlGenericErrorFunc)err;
	initGenericErrorDefaultFunc(&handler);
	
	//determine if path is tsv or xgmml
	xmlDocPtr doc = xmlParseFile(path.c_str());
	if (!doc)
	{
		return readTSVGraph(path);	
	}
	else
	{
		xmlFreeDoc(doc);
		return readXGMMLGraph(path);
	}
}

Graph& BoostGraphIO::readXGMMLGraph(string path)
{
	//assume the user sends a path that is xgmml
	xmlDocPtr doc;
	xmlNodePtr nodeLevel1, nodeLevel2;
	NodeProperties vp;
	LinkProperties link;
	string str1, str2, weight;
	map<string, Graph::Node> idVertMap;
	//clock_t start = clock();
	
	doc = xmlParseFile(path.c_str());
	if (!doc)
		{xmlFreeDoc(doc);cout << "Invalid XML file "<< path << endl;exit(1);}

	graphPtr->Clear();

	nodeLevel1 = doc->children;
	while (strcmp((const char*)nodeLevel1->name, "graph") && nodeLevel1 != NULL)
		nodeLevel1 = nodeLevel1->next;
	
	if (nodeLevel1 == NULL)
	{xmlFreeDoc(doc);cout << "Invalid XGMML format in "<< path <<endl;exit(1);}
	
	for(nodeLevel2 = nodeLevel1->children;
		nodeLevel2 != NULL;
		nodeLevel2 = nodeLevel2->next)
	{
		if (!strcmp((const char*)nodeLevel2->name, "node"))
		{
			str1 = (char*)xmlGetProp(nodeLevel2, (xmlChar *)"id");
			if (idVertMap[str1] == NULL)
			{
				vp.geneId = (char*)xmlGetProp(nodeLevel2, (xmlChar *)"label");
				to_upper(vp.geneId);
				idVertMap[str1] = graphPtr->AddNode(vp);
			}
		}
		if (!strcmp((const char*)nodeLevel2->name, "edge"))
		{
			str1 = (char*)xmlGetProp(nodeLevel2, (xmlChar *)"source");
			str2 = (char*)xmlGetProp(nodeLevel2, (xmlChar *)"target");
			link.weight = cutoffScore;
			if ((char*)xmlGetProp(nodeLevel2, (xmlChar *)"weight") != NULL && useCutoff)
			{
				weight = (char*)xmlGetProp(nodeLevel2, (xmlChar *)"weight");
				link.weight = atof(weight.c_str());
				if (link.weight >= cutoffScore)
					graphPtr->AddLink(idVertMap[str1], idVertMap[str2], link);	
			}
			else{	
				graphPtr->AddLink(idVertMap[str1], idVertMap[str2], link);
			}
		}
	}
	
	xmlFreeDoc(doc);
	
	//clear links with 0 connections
	for (Graph::node_range_t vr = graphPtr->getNodes(); vr.first != vr.second; vr.first++)
	{
		if (graphPtr->getNodeDegree(*vr.first) == 0)
			graphPtr->RemoveNode(*vr.first);
	}
/*	
#if USE_BIDIRECTIONAL
	printf("Loaded %d links between %d nodes in %f seconds.\n", graphPtr->getLinkCount()/2, graphPtr->getNodeCount(), ((clock()-start)+0.0)/CLOCKS_PER_SEC);	
#else 
	printf("Loaded %d links between %d nodes in %f seconds.\n", graphPtr->getLinkCount(), graphPtr->getNodeCount(), ((clock()-start)+0.0)/CLOCKS_PER_SEC);	
#endif
*/

	return (*graphPtr);
}

Graph& BoostGraphIO::readTSVGraph(string path)
{
	ifstream file;
	string line;
	vector<string> lineVals;
	
	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
	getline(file, line);
	file.close();
	
	split(lineVals, line, is_any_of("\t "));
	//cout << "size" << lineVals.size()<<endl;
	if (lineVals.size() > SIMPLE_TSV)
	{
		return readFunCoupTSVGraph(path);
	}
	else if (lineVals.size() <= SIMPLE_TSV && lineVals.size() > 1)
	{
		return readSimpleTSVGraph(path);
	}
	else
	{
		cout << "Invalid TSV format in "<< path <<endl;
		exit(1);
	}
	return graph;
}

Graph& BoostGraphIO::readSimpleTSVGraph(string path)
{
	return __readTSVGraph__(path, SIMPLE_TSV);
}

Graph& BoostGraphIO::readFunCoupTSVGraph(string path)
{
	return __readTSVGraph__(path, FUNCOUP_TSV);
}

Graph& BoostGraphIO::__readTSVGraph__(string path, int flag)
 {
	//assume the user sends a path that is a tsv
	ifstream file;
	string line, cell, header;
	NodeProperties vp;
	LinkProperties link;
	Graph::Node v1, v2;
	int countCell = 0;
	float maxScore = 0;
	bool hasFirst, hasSecond;
	string first, second;
	stringstream lineStream;
	map<string, Graph::Node> geneVertMap;
	vector<string> lineVals;
	
	//clock_t start = clock();

	file.open(path.c_str());
	if (!file)
	{
		cout << "Error in opening "<< path << endl;
		exit(1);
	}
	
	if (useCutoff)
		cout << "Reading network from "<< path << " using link weight cutoff >= "<< cutoffScore <<" ..." << endl;
	else
		cout << "Reading network from "<< path << " ..." << endl;
	
	if(flag == FUNCOUP_TSV)getline(file, header);

	//could do some file type validation here...
	graphPtr->Clear();
		
#if DO_LINE_LIMIT
	int countLines = 0;
	while(getline(file, line) && countLines < MAX_LINES)
	{
		countLines++;
#else
	while(getline(file, line))
	{
#endif
		switch(flag)
		{
			case SIMPLE_TSV:
				split(lineVals, line, is_any_of("\t "));
				switch(lineVals.size())
				{
					case 2: //must be protein1\tprotein2
						first = lineVals[0];
						second = lineVals[1];
						maxScore = cutoffScore+100.0; //no weight on link so add a little to pass the test below
						break;
					case 3: //must be protein1\tprotein2\tscore
						first = lineVals[0];
						second = lineVals[1];
						maxScore = atof(lineVals[2].c_str());
						break;
					default:
						//cout << "Invalid TSV format in "<<path<<endl;
						//exit(1);
						continue;
						break;
				}
			
				break;
			case FUNCOUP_TSV:		
				//parse the line
				lineStream.str(line);
				countCell = 0;
				while(getline(lineStream, cell, '\t') && countCell < 7)
				{
					if (countCell == FUNCOUP_MAX_SCORE)
						maxScore = atof(cell.c_str());
					else if (countCell == FUNCOUP_PROTEIN1)
						first = cell;
					else if (countCell == FUNCOUP_PROTEIN2)
						second = cell;
					countCell++;
				}
				break;
			default:
				break;
		}
			
		if (!first.size() || !second.size())
			continue;
			
		if ((maxScore >= cutoffScore && useCutoff) || !useCutoff)
		{
			to_upper(first);
			to_upper(second);
			first.erase(first.find_last_not_of(" \n\r\t")+1);
			second.erase(second.find_last_not_of(" \n\r\t")+1);
		
			hasFirst = (geneVertMap[first] != NULL);
			hasSecond = (geneVertMap[second] != NULL);

			//The following logic is based upon the assumption that the 
			//results from FunCoup are not directional.
			if (!hasFirst && !hasSecond)
			{
				vp.geneId = first;
				v1 = graphPtr->AddNode(vp);
				vp.geneId = second;
				v2 = graphPtr->AddNode(vp);
			
				geneVertMap[first] = v1;
				geneVertMap[second] = v2;
				
				link.weight = maxScore;
				graphPtr->AddLink(v1, v2, link, link);
				
			}
			else if(hasFirst && !hasSecond)
			{		
				v1 = geneVertMap[first];
				vp.geneId = second;
				v2 = graphPtr->AddNode(vp);
			
				geneVertMap[second] = v2;
			
				link.weight = maxScore;
				graphPtr->AddLink(v1, v2, link, link);
			}
			else if (hasSecond && !hasFirst)
			{
				vp.geneId = first;
				v1 = graphPtr->AddNode(vp);
				v2 = geneVertMap[second];
				
				geneVertMap[first] = v1;
			
				link.weight = maxScore;
				graphPtr->AddLink(v1, v2, link, link);
			}
			else
			{
				v1 = geneVertMap[first];
				v2 = geneVertMap[second];
				
				link.weight = maxScore;
				graphPtr->AddLink(v1, v2, link, link);
			}
		}	
	}
	file.close();

/*		
#if USE_BIDIRECTIONAL
	printf("Loaded %d links between %d nodes in %f seconds.\n", graphPtr->getLinkCount()/2, graphPtr->getNodeCount(), ((clock()-start)+0.0)/CLOCKS_PER_SEC);	
#else 
	printf("Loaded %d links between %d nodes in %f seconds.\n", graphPtr->getLinkCount(), graphPtr->getNodeCount(), ((clock()-start)+0.0)/CLOCKS_PER_SEC);	
#endif
*/
	return (*graphPtr);
}

//TODO: writers		
void BoostGraphIO::writeXGMMLGraph(string path){}

void BoostGraphIO::writeTSVGraph(string path)
{
	ofstream file;
	Graph::Node v1, v2;
	Graph::Link e;
	file.open(path.c_str());
	for (Graph::link_range_t er = graphPtr->getLinks(); er.first != er.second; er.first++)
	{
		e = *er.first;
		graphPtr->getNodesByLink(e, v1, v2);
		file << graphPtr->properties(v1).geneId << "\t" << graphPtr->properties(v2).geneId << "\t" << graphPtr->properties(e).weight << endl;
	}
	file.close();
}





