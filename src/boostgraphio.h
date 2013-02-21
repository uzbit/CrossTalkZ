/*

This file contains BoostGraphIO class declaration that is used to load
and write networks. It allows the user to easily load XGMML, 
simple TSV, and funcoup TSV graphs. It currently only writes 
simple TSV graphs. 

Written by Ted McCormack.

*/

#ifndef __BOOSTGRAPHIO_H__
#define __BOOSTGRAPHIO_H__

#include "defines.h"
#include "types.h"

using namespace std;

extern float cutoffScore;
extern bool useCutoff;

class BoostGraphIO
{
public:
	BoostGraphIO(){graphPtr = &graph;}
	
	~BoostGraphIO(){}
	
	Graph& readGraph(string path);
	
	//See: http://en.wikipedia.org/wiki/XGMML
	Graph& readXGMMLGraph(string path);
	
	//readTSVGraph calls readSimpleTSVGraph or readFunCoupTSVGraph
	Graph& readTSVGraph(string path);	 
	
	//readSimpleTSVGraph file with formats:
	//protein1\tprotein2\tscore OR
	//protein1\tprotein2
	Graph& readSimpleTSVGraph(string path); 
	
	//readFunCoupTSVGraph FunCoup file: http://funcoup.sbc.su.se/
	Graph& readFunCoupTSVGraph(string path);
	
	void setGraph(Graph *g){graphPtr = g; graph = *graphPtr;}
	void setGraph(Graph &g){graph = g;}
	Graph& getGraph(){return graph;}
	
	void writeXGMMLGraph(string path);
	void writeTSVGraph(string path);
	
	void clearGraph(){graphPtr->Clear();}
	
protected:
	Graph& __readTSVGraph__(string path, int flag);
	Graph *graphPtr;
	Graph graph;
	
};




#endif



