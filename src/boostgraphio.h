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
This file contains BoostGraphIO class declaration that is used to load
and write networks. It allows the user to easily load XGMML, 
simple TSV, and funcoup TSV graphs. It currently only writes 
simple TSV graphs. 

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



