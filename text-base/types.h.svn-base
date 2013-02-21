#ifndef __TYPES_H__
#define __TYPES_H__

#include "boostgraph.h"

using namespace boost;
using namespace std;


class NodeProperties 
{
	public:
		NodeProperties(){}
		~NodeProperties(){}
	
		string geneId;
		vector<int> connectedDegrees;
};

class LinkProperties 
{
	public:
		LinkProperties(){}
		~LinkProperties(){}
	
		float weight;
};

typedef BoostGraph<NodeProperties, LinkProperties> Graph;

class GeneGroup
{
	public:
		GeneGroup(){}
		~GeneGroup(){}
	
		string groupId;
		string groupSpe;
		string groupSys;
		string groupDesc;
		vector<string> groupGenes;
		
		string inputFilePath;
		
		GeneGroup &operator=(const GeneGroup &rhs)
		{	
			this->groupId = rhs.groupId;
			this->groupSpe = rhs.groupSpe;
			this->groupSys = rhs.groupSys;
			this->groupDesc = rhs.groupDesc;
			this->groupGenes = rhs.groupGenes;
			this->inputFilePath = rhs.inputFilePath;
			
			return *this;
		}
};


class Stats	
{
	public:
		Stats(){this->clear();}
		~Stats(){}
	
		vector<int> linkCount;
		vector<float> clusteringCoeff;
		float expectedLinks;
		float observedLinks;
		float zScore; 
		long double pValue;
		float stdDev; 
		float chiSqr;
	
		void clear()
		{
			linkCount.clear();
			clusteringCoeff.clear();
			expectedLinks = observedLinks = pValue = stdDev = chiSqr = 0.0;
			zScore = 0.0;
		}
		
		Stats &operator=(const Stats &rhs)
		{	
			this->linkCount = rhs.linkCount;
			this->clusteringCoeff = rhs.clusteringCoeff;
			this->expectedLinks = rhs.expectedLinks;
			this->observedLinks = rhs.observedLinks;
			this->zScore = rhs.zScore;
			this->pValue = rhs.pValue;
			this->stdDev = rhs.stdDev;
			this->chiSqr = rhs.chiSqr;

			return *this;
		}
};


struct Record
{
	Graph::Node node;
	int degree;
	//int smetric;
};

#endif
