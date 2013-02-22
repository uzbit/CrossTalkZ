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
Preprocessor definitions utilized throughout the CrossTalkZ program.

*/

#ifndef __DEFINES_H__
#define __DEFINES_H__

//The following are definitions for column numbers starting from 0..N
//tab separated file should have format: Protein1\tProtein2[\tScore]
#define FUNCOUP_MAX_SCORE		0
#define FUNCOUP_PROTEIN1		5
#define FUNCOUP_PROTEIN2		6

//A groups file should have format: Protein\tGroup\t[System, Species, Description]
#define GROUP_GENE	0
#define GROUP_ID		1
#define GROUP_SYS		2
#define GROUP_SPE		3
#define GROUP_DESC	4


//if TRIM_GROUPS == 1, remove groups that do not have sufficient
//(< minimumGenesForGroup) representation in network.  
#define TRIM_GROUPS		1 

#define MODE_0			0 //Link isn't counted if either gene belongs to both groups. 
#define MODE_1			1 //Link isn't counted if both genes belong to both groups.

#define METHOD_LINKSWAP		0
#define METHOD_ASSIGN			1
#define METHOD_ASSIGN_SECOND			2
#define METHOD_LABELSWAP		3
#define METHOD_DEFAULT		METHOD_ASSIGN_SECOND

#define VERSION		"1.3.3"

#define DEBUG		0
#if DEBUG							  //debug defines
	#define RAND_GROUPS		0 //Use this to randomize the group members (false positive test).
	#define DO_LINE_LIMIT			0 //0 for no line limit 1 for line limit of MAX_LINES
	#define MAX_LINES		1000000
	#define VERBOSE		0 //prints extra information if 1
	#define SPLIT_GROUPS		0 //Use this to split the groups into random halves. (part of false negative test).
#endif


#endif

