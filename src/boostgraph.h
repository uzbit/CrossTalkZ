/*

This file contains a Graph container based on the boost::graph library
was taken from the web at:
http://stackoverflow.com/questions/671714/modifying-vertex-properties-in-a-boostgraph
with the exception of a few modifications/additions

*/

#ifndef __BOOST_GRAPH_H__
#define __BOOST_GRAPH_H__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/random.hpp>
//#include <boost/graph/graphml.hpp>
#include <boost/random.hpp>
#include <iostream>
#include <time.h>

#define LINK_NONE			0
#define LINK_12				1
#define LINK_21				2
#define LINK_BOTH			3

#define USE_BIDIRECTIONAL	0
#define USE_VECTOR			0

using namespace boost;
using namespace std;

/*  */

/* definition of basic boost::graph properties */
enum vertex_properties_t { vertex_properties };
enum edge_properties_t { edge_properties };

namespace boost 
{
	BOOST_INSTALL_PROPERTY(vertex, properties);
	BOOST_INSTALL_PROPERTY(edge, properties);
}

/* the graph base class template */
template < typename NODEPROPERTIES, typename LINKPROPERTIES >
class BoostGraph
{
public:

	/* an adjacency_list like we need it */
	typedef adjacency_list<
		setS, // setS disallows parallel edges
#if USE_VECTOR
		vecS, // vertex container
#else
		listS,
#endif 
		undirectedS, // undirected graph
		property<vertex_properties_t, NODEPROPERTIES>,
		property<edge_properties_t, LINKPROPERTIES>
	> GraphContainer;


	/* a bunch of graph-specific typedefs */
	typedef typename graph_traits<GraphContainer>::vertex_descriptor Node;
	typedef typename graph_traits<GraphContainer>::edge_descriptor Link;
	typedef std::pair<Link, Link> LinkPair;

	typedef typename graph_traits<GraphContainer>::vertex_iterator node_iter;
	typedef typename graph_traits<GraphContainer>::edge_iterator link_iter;
	typedef typename graph_traits<GraphContainer>::adjacency_iterator adjacency_iter;
	typedef typename graph_traits<GraphContainer>::out_edge_iterator out_link_iter;

	typedef typename graph_traits<GraphContainer>::degree_size_type degree_t;

	typedef std::pair<adjacency_iter, adjacency_iter> adjacency_node_range_t;
	typedef std::pair<out_link_iter, out_link_iter> out_link_range_t;
	typedef std::pair<node_iter, node_iter> node_range_t;
	typedef std::pair<link_iter, link_iter> link_range_t;

	int id;

	/* constructors etc. */
	BoostGraph(){}

	//Use this constructor to make a new copy of a BoostGraph
	BoostGraph(const BoostGraph& g) :
		graph(g.graph){}

	virtual ~BoostGraph(){}

	/* structure modification methods */
	void Clear()
	{
		graph.clear();
	}

	Node AddNode(const NODEPROPERTIES& prop)
	{
		Node v = add_vertex(graph);
		properties(v) = prop;
		return v;
	}

	void RemoveNode(const Node& v)
	{
		clear_vertex(v, graph);
		remove_vertex(v, graph);
	}

	void RemoveLink(const Node &v1, const Node &v2)
	{
		remove_edge(v1, v2, graph);
		remove_edge(v2, v1, graph);
	}
	
	void RemoveAllLinks()
	{
		Link e;
		for (link_range_t er = getLinks(); er.first != er.second;)
		{
			e = *er.first;
			er.first++;
			remove_edge(e, graph);
		}
	}

	Link AddLink(const Node& v1, const Node& v2, const LINKPROPERTIES& prop_12)
	{
		/* TODO: maybe one wants to check if this edge could be inserted */
		Link addedLink1 = add_edge(v1, v2, graph).first;
		
		properties(addedLink1) = prop_12;
		
		return addedLink1;
	}

	LinkPair AddLink(const Node& v1, const Node& v2, const LINKPROPERTIES& prop_12, const LINKPROPERTIES& prop_21)
	{
#if USE_BIDIRECTIONAL
		/* TODO: maybe one wants to check if this edge could be inserted */
		Link addedLink1 = add_edge(v1, v2, graph).first;
		Link addedLink2 = add_edge(v2, v1, graph).first;

		properties(addedLink1) = prop_12;
		properties(addedLink2) = prop_21;
#else
		Link addedLink1 = add_edge(v1, v2, graph).first;

		properties(addedLink1) = prop_12;

		Link addedLink2 = addedLink1;
#endif
		return LinkPair(addedLink1, addedLink2);
	}


	/* property access */
	NODEPROPERTIES& properties(const Node& v)
	{
		typename property_map<GraphContainer, vertex_properties_t>::type param = get(vertex_properties, graph);
		return param[v];
	}

	const NODEPROPERTIES& properties(const Node& v) const
	{
		typename property_map<GraphContainer, vertex_properties_t>::const_type param = get(vertex_properties, graph);
		return param[v];
	}

	LINKPROPERTIES& properties(const Link& v)
	{
		typename property_map<GraphContainer, edge_properties_t>::type param = get(edge_properties, graph);
		return param[v];
	}

	const LINKPROPERTIES& properties(const Link& v) const
	{
		typename property_map<GraphContainer, edge_properties_t>::const_type param = get(edge_properties, graph);
		return param[v];
	}


	/* selectors and properties */
	const GraphContainer& getGraph() const
	{
		return graph;
	}

	node_range_t getNodes() const
	{
		return vertices(graph);
	}

	link_range_t getLinks() const
	{
		return edges(graph);
	}
	
	Node getRandomNode()
	{
		mt19937 gen(time(0));
		return (Node)random_vertex(graph, gen);
	}
		  
	Link getRandomLink()
	{
		mt19937 gen(time(0));
		return (Link)random_edge(graph, gen);
	}
	
	Link getLinkByIndex(int index)
	{
		link_range_t er = getLinks();
		advance(er.first, index);
		return (Link)(*(er.first));
	}
	
	void getNodesByLink(const Link &e, Node &v1, Node &v2) const 
	{
		v1 = source(e, graph);	
		v2 = target(e, graph);
	}
		  
	//Results in ep containing a pair of edges connecting v1 and v2 vertices
	//Return value is used to determine validity of edges in the pair by:
	//if no edges: LINK_NONE
	//if only an edge points from v1 -> v2: LINK_12
	//if only an edge points from v2 -> v1: LINK_21
	//if both edges are present: LINK_BOTH
	int getLinkPair(const Node &v1, const Node &v2, LinkPair &ep) const
	{
		bool found1, found2;
		tie(ep.first, found1) = edge(v1, v2, graph);
		tie(ep.second, found2) = edge(v2, v1, graph);
		return ((found1?LINK_12:LINK_NONE) + (found2?LINK_21:LINK_NONE));
	}

	bool hasDirectedLink(const Node &v1, const Node &v2)
	{
		return (edge(v1, v2, graph).second);
	}
	
	bool hasLink(const Node &v1, const Node &v2) const
	{
#if USE_BIDIRECTIONAL
		return (edge(v1, v2, graph).second || edge(v2, v1, graph).second);
#else
		return (edge(v1, v2, graph).second);
#endif
	}
	
	adjacency_node_range_t getAdjacentNodes(const Node& v) const
	{
		return adjacent_vertices(v, graph);
	}

	int getNodeCount() const
	{
		return num_vertices(graph);
	}
	
	int getLinkCount() const
	{
		return num_edges(graph);
	}
	
	int getNodeDegree(const Node& v) const
	{
		return out_degree(v, graph);
	}

	/*other functions*/

/*	void writeGraphML(string path)
	{
#if USE_VECTOR
		//TODO: write_graphml()
		typename property_map<GraphContainer, vertex_properties_t>::type paramv = get(vertex_properties, graph);
		typename property_map<GraphContainer, edge_properties_t>::type parame = get(edge_properties, graph);
		dynamic_properties dp; 
		
		dp.property("vertex_property", paramv); 
		dp.property("edge_property", parame); 
		ofstream file;

		file.open(path.c_str());	
		write_graphml(file, graph, dp, false);
		file.close();
#endif
	}*/

	/* operators */
	BoostGraph& operator=(const BoostGraph &rhs)
	{
		graph = rhs.graph;	
		return *this;
	}

protected:
	GraphContainer graph;
};

#endif
