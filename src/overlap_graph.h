#ifndef __OVERLAP_GRAPH_H__
#define __OVERLAP_GRAPH_H__

#include "directed_graph.h"
#include "hyper_set.h"
#include "config.h"

using namespace std;

typedef pair<bool, int> PBI;

class overlap_graph : public directed_graph
{
public:
	overlap_graph();
	overlap_graph(const hyper_set &hs);
	virtual ~overlap_graph();

public: 
	int contract();
	int extract_paths();
	int clear();

public:
	VVI v2el;            // vertex corresponding paired-path's edge-list
	vector<double> vwrt; // graph vertex weight
	MED ewrt;            // graph edge weight
	MISI ppe2v;          // index: from paired-path's edge to vertices
	

	vector<path> paths;   // extended paired paths

public:
	int set_vertex_edge_list(int v, const vector<int> &el);
	int set_vertex_weight(int v, double w);
	double get_vertex_weight(int v) const;
	int set_edge_weight(edge_descriptor e, double w);
	double get_edge_weight(edge_descriptor e) const;

	int add_vertices(const hyper_set & hs);
	int add_necessary_vertices(const hyper_set & hs);
	int build_index();
	int add_edges();
	PBI compatible(int x, int y);
	
	PEB find_seed_edge(set<edge_descriptor> &used);
	int left_extend(int v);
	int right_extend(int v);
	double get_path_abundance(const vector<int> &v);
	int update(vector<int> &v, double w, set<edge_descriptor> &used);

	int contract_from(int x);

	int write(ostream &out);
};
#endif // !__OVERLAP_GRAPH_H__

