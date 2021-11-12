#ifndef __PHASING_GRAPH_H__
#define __PHASING_GRAPH_H__

using namespace std;

#include "directed_graph.h"
#include "subgraph.h"
#include "path.h"
#include "hyper_set.h"

class phasing_graph : public directed_graph
{
public:
	phasing_graph();
	phasing_graph(splice_graph &gr, VE &i2e, MEI &e2i, vector<path> &pp, hyper_set &hs);
	virtual ~phasing_graph();

public:
	int extract_paths(splice_graph &gr, VE &i2e);
	int write(ostream &out);

public:
	VVI v2el;            // node corresponding edge or edge-list
	vector<double> vwrt; // corresponding splicing graph edge weight or paired-path cov
	MEI ewrt;            // 0 or 1
	int pos;             // split position (point at first paired-path node)

	vector<path> paths;  // paths of edges of splicing graph
public:
	int set_vertex_edge_list(int v, vector<int> &el);
	int set_vertex_weight(int v, double w);
	double get_vertex_weight(int v) const;
	int set_edge_weight(edge_descriptor e, int w);
	int get_edge_weight(edge_descriptor e) const;

	int add_normal_nodes(splice_graph &gr, VE &i2e, MEI &e2i, vector<path> &pp, hyper_set &hs);
	int add_contracted_nodes(vector<path> &pp);

	int find_seed(set<int> &c_used, set<int> &n_used);
	int left_extend(int v);
	int right_extend(int v);
	double get_path_abundance(const vector<int> &v);
	int update(const vector<int> &v, double w, set<int> &c_used, set<int> &n_used, const VE &i2e, int vn);
};
#endif // !__PHASING_GRAPH_H__

