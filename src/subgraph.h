#ifndef __SUBGRAPH_H__
#define __SUBGRAPH_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "gurobi_mip.h"
#include "transcript.h"
#include "config.h"
#include "path.h"
#include "scallop.h"
#include "overlap_graph.h"
#include "phasing_graph.h"
#include <stack>

typedef pair<double, bool> PDB;
typedef pair<int, PDB> PIDB;
typedef map<int, PDB> MIDB;
typedef MIDB::iterator MIDBI;
typedef vector<MIDB> VMIDB;

typedef pair<double, vector<int>> PDVI;
typedef map<double, vector<int>> MDVI;
typedef map<double, vector<int>>::iterator MDVII;

typedef unsigned long UL;
typedef pair<double, UL> PDUL;
typedef vector<PDUL> VDUL;
typedef pair<int, bool> PIB;
typedef map<int, bool> MIB;
typedef vector<MIB> VMIB;
typedef pair<UL, MIB> PULMIB;
typedef map<double, PULMIB> MDULMIB;
typedef vector<MDULMIB> VMDULMIB;
typedef MIB::iterator MIBI;
typedef MDULMIB::iterator MDULMIBI;
typedef set<double> SD;
typedef SD::iterator SDI;
typedef vector<set<double>> VSD;
typedef VSD::iterator VSDI;
typedef pair<double, double> PDD;
typedef map<int, double> MID;
typedef pair<double, MID> PDMID;
typedef vector<PDMID> VDMID;


class subgraph
{
public:
	subgraph();
	subgraph(const splice_graph &g, const hyper_set &h);
	~subgraph();

public:
	int assemble();
	int write(ostream &out, int index);

	// access functions
	int get_edge_index(edge_descriptor e) const;
	edge_descriptor get_edge(int v) const;

public:
	splice_graph gr;
	hyper_set hs;
	MEI e2i;                  // edge map, from edge to index. 
	VE i2e;                   // edge map, from index to edge.

	vector<path> backup;      // backup when lp has no solution
	vector<transcript> trsts; // predicted transcripts

private:
	VVB B;
	MDVI B_covs;
	vector<int> lp_x;
	vector<double> lp_t;

private:
	// find paths
	int specified_find_paths(bool use_dump);
	int dfs(bool use_dump);
	int dfs_heaviest_paths(bool use_dump);
	int dfs_backwards_filter_paths(bool use_dump);
	int hs_extend_paths(bool use_dump);
	int scallop_find_paths(bool save2backup, bool use_dump);
	int heaviest_compatible_paths(bool use_dump);
	int ipac_paired_paths(bool use_dump);

	// find paths use
	int path2line(const path &p, vector<bool> &es);
	int vertex_path(const path & p1, path & p2);
	bool dump_path(const path &p) const;
	int B_add_new_path(const path &p, const vector<bool> &es, int bottom);
	int B_add_path(path & p, int bottom, bool use_dump);
	int stack2path(stack<int> s, path &p);
	int stack2path_forward(stack<int> s, path & p);
	int edges2vertices(vector<int>& e, vector<int>& v);
	int write_stack(ostream &out, stack<int> s);
	int write_VMIDB(ostream &out, const VMIDB &adj);
	int write_VMDULMIB(ostream &out, const VMDULMIB & vps);
	int write_VSD(ostream & out, const VSD & vc);
	// dfs_heaviest_paths use
	VMIDB get_VMIDB_adj();
	int set_path_abd(VMIDB &adj, path &p);
	int extract_path_and_update(VMIDB &adj, path &p);
	int extract_path_and_update(VMIDB &adj, vector<double> &vwrt, path &p);
	int find_heaviest_stackable(VMIDB &adj, int s, vector<bool> &stacked) const;
	int reset_visit(VMIDB &adj, int s);
	int dfs_heaviest_path(VMIDB &adj, path &p);
	// hs_extend_paths use
	int get_heaviest_inner_edge(VMIDB &adj);
	int get_heaviest_predecessor(VMIDB &adj, int ei);
	int get_heaviest_successor(VMIDB &adj, int ei);
	set<int> get_remain_edges(VMIDB &adj);
	int hs_extend_path(VMIDB &adj, hyper_set &hs_copy, path &p);
	int hs_extend_path(VMIDB &adj, path &p);
	int get_heaviest_compatible_predecessor(const MEI &secnts, set<int> &eis, vector<int> &e);
	int get_heaviest_compatible_successor(const MEI & secnts, set<int>& eis, vector<int>& e);
	// heaviest_compatible_paths use
	int heaviest_compatible_path(VMIDB &adj, vector<double> &vwrt, const MEI &secnts, path &p);
	// dfs_backwards_filter_paths use
	VMDULMIB compute_vertices_path_info();
	int dfs_backwards(VMDULMIB & vps, double cov, UL pn, int bottom, bool use_dump);

	//lp
	vector<double> get_edge_covs();
    int compute_hs_span_in_paths(VVB &hss, vector<int> &ucnts, vector<bool> &keep);
	int lp(ostream &out);

	// output trsts
	int single_node_sub_add_trst();
	int Bline2path(const vector<bool>& bi, path & p);
	int reset_lp_xt();
	int multi_nodes_sub_add_trsts(int pn, int lp_score, bool use_backup);

	// write
	int write_all_paths();
	int write_matrix(ostream &out, VVB &M);
	int write_kept_hs(ostream &out, const vector<bool> &keep);
	int write_hss(ostream &out, const VVB &hss, const vector<int> &ucnts);
};

#endif // !__SUBGRAPH_H__


