#include "subgraph.h"

#include <iostream>  
#include <algorithm> 
#include <float.h>  

using namespace std;
int max_remain_paths = 0;
bool PDUL_cmp(PDUL &a, PDUL &b) { return a.first > b.first; }

subgraph::subgraph()
{}

subgraph::subgraph(const splice_graph &g, const hyper_set &h)
	:gr(g), hs(h)
{
	if (output_tex_files) gr.draw(output_prefix + "/" + gr.gid + ".tex");
	gr.get_edge_indices(i2e, e2i);
	//gr.get_inner_edge_indices(i2e, e2i);
	hs.build(gr, e2i);
}

subgraph::~subgraph()
{}

int subgraph::assemble()
{
	ofstream fout;
	if (verbose == -1)
	{
		fout.open(log_file.c_str(), ios_base::app);
		if (fout.fail()) exit(2);
	}

	int vn = gr.num_vertices();
	int en = gr.num_edges();
	int ien = gr.num_inner_edges();
	long pn = gr.compute_num_paths();
	if (verbose >= 1) printf("all path num=%ld\n", pn);
	if (verbose == -1) fout << "all path num=" << pn << "\n";
	if (ien <= 0)
	{
		assert(en == 2);
		assert(vn == 3);
		single_node_sub_add_trst();
		return 0;
	}
	else
	{
		assert(vn >= 4);
		max_remain_paths = max_AB_element / (en);
		if (verbose >= 1) printf("max_remain_paths=%d\n", max_remain_paths);
		if (verbose == -1) fout << "max_remain_paths=" << max_remain_paths << "\n";
		if (pn <= max_remain_paths) B.reserve(pn);
		else B.reserve(max_remain_paths);

		if (output_all_paths) specified_find_paths(true);
		else
		{
			//!output_all_paths
			assert(simple_path_num >= small_path_num);
			if (pn <= simple_path_num) dfs(true);
			else specified_find_paths(true);
		}

		if (output_all_paths)
		{
			write_all_paths();
			return 0;
		}
		B.shrink_to_fit();
		if (B.size() >= 1)
		{
			int lp_score = lp(fout);
			multi_nodes_sub_add_trsts(pn, lp_score, use_backup_paths);
		}
		if (verbose >= 2)
		{
			cout << "B";
			write_matrix(cout, B);
			if (verbose >= 3)
			{
				int in;
				cout << "please input an integer:" << endl;
				cin >> in;
			}
		}
		if (verbose == -1)
		{
			fout << "B";
			write_matrix(fout, B);
		}
	}
	//write_all_paths();//若注释此，注意取消注释assembler.assemble()的write()
	return 0;
}

int subgraph::get_edge_index(edge_descriptor e) const
{
	MEI::const_iterator it = e2i.find(e);
	assert(it != e2i.end());
	return it->second;
}

edge_descriptor subgraph::get_edge(int v) const
{
	assert(v >= 0 && v < i2e.size());
	return i2e[v];
}

// find paths
int subgraph::specified_find_paths(bool use_dump)
{
	for (int i = 0; i < path_options.length(); i++)
	{
		char option = path_options[i];
		//cout << "current option=" << option << endl;
		switch (option)
		{
		case 'a':
			dfs(use_dump);                                 
			break;
		case 'b':
			dfs_heaviest_paths(use_dump);                 
			break;
		case 'c':
			dfs_backwards_filter_paths(use_dump);          
			break;
		case 'd':
			scallop_find_paths(use_backup_paths, use_dump);
			break;
		case 'e':
			heaviest_compatible_paths(use_dump);           
			break;
		case 'f':
			hs_extend_paths(use_dump);                     
			break;
		case 'g':
			//hs_compatible_extend_paths();
			break;
		case 'h':
			ipac_paired_paths(use_dump);
			break;
		case 'i':
			break;
		default:
			break;
		}
	}
	return 0;
}
int subgraph::dfs(bool use_dump)
{
	if (verbose >= 4) cout << "into dfs:" << endl;
	int vn = gr.num_vertices();
	assert(vn > 3);
	int ien = gr.num_inner_edges();
	assert(ien > 0);
	int bottom = (int)B.size() - 1;

	VMIB ava;//adjacent vertices accessed
	ava.reserve(vn);
	vector<bool> stacked;
	stacked.resize(vn, false);
	for (int i = 0; i < vn; i++)
	{
		MIB adj;
		edge_iterator it1, it2;
		PEEI pei;
		for (pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int t = (*it1)->target();
			assert(adj.find(t) == adj.end());
			adj[t] = false;
		}
		ava.push_back(adj);
	}

	//dfs
	stack<int> s1;
	s1.push(0);
	stacked[0] = true;
	while (!s1.empty())
	{
		if (verbose >= 4) write_stack(cout, s1);
		if (s1.top() == vn - 1)
		{
			assert((int)B.size() <= max_remain_paths + 1);
			path p;
			stack2path(s1, p);
			int a = B_add_path(p, bottom, use_dump);
			if (verbose >= 4)
			{
				write_matrix(cout, B);
				int in;
				cout << "please input an integer:" << endl;
				cin >> in;
			}
			if (a == -1) return 0;
		}
		MIB &p = ava[s1.top()];
		MIBI it;
		for (it = p.begin(); it != p.end(); it++)
		{
			if (it->second) continue;
			s1.push(it->first);
			stacked[it->first] = true;//in stack
			it->second = true;        //accessed
			break;
		}

		if (it == p.end())
		{
			for (it = p.begin(); it != p.end(); it++) {
				it->second = false;
			}
			stacked[s1.top()] = false;
			s1.pop();
		}
	}
	return 0;
}
int subgraph::dfs_heaviest_paths(bool use_dump)
{
	if (verbose >= 4) cout << "into dfs_heaviest_paths:" << endl;
	VMIDB adj = get_VMIDB_adj();
	if (verbose >= 4) write_VMIDB(cout, adj);
	int bottom = (int)B.size() - 1;
	int iter = 0;
	while (true)
	{
		if (iter >= 100000) break;
		if (verbose >= 4) cout << "iter-" << iter << ":" << endl;
		path p;
		int a = dfs_heaviest_path(adj, p);
		assert(a == 0 || a == -1);
		if (a == -1) break;
		if (verbose >= 4) cout << "not break" << endl;
		int rt = B_add_path(p, bottom, use_dump);
		if (verbose >= 4)
		{
			write_matrix(cout, B);
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
		if (rt == -1) break;
		iter += 1;
		//system("read");
	}

	return 0;
}
int subgraph::dfs_backwards_filter_paths(bool use_dump)
{
	if (verbose >= 4) cout << "into dfs_backwards_filter_paths:" << endl;
	int vn = gr.num_vertices();
	VMDULMIB vps = compute_vertices_path_info();
	if (verbose >= 4) write_VMDULMIB(cout, vps);
	MDULMIB &tt = vps[vn - 1];
	MDULMIBI it;
	VDUL ttp;
	ttp.reserve(tt.size());
	for (it = tt.begin(); it != tt.end(); ++it)
	{
		ttp.push_back(PDUL(it->first, it->second.first));
	}
	sort(ttp.begin(), ttp.end(), PDUL_cmp);

	int bottom = (int)B.size() - 1;
	for (int i = 0; i < ttp.size(); i++)
	{
		int rt = dfs_backwards(vps, ttp[i].first, ttp[i].second, bottom, use_dump);
		if (rt == -1) break;
	}

	return 0;
}
int subgraph::hs_extend_paths(bool use_dump)
{
	if (verbose >= 4) cout << "into hs_extend_paths:" << endl;
	VMIDB adj = get_VMIDB_adj();
	if (verbose >= 4) write_VMIDB(cout, adj);
	int bottom = (int)B.size() - 1;
	int iter = 0;
	while (true)
	{
		if (verbose >= 4) cout << "iter-" << iter << ":" << endl;
		path p;
		int rtn = hs_extend_path(adj, p);
		assert(rtn == 0 || rtn == -1);
		if (rtn == -1) break;
		rtn = B_add_path(p, bottom, use_dump);
		if (verbose >= 4)
		{
			write_matrix(cout, B);
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
		if (rtn == -1) break;
		iter += 1;
	}
	return 0;
}
int subgraph::scallop_find_paths(bool save2backup, bool use_dump)
{
	if (verbose >= 4) cout << "into scallop_find_paths:" << endl;
	int bottom = (int)B.size() - 1;
	scallop sc(gr, hs);
	sc.assemble();
	vector<path> &paths = sc.paths;
	if (save2backup) backup = paths;
	for (int i = 0; i < paths.size(); i++)
	{
		int rtn = B_add_path(paths[i], bottom, use_dump);
		if (rtn == -1) break;
	}
	return 0;
}
int subgraph::heaviest_compatible_paths(bool use_dump)
{
	if (verbose >= 4) cout << "into heaviest_compatible_paths:" << endl;
	VMIDB adj = get_VMIDB_adj();
	vector<double> vwrt = gr.vwrt;
	MEI secnts = hs.get_single_edge_counts(gr, i2e);
	if (verbose >= 4)
	{
		cout << "secnts: ";
		MEI::iterator it;
		for (it = secnts.begin(); it != secnts.end(); ++it) { cout << get_edge_index(it->first) << "=" << it->second << ", "; }
		cout << endl;
	}
	int bottom = (int)B.size() - 1;
	int iter = 0;
	while (true)
	{
		if (verbose >= 4) cout << "iter-" << iter << ":" << endl;
		path p;
		int rtn = heaviest_compatible_path(adj, vwrt, secnts, p);
		assert(rtn == 0 || rtn == -1);
		if (rtn == -1) break;
		rtn = B_add_path(p, bottom, use_dump);
		if (verbose >= 4)
		{
			write_matrix(cout, B);
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
		if (rtn == -1) break;
		iter += 1;
	}
	return 0;
}
int subgraph::ipac_paired_paths(bool use_dump)
{
	if (verbose >= 5) cout << "into ipac_paired_paths:" << endl;
	int bottom = (int)B.size() - 1;
	overlap_graph og(hs);
	og.extract_paths();
	if (verbose >= 6)
	{
		int in;
		cout << "please input an integer:" << endl;
		cin >> in;
	}
	phasing_graph pg(gr, i2e, e2i, og.paths, hs);
	pg.extract_paths(gr, i2e);
	vector<path> &paths = pg.paths;
	for (int i = 0; i < paths.size(); i++)
	{
		path p;
		vertex_path(paths[i], p);
		int rtn = B_add_path(p, bottom, use_dump);
		if (rtn == -1) break;
	}
	return 0;
}

// add paths use
int subgraph::path2line(const path &p, vector<bool> &es)
{
	int vn = gr.num_vertices();
	int en = gr.num_edges();
	int ien = gr.num_inner_edges();
	assert(vn >= 4);
	assert(ien > 0);

	const vector<int> &v = p.v;
	assert(v.size() >= 4);
	assert(v[0] == 0);
	assert(v.back() == vn - 1);
	es.resize(en, false);
	for (int i = 0; i < (int)v.size() - 1; i++)
	{
		PEB e = gr.edge(v[i], v[i + 1]);
		assert(e.second);
		int ei = get_edge_index(e.first);
		es[ei] = true;
	}
	return 0;
}
int subgraph::vertex_path(const path &p1, path &p2)
{
	const vector<int> &v = p1.v;
	set<int> s;
	edge_descriptor e = get_edge(v[0]);
	int ss = e->source();
	assert(ss == 0);
	s.insert(ss);
	for (int i = 0; i < v.size(); i++)
	{
		edge_descriptor e = get_edge(v[i]);
		s.insert(e->target());
	}
	p2.v.assign(s.begin(), s.end());
	p2.abd = p1.abd;
	return 0;
}
bool subgraph::dump_path(const path &p) const 
{
	if (verbose >= 4) cout << "dump_path: p.abd=" << p.abd << ", p.length=" << p.length << ", p.exon_num=" << gr.num_exons(p) << endl;
	bool dump = false;
	int en = gr.num_exons(p);
	int minl = min_transcript_length_base + en*min_transcript_length_increase;
	for (int i = 0; i < filter_path.length(); i++)
	{
		switch (filter_path[i])
		{
		case 'a'://cov<1.01
			if (en >= 2 && p.abd < min_transcript_coverage) dump = true;
			break;
		case 'b':
			if (en == 1 && p.abd < min_single_exon_coverage) dump = true;
			break;
		case 'c'://length<150+*50
			if (p.length < minl) dump = true;
			break;
		case 'd'://length<200
			if (p.length < 200) dump = true;
			break;
		default:
			break;
		}
	}
	return dump;
}
int subgraph::B_add_new_path(const path &p, const vector<bool> &es, int bottom) 
{
	MDVII it = B_covs.find(p.abd);
	if (it == B_covs.end())
	{
		B.push_back(es);
		int index = B.size() - 1;
		B_covs.insert(PDVI(p.abd, { index }));
		return 0;
	}
	else
	{
		vector<int> &v = it->second;
		bool repeat = false;
		for (int i = 0; i < v.size(); i++)
		{
			if (v[i] > bottom) break;
			if (es == B[v[i]])
			{
				repeat = true;
				break;
			}
		}
		if (!repeat)
		{
			B.push_back(es);
			int index = B.size() - 1;
			it->second.push_back(index);
			return 0;
		}
	}
	return -1;
}
int subgraph::B_add_path(path &p, int bottom, bool use_dump) 
{
	if (B.size() >= max_remain_paths) return -1;
	gr.set_path_abd(p);
	gr.set_path_length(p);
	if (use_dump)
	{
		bool dump = dump_path(p);
		if (dump&&verbose >= 4) cout << "dump this path." << endl;
		if (dump) return 0;
	}
	if (verbose >= 4) cout << "didn't dump this path." << endl;
	vector<bool> es;
	path2line(p, es);
	int a = B_add_new_path(p, es, bottom);
	if (verbose >= 4 && a == 0) cout << "add this path to AB." << endl;
	if (verbose >= 4 && a == -1) cout << "didn't add this path to AB." << endl;
	//PE-support?
	
	return 0;
}
int subgraph::stack2path(stack<int> s, path &p)
{
	int n = s.size();
	vector<int> &v = p.v;
	assert(v.empty());
	v.resize(n);
	assert(v.size() == v.capacity());
	for (int i = n - 1; i >= 0; i--)
	{
		v[i] = s.top();
		s.pop();
	}
	return 0;
}
int subgraph::stack2path_forward(stack<int> s, path &p)
{
	int n = s.size();
	vector<int> &v = p.v;
	assert(v.empty());
	v.resize(n);
	for (int i = 0; i < n; i++)
	{
		v[i] = s.top();
		s.pop();
	}
	return 0;
}
int subgraph::edges2vertices(vector<int> &e, vector<int> &v)
{
	v.clear();
	v.push_back(i2e[e[0]]->source());
	for (int i = 0; i < e.size(); i++) { v.push_back(i2e[e[i]]->target()); }
	return 0;
}
int subgraph::write_stack(ostream &out, stack<int> s)
{
	while (!s.empty())
	{
		out << s.top() << ",";
		s.pop();
	}
	out << endl;
	return 0;
}
int subgraph::write_VMIDB(ostream &out, const VMIDB &adj)
{
	for (int i = 0; i < adj.size(); i++)
	{
		const MIDB &vs = adj[i];
		MIDB::const_iterator it;
		out << i << "={";
		for (it = vs.begin(); it != vs.end(); ++it)
		{
			out << "<" << it->first << ", " << it->second.first << ", " << it->second.second << "> ";
		}
		out << "}\n";
	}
	return 0;
}
int subgraph::write_VMDULMIB(ostream &out, const VMDULMIB &vps)
{
	for (int i = 0; i < vps.size(); i++)
	{
		const MDULMIB &v = vps[i];
		MDULMIB::const_iterator it1;
		MIB::const_iterator it2;
		out << i << "={";
		for (it1 = v.begin(); it1 != v.end(); ++it1)
		{
			out << "<" << it1->first << "," << it1->second.first << ",(";
			const MIB &vi = it1->second.second;
			for (it2 = vi.begin(); it2 != vi.end(); ++it2)
			{
				out << it2->first << "=" << it2->second << ",";
			}
			out << ")> ";
		}
		out << "}\n";
	}
	return 0;
}
int subgraph::write_VSD(ostream &out, const VSD &vc)
{
	SDI it;
	for (int i = 0; i < vc.size(); i++)
	{
		const SD &v = vc[i];
		out << i << "={";
		for (it = v.begin(); it != v.end(); ++it)
		{
			out << (*it) << ",";
		}
		out << "}, ";
	}
	out << "\n";
	return 0;
}
// dfs_heaviest_paths use
VMIDB subgraph::get_VMIDB_adj()
{
	int vn = gr.num_vertices();
	VMIDB adj;
	adj.reserve(vn);
	for (int i = 0; i < vn; i++)
	{
		MIDB vadj;
		edge_iterator it1, it2;
		PEEI pei;
		for (pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; ++it1)
		{
			int t = (*it1)->target();
			double ew = gr.get_edge_weight(*it1);
			assert(vadj.find(t) == vadj.end());
			vadj[t] = PDB(ew, false);
		}
		adj.push_back(vadj);
	}
	return adj;
}
int subgraph::set_path_abd(VMIDB &adj, path &p)
{
	double min = 0;
	bool first = true;
	vector<int> &v = p.v;
	int pl = v.size();
	assert(pl >= 2);
	for (int i = 1; i < pl; i++)
	{
		int s = v[i - 1];
		int t = v[i];
		double cov = adj[s][t].first;
		if (first)
		{
			min = cov;
			first = false;
		}
		else if (cov < min) min = cov;
	}
	assert(min > 0);
	p.abd = min;
	return 0;
}
int subgraph::extract_path_and_update(VMIDB &adj, path &p)
{
	// 1. get path cov
	set_path_abd(adj, p);
	vector<int> &v = p.v;

	// 2. extract path
	for (int i = 1; i < v.size(); i++)
	{
		int s = v[i - 1];
		int t = v[i];
		PDB &st = adj[s][t];
		st.first -= p.abd;
		if (st.first <= 0) adj[s].erase(t);
	}

	// 3. backwards find vertices without out-edge(except tt), and delete in-edges to them.
	int vn = gr.num_vertices();
	for (int i = vn - 2; i >= 0; i--)
	{
		if (!adj[i].empty()) continue;
		for (int j = 0; j < i; j++)
		{
			MIDB &s = adj[j];
			if (s.find(i) == s.end()) continue;
			s.erase(i);
		}
	}

	// 4. forwards find vertices without in-edge(except ss), and delete out-edges from them.
	for (int i = 1; i < vn; i++)
	{
		bool has_in_edge = false;
		int j;
		for (j = 0; j < i; j++)
		{
			MIDB &s = adj[j];
			if (s.find(i) != s.end()) break;
		}
		if (j == i) adj[i].clear();
	}
	// 3,4后，adj应该所有点都既有入边，又有出边（除sstt）,或孤岛
	for (int i = 0; i < vn; i++)
	{
		if (i == 0 || i == vn - 1) continue;
		bool in, out;
		(adj[i].empty()) ? out = false : out = true;
		int j;
		for (j = 0; j < i; j++)
		{
			MIDB &s = adj[j];
			if (s.find(i) != s.end()) break;
		}
		(j < i) ? in = true : in = false;
		assert((in&&out) || (!in && !out));
	}

	// 5. reset visit-tag
	for (int i = 0; i < vn; i++)
	{
		MIDB &vs = adj[i];
		MIDBI it;
		for (it = vs.begin(); it != vs.end(); ++it) it->second.second = false;
	}

	if (verbose >= 4)
	{
		cout << "extract path: ";
		p.print(0);
		write_VMIDB(cout, adj);
	}
	return 0;
}
int subgraph::extract_path_and_update(VMIDB &adj, vector<double> &vwrt, path &p)
{
	// 必须跟在上一函数之后
	vector<int> &v = p.v;
	for (int i = 0; i < v.size(); i++) { vwrt[v[i]] -= p.abd; }
	for (int i = 0; i < adj.size(); i++) { if (adj[i].empty()) vwrt[i] = 0; }
	return 0;
}
int subgraph::find_heaviest_stackable(VMIDB &adj, int s, vector<bool> &stacked) const
{
	int t = -1;
	double max = 0;
	MIDB &svs = adj[s];
	MIDBI it;
	for (it = svs.begin(); it != svs.end(); ++it)
	{
		if (stacked[it->first]) continue;
		if (it->second.second) continue;
		if (it->second.first > max)
		{
			max = it->second.first;
			t = it->first;
		}
	}
	return t;
}
int subgraph::reset_visit(VMIDB &adj, int s)
{
	MIDB &svs = adj[s];
	MIDBI it;
	for (it = svs.begin(); it != svs.end(); ++it)
	{
		it->second.second = false;
	}
	return 0;
}
int subgraph::dfs_heaviest_path(VMIDB &adj, path &p)
{
	int vn = gr.num_vertices();
	assert(p.v.empty());
	vector<bool> stacked(vn, false);

	stack<int> s1;
	s1.push(0);
	stacked[0] = true;
	while (!s1.empty())
	{
		if (verbose >= 4) write_stack(cout, s1);
		if (verbose >= 4) write_VMIDB(cout, adj);
		int s = s1.top();
		if (s == vn - 1)
		{
			stack2path(s1, p);
			extract_path_and_update(adj, p);
			return 0;
		}

		int t = find_heaviest_stackable(adj, s, stacked);
		if (t != -1)
		{
			adj[s][t].second = true;
			s1.push(t);
			stacked[t] = true;
		}
		else
		{
			stacked[s] = false;
			reset_visit(adj, s);
			s1.pop();
		}
	}
	return -1;
}
// hs_extend_paths use
int subgraph::get_heaviest_inner_edge(VMIDB &adj)
{
	int ei = -1;
	double max = 0;
	int vn = gr.num_vertices();
	for (int i = 1; i <= vn - 3; i++)
	{
		MIDB &vs = adj[i];
		MIDBI it;
		for (it = vs.begin(); it != vs.end(); ++it)
		{
			double w = it->second.first;
			if (w <= max) continue;
			max = w;
			PEB e = gr.edge(i, it->first);
			assert(e.second);
			ei = e2i[e.first];
		}
	}
	return ei;
}
int subgraph::get_heaviest_predecessor(VMIDB &adj, int ei)
{
	int lei = -1;
	double max = 0;
	int s = i2e[ei]->source();
	for (int i = 0; i < s; i++)
	{
		MIDB &vs = adj[i];
		MIDBI it = vs.find(s);
		if (it == vs.end()) continue;
		if (it->second.first <= max) continue;
		max = it->second.first;
		PEB le = gr.edge(i, s);
		assert(le.second);
		lei = get_edge_index(le.first);
	}
	return lei;
}
int subgraph::get_heaviest_successor(VMIDB &adj, int ei)
{
	int vn = gr.num_vertices();
	int rei = -1;
	double max = 0;
	int t = i2e[ei]->target();
	if (t == vn - 1) return -1;
	MIDB &vs = adj[t];
	MIDBI it;
	for (it = vs.begin(); it != vs.end(); ++it)
	{
		if (it->second.first <= max) continue;
		max = it->second.first;
		PEB re = gr.edge(t, it->first);
		assert(re.second);
		rei = get_edge_index(re.first);
	}
	return rei;
}
set<int> subgraph::get_remain_edges(VMIDB &adj)
{
	set<int> eis;
	for (int i = 0; i < adj.size(); i++)
	{
		if (adj[i].empty()) continue;
		MIDB &vs = adj[i];
		MIDBI it;
		for (it = vs.begin(); it != vs.end(); ++it)
		{
			PEB e = gr.edge(i, it->first);
			assert(e.second);
			int ei = get_edge_index(e.first);
			eis.insert(ei);
		}
	}
	return eis;
}
int subgraph::hs_extend_path(VMIDB &adj, hyper_set &hs_copy, path &p)
{
	//extract_path_and_update(hs_copy, p);//怎么更新hs呢？
	return 0;
}
int subgraph::hs_extend_path(VMIDB &adj, path &p)
{
	int vn = gr.num_vertices();
	vector<int> &v = p.v;
	assert(v.empty());
	int seed = get_heaviest_inner_edge(adj);
	if (seed == -1) return -1;
	int s = i2e[seed]->source();
	int t = i2e[seed]->target();
	v.push_back(s);
	v.push_back(t);
	set<int> eis = get_remain_edges(adj);
	if (verbose >= 4)
	{
		cout << "current p.v: ";
		writev(cout, v);
		cout << "\nremain edges: ";
		writes(cout, eis);
		cout << endl;
	}
	// left extend
	int ei = seed;
	while (true)
	{
		int lei = hs.get_heaviest_predecessor(eis, ei);
		if (lei == -1) lei = get_heaviest_predecessor(adj, ei);
		if (lei == -1) break;
		assert(i2e[lei]->target() == i2e[ei]->source());
		int les = i2e[lei]->source();
		v.insert(v.begin(), les);
		ei = lei;
		if (verbose >= 4)
		{
			cout << "left extend: ";
			writev(cout, v);
			cout << endl;
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	// right extend
	ei = seed;
	while (true)
	{
		int rei = hs.get_heaviest_successor(eis, ei);
		if (rei == -1) rei = get_heaviest_successor(adj, ei);
		if (rei == -1) break;
		assert(i2e[rei]->source() == i2e[ei]->target());
		int ret = i2e[rei]->target();
		v.push_back(ret);
		ei = rei;
		if (verbose >= 4)
		{
			cout << "right extend: ";
			writev(cout, v);
			cout << endl;
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	assert(v[0] == 0 && v.back() == vn - 1);
	extract_path_and_update(adj, p);
	return 0;
}
// heaviest_compatible_paths use
int subgraph::get_heaviest_compatible_predecessor(const MEI &secnts, set<int> &eis, vector<int> &e)
{
	MI m = hs.get_compatible_predecessors(e);
	MEI::const_iterator it;
	for (it = secnts.begin(); it != secnts.end(); ++it)
	{
		if (it->first->target() != i2e[e[0]]->source()) continue;
		int ei = get_edge_index(it->first);
		if (m.find(ei) == m.end()) m.insert(PI(ei, it->second));
		else m[ei] += it->second;
	}
	return hs.get_heaviest_exist_compatible(m, eis, e);
}
int subgraph::get_heaviest_compatible_successor(const MEI &secnts, set<int> &eis, vector<int> &e)
{
	MI m = hs.get_compatible_successors(e);
	MEI::const_iterator it;
	for (it = secnts.begin(); it != secnts.end(); ++it)
	{
		if (it->first->source() != i2e[e.back()]->target()) continue;
		int ei = get_edge_index(it->first);
		if (m.find(ei) == m.end()) m.insert(PI(ei, it->second));
		else m[ei] += it->second;
	}
	return hs.get_heaviest_exist_compatible(m, eis, e);
}
int subgraph::heaviest_compatible_path(VMIDB &adj, vector<double> &vwrt, const MEI &secnts, path &p)
{
	vector<int> ep;//edge path
	vector<int> &v = p.v;
	// find heaviest vertex
	vector<double>::iterator it = max_element(vwrt.begin(), vwrt.end());
	if (*it <= 0) return -1;
	int vi = distance(vwrt.begin(), it);
	assert(vi >= 0 && vi < vwrt.size());

	// find first heaviest compatible edge
	int lvi = -1;
	double max = 0;
	for (int i = 0; i < vi; i++)
	{
		MIDB &vadj = adj[i];
		if (vadj.find(vi) == vadj.end()) continue;
		if (vadj[vi].first <= max) continue;
		max = vadj[vi].first;
		lvi = i;
	}
	assert(lvi != -1);
	PEB e = gr.edge(lvi, vi);
	assert(e.second);
	int seed = get_edge_index(e.first);
	ep.push_back(seed);
	set<int> eis = get_remain_edges(adj);
	if (verbose >= 4)
	{
		cout << "vwrt: ";
		writev(cout, vwrt);
		cout << "\nadj: ";
		write_VMIDB(cout, adj);
		cout << "ep: ";
		writev(cout, ep);
		cout << "\nremain edges: ";
		writes(cout, eis);
		cout << endl;
	}
	// left extend
	int ei = seed;
	while (true)
	{
		int lei = get_heaviest_compatible_predecessor(secnts, eis, ep);
		if (lei == -1) lei = get_heaviest_predecessor(adj, ei);
		if (lei == -1) break;
		assert(i2e[lei]->target() == i2e[ei]->source());
		ep.insert(ep.begin(), lei);
		ei = lei;
		if (verbose >= 4)
		{
			cout << "left extend(edges): ";
			writev(cout, ep);
			cout << endl;
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	//right extend
	ei = seed;
	while (true)
	{
		int rei = get_heaviest_compatible_successor(secnts, eis, ep);
		if (rei == -1) rei = get_heaviest_successor(adj, ei);
		if (rei == -1) break;
		assert(i2e[rei]->source() == i2e[ei]->target());
		ep.push_back(rei);
		ei = rei;
		if (verbose >= 4)
		{
			cout << "right extend(edges): ";
			writev(cout, ep);
			cout << endl;
			int in = 0;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	assert(i2e[ep[0]]->source() == 0);
	assert(i2e[ep.back()]->target() == (int)adj.size() - 1);
	edges2vertices(ep, v);
	extract_path_and_update(adj, p);
	extract_path_and_update(adj, vwrt, p);
	return 0;
}
// dfs_backwards_filter_paths use
VMDULMIB subgraph::compute_vertices_path_info()
{
	int vn = gr.num_vertices();
	VMDULMIB vps(vn);
	MIB pre;
	pre[-1] = false;
	vps[0][DBL_MAX] = PULMIB(1, pre);
	for (int i = 1; i < vn; i++)
	{
		MDULMIBI it;
		MDULMIB &t = vps[i];
		PEEI pei = gr.in_edges(i);
		edge_iterator it1, it2;
		for (it1 = pei.first, it2 = pei.second; it1 != it2; ++it1)
		{
			int si = (*it1)->source();
			int ti = (*it1)->target();
			assert(ti == i);
			double cov = gr.get_edge_weight((*it1));
			MDULMIB &s = vps[si];
			for (it = s.begin(); it != s.end(); ++it)
			{
				if (cov > it->first)
				{
					if (t.find(it->first) == t.end())
					{
						MIB tpre;
						tpre[si] = false;
						t[it->first] = PULMIB(it->second.first, tpre);
					}
					else
					{
						PULMIB &s_info = it->second;
						PULMIB &t_info = t[it->first];
						t_info.first += s_info.first;
						t_info.second.insert(PIB(si, false));
					}
				}
				else
				{
					if (t.find(cov) == t.end())
					{
						MIB tpre;
						tpre[si] = false;
						t[cov] = PULMIB(it->second.first, tpre);
					}
					else
					{
						PULMIB &s_info = it->second;
						PULMIB &t_info = t[cov];
						t_info.first += s_info.first;
						t_info.second.insert(PIB(si, false));
					}
				}
			}
		}
	}

	return vps;
}
int subgraph::dfs_backwards(VMDULMIB &vps, double cov, UL pn, int bottom, bool use_dump)
{
	if (verbose >= 4) cout << "find paths with cov " << cov << ":" << endl;
	int vn = gr.num_vertices();
	UL cur_pn = 0;
	VSD vc(vn);//used cov for each vertex

	stack<int> s1;
	stack<PDD> pre_covs;//t之前连到点的<边的厚度,实际取得厚度>
	s1.push(vn - 1);
	pre_covs.push(PDD(DBL_MAX, cov));
	while (!s1.empty())
	{
		if (verbose >= 4)
		{
			write_stack(cout, s1);
			write_VMDULMIB(cout, vps);
			write_VSD(cout, vc);
		}
		int t_idx = s1.top();
		int s_idx = -1;
		double s_cov = 0;
		double pre_e_cov = pre_covs.top().first;
		double pre_t_cov = pre_covs.top().second;
		MDULMIB &t = vps[t_idx];
		SD &vct = vc[t_idx];
		MDULMIBI it;
		MIBI it1;
		if (t_idx == 0)
		{
			cur_pn += 1;
			path p;
			stack2path_forward(s1, p);
			int a = B_add_path(p, bottom, use_dump);
			if (verbose >= 4)
			{
				write_matrix(cout, B);
				int in = 0;
				cout << "please input an integer:" << endl;
				cin >> in;
			}
			if (a == -1) return -1;
		}

		if (pre_e_cov == pre_t_cov)
		{
			assert(t_idx != vn - 1);
			for (it = t.begin(); it != t.end(); ++it)
			{
				s_cov = it->first;
				if (s_cov < pre_t_cov) continue;
				if (vct.find(s_cov) != vct.end()) continue;
				MIB &tss = it->second.second;
				for (it1 = tss.begin(); it1 != tss.end(); ++it1)
				{
					if (it1->second) continue;
					s_idx = it1->first;
					it1->second = true;
					if (++it1 == tss.end())
					{
						vct.insert(s_cov);
						assert(vct.size() > 0);
						--it1;
					}
					break;
				}
				if (it1 != tss.end()) break;
			}
		}
		else
		{
			assert(pre_e_cov > pre_t_cov);
			assert(t.find(pre_t_cov) != t.end());
			if (vct.find(pre_t_cov) == vct.end())
			{
				MIB &tss = t[pre_t_cov].second;
				for (it1 = tss.begin(); it1 != tss.end(); ++it1)
				{
					if (it1->second) continue;
					s_idx = it1->first;
					s_cov = pre_t_cov;
					it1->second = true;
					if (++it1 == tss.end())
					{
						vct.insert(s_cov);
						assert(vct.size() > 0);
						--it1;
					}
					break;
				}
			}
		}

		if (s_idx != -1)
		{
			assert(s_cov > 0);
			s1.push(s_idx);
			PEB e = gr.edge(s_idx, t_idx);
			assert(e.second);
			double st = gr.get_edge_weight(e.first);
			pre_covs.push(PDD(st, s_cov));
		}
		else
		{
			SDI it2;
			assert(vct.size() > 0);
			for (it2 = vct.begin(); it2 != vct.end(); ++it2)
			{
				assert(t.find(*it2) != t.end());
				MIB &tss = t[*it2].second;
				for (it1 = tss.begin(); it1 != tss.end(); ++it1)
				{
					it1->second = false;
				}
			}
			s1.pop();
			pre_covs.pop();
			vct.clear();
		}
	}

	assert(pn == cur_pn);
	return 0;
}

// lp
vector<double> subgraph::get_edge_covs()
{
	vector<double> i2w;
	i2w.reserve(i2e.size());
	for (int i = 0; i < i2e.size(); i++)
	{
		assert(gr.ewrt.find(i2e[i]) != gr.ewrt.end());
		double ew = gr.ewrt.find(i2e[i])->second;
		i2w.push_back(ew);
	}
	return i2w;
}
int subgraph::compute_hs_span_in_paths(VVB &hss, vector<int> &ucnts, vector<bool> &keep)
{
	VVI &edges = hs.edges;
	vector<int> &eucnts = hs.eucnts;
	int pn = B.size();
	int en = edges.size();
	assert(hss.empty());
	hss.reserve(en);
	ucnts.reserve(en);
	keep.resize(en, false);
	for (int i = 0; i < edges.size(); i++)
	{
		for (int j = 0; j < edges[i].size(); j++) assert(edges[i][j] >= 0);
		if (eucnts[i] <= 0) continue;
		vector<bool> els;
		els.reserve(pn);
		vector<int> &el = edges[i];
		int el_sum = 0;
		for (int j = 0; j < pn; j++)
		{
			bool intersection = true;
			vector<bool> &bl = B[j];
			for (int k = 0; k < el.size(); k++)
			{
				intersection = (intersection&&bl[el[k]]);
			}
			if (intersection) el_sum += 1;
			els.push_back(intersection);
		}
		if (el_sum <= 0) continue;
		hss.push_back(els);
		ucnts.push_back(eucnts[i]);
		keep[i] = true;
	}
	hss.shrink_to_fit();
	ucnts.shrink_to_fit();
	assert(hss.size() == ucnts.size());
	return 0;
}
int subgraph::lp(ostream &out)
{
	VVB hss;
	vector<int> ucnts;
	vector<bool> keep;
	compute_hs_span_in_paths(hss, ucnts, keep);
	if (verbose >= 3)
	{
		cout << "kept hyper-edges:\n";
		write_kept_hs(cout, keep);
		write_hss(cout, hss, ucnts);
	}
	if (verbose == -1)
	{
		out << "kept hyper-edges:\n";
		write_kept_hs(out, keep);
		write_hss(out, hss, ucnts);
	}
	vector<double> i2w = get_edge_covs();
	gurobi grb(B, i2w, hss, ucnts);
	grb.gid = gr.gid;
	grb.lp();
	lp_t = grb.ts;
	lp_x = grb.xs;
	return grb.score;
}

// output trsts
int subgraph::single_node_sub_add_trst()
{
	double abd = gr.get_vertex_weight(1);
	//if (abd < min_single_exon_coverage) return 0;
	if (abd < min_single_node_coverage) return 0;
	vertex_info vi = gr.get_vertex_info(1);
	int vl = vi.rpos - vi.lpos;
	int minl = min_transcript_length_base + min_transcript_length_increase;
	if (vl < minl) return 0;

	vector<path> paths;
	path p;
	p.v = { 0,1,2 };
	p.abd = abd;
	paths.push_back(p);
	trsts.clear();
	gr.output_transcripts(trsts, paths);
	return 0;
}
int subgraph::Bline2path(const vector<bool> &bi, path &p)
{
	assert(p.v.empty());
	set<int> s;
	for (int j = 0; j < bi.size(); j++)
	{
		if (!bi[j]) continue;
		edge_descriptor e = i2e[j];
		s.insert(e->source());
		s.insert(e->target());
	}
	vector<int> v(s.begin(), s.end());
	sort(v.begin(), v.end());
	int vn = gr.num_vertices();
	assert(v[0] == 0 && v.back() == vn - 1);
	p.v = v;
	return 0;
}
int subgraph::reset_lp_xt()
{
	if (B.size() <= 0) return 0;
	lp_x.clear();
	lp_t.clear();
	lp_x.resize(B.size(), 1);
	lp_t.resize(B.size(), 1.01);
	/*for (int i = 0; i < B.size(); i++)
	{
		path p;
		Bline2path(B[i], p);
		gr.set_path_abd(p);
		lp_t[i] = (p.abd > 1.01) ? p.abd : 1.01;
	}*/
	return 0;
}
int subgraph::multi_nodes_sub_add_trsts(int pn, int lp_score, bool use_backup)
{
	if (pn > small_path_num && lp_score == -1)
	{
		if (use_backup)
		{
			trsts.clear();
			gr.output_transcripts(trsts, backup);
		}
		return 0;
	}

	switch (small_path_mode)
	{
	case 1:
		if (verbose >= 5) cout << "into small_path_mode==1" << endl;
		if (pn <= small_path_num) reset_lp_xt();
		break;
	case 2:
		if (verbose >= 5) cout << "into small_path_mode==2" << endl;
		if (pn <= small_path_num && lp_score <= 0)
		{
			B.clear();
			dfs(true);
			reset_lp_xt();
		}
		break;
	case 3:
		if (verbose >= 5) cout << "into small_path_mode==3" << endl;
		if (pn <= small_path_num && lp_score <= 0) reset_lp_xt();
		break;
	default:
		break;
	}
	
	vector<path> paths;
	for (int i = 0; i < B.size(); i++)
	{
		if (lp_x[i] == 0) continue;
		path p;
		Bline2path(B[i], p);
		p.abd = lp_t[i];
		paths.push_back(p);
	}
	trsts.clear();
	gr.output_transcripts(trsts, paths);//path必须[0,num_vertices()-1]
	return 0;
}

// write
int subgraph::write(ostream &out, int index)
{
	int vn = gr.num_vertices();
	int en = gr.num_edges();
	int ien = gr.num_inner_edges();
	out << "subgraph-" << index << ": " << vn << " vertices, " << en << " edges, " << ien << " inner-edges\n";
	if (verbose <= 1 && verbose != -1) return 0;
	PEEI pei = gr.edges();
	edge_iterator it1 = pei.first, it2 = pei.second;
	for (; it1 != it2; ++it1)
	{
		int ei = get_edge_index(*it1);
		double ew = gr.get_edge_weight(*it1);
		out << "edge " << ei << ":" << (*it1)->source() << "->" << (*it1)->target() << ":" << ew << "\n";
	}
	hs.write_edges(cout);
	return 0;
}

int subgraph::write_all_paths()
{
	//cout << "into write_all_paths---------------------------------" << endl;
	ofstream fout(output_file.c_str(), ios_base::app);
	if (fout.fail()) exit(1);

	vector<transcript> ts;
	vector<path> paths;
	for (int i = 0; i < B.size(); i++)
	{
		path p;
		Bline2path(B[i], p);
		p.abd = 1;
		paths.push_back(p);
	}
	gr.output_transcripts(ts, paths);
	for (int i = 0; i < ts.size(); i++)
	{
		transcript &t = ts[i];
		t.assign_RPKM(1);
		t.write(fout);
		//cout << "write 1 transcript----------------------------------" << endl;
	}
	fout.close();
	return 0;
}

int subgraph::write_matrix(ostream &out, VVB &M)
{
	if (M.size() == 0)
	{
		out << "(0*0):\n";
		return 0;
	}

	out << "(" << M.size() << "*" << M[0].size() << "): " << endl;
	out << "\t";
	for (int i = 0; i < M[0].size(); i++) {
		if (i % 10 == 0) {
			out << i / 10 << "          ";//10个空格
		}
	}
	out << '\n';
	for (int i = 0; i < M.size(); i++) {
		if (i % 10 == 0) {
			out << i / 10 << '\t';
		}
		else {
			out << '\t';
		}
		for (int j = 0; j < M[0].size(); j++) {
			if (j % 10 == 0) {
				out << ' ';
			}
			out << M[i][j];
		}
		out << '\n';
	}
	return 0;
}

int subgraph::write_kept_hs(ostream &out, const vector<bool> &keep)
{
	VVI &edges = hs.edges;
	vector<int> ecnts = hs.ecnts;
	vector<int> eucnts = hs.eucnts;
	assert(edges.size() == ecnts.size());
	assert(edges.size() == eucnts.size());
	assert(edges.size() == keep.size());
	for (int i = 0; i < keep.size(); i++)
	{
		if (!keep[i]) continue;
		vector<int> &el = edges[i];
		out << "hyper-edge (edges), counts=" << ecnts[i] << ", unique_counts=" << eucnts[i] << ", list=(";
		writev(out, el);
		out << ")\n";
	}
	return 0;
}

int subgraph::write_hss(ostream &out, const VVB &hss, const vector<int> &ucnts)
{
	for (int i = 0; i < hss.size(); i++)
	{
		for (int j = 0; j < hss[i].size(); j++)
		{
			if (hss[i][j]) out << "1";
			else out << "0";
		}
		out << ", unique_counts=" << ucnts[i] << "\n";
	}
	return 0;
}

