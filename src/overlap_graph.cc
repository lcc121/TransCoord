#include "overlap_graph.h"

#include <assert.h>

overlap_graph::overlap_graph()
{}

overlap_graph::overlap_graph(const hyper_set &hs)
{
	clear();
	// add vertices
	const VVI &edges = hs.edges;
	const vector<int> &uppcnts = hs.uppcnts;
	int j = 0;
	for (int i = 0; i < edges.size(); i++)
	{
		const vector<int> &el = edges[i];
		double count = uppcnts[i];
		if (count <= 0) continue;
		add_vertex();
		for (int k = 0; k < el.size(); k++) assert(el[k] != -1);
		set_vertex_edge_list(j, el);
		set_vertex_weight(j, count);
		j++;
	}
	build_index();
	add_edges();
	if (verbose >= 2) write(cout);
	contract();
	if (verbose >= 2)
	{
		cout << "contracted ";
		write(cout);
	}
}

overlap_graph::~overlap_graph()
{}

int overlap_graph::contract()
{
	// simplify overlap_graph by contracting consecutive nodes in a linear path.
	// find node with only 1 out-edge
	if (verbose >= 5) cout << "into overlap_graph contract." << endl;
	int vn = num_vertices();
	for (int i = 0; i < vn; i++)
	{
		int od = out_degree(i);
		if (od != 1) continue;
		contract_from(i);
	}
	return 0;
}

int overlap_graph::extract_paths()
{
	paths.clear();
	set<edge_descriptor> used;
	int iter = 0;
	while (true)
	{
		iter += 1;
		PEB seed = find_seed_edge(used);
		if (seed.second == false) break;
		int s = seed.first->source();
		int t = seed.first->target();
		//PEB e = edge(s, t);
		vector<int> v;
		v.push_back(s);
		v.push_back(t);
		if (verbose >= 5)
		{
			cout << "overlap_graph seed=";
			writev(cout, v);
			cout << endl;
		}
		while (true)
		{
			s = left_extend(s);
			if (s == -1) break;
			v.insert(v.begin(), s);
		}
		while (true)
		{
			t = right_extend(t);
			if (t == -1) break;
			v.push_back(t);
		}
		assert(v.size() <= vv.size());
		//set<int> es;
		vector<int> es;
		for (int i = 0; i < v.size(); i++)
		{
			vector<int> &el = v2el[v[i]];
			//es.insert(el.begin(), el.end());
			if (es.empty()) es.insert(es.end(), el.begin(), el.end());
			else 
			{
				int j;
				for (j = 0; j < el.size(); j++) { if (el[j] == es.back()) break; }
				for (j = j + 1; j < el.size(); j++) { es.push_back(el[j]); }
			}
		}
		path p;
		p.v.assign(es.begin(), es.end());
		p.abd = get_path_abundance(v);
		if (verbose >= 5)
		{
			cout << "seed_lr_extend=";
			writev(cout, v);
			cout << "=" << p.abd << endl;
		}
		paths.push_back(p);
		update(v, p.abd, used);
		if (verbose >= 6)
		{
			write(cout);
			int in;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	if (verbose >= 2) cout << "overlap_graph iterated " << iter << " times." << endl;
	if (verbose >= 5)
	{
		cout << "found " << paths.size() << " paths:\n";
		for (int i = 0; i < paths.size(); i++) { paths[i].print(i); }
	}
	return 0;
}

int overlap_graph::clear()
{
	directed_graph::clear();
	v2el.clear();
	vwrt.clear();
	ewrt.clear();
	ppe2v.clear();
	return 0;
}


int overlap_graph::set_vertex_edge_list(int v, const vector<int> &el)
{
	assert(v >= 0 && v < vv.size());
	if (v2el.size() != vv.size()) v2el.resize(vv.size());
	v2el[v] = el;
	return 0;
}

int overlap_graph::set_vertex_weight(int v, double w)
{
	assert(v >= 0 && v < vv.size());
	if (vwrt.size() != vv.size()) vwrt.resize(vv.size());
	vwrt[v] = w;
	return 0;
}

double overlap_graph::get_vertex_weight(int v) const
{
	assert(v >= 0 && v < vwrt.size());
	return vwrt[v];
}

int overlap_graph::set_edge_weight(edge_descriptor e, double w)
{
	if (ewrt.find(e) != ewrt.end()) ewrt[e] = w;
	else ewrt.insert(PED(e, w));
	return 0;
}

double overlap_graph::get_edge_weight(edge_descriptor e) const
{
	MED::const_iterator it = ewrt.find(e);
	assert(it != ewrt.end());
	return it->second;
}


int  overlap_graph::add_vertices(const hyper_set &hs)
{

}

int overlap_graph::add_necessary_vertices(const hyper_set &hs)
{

}

int overlap_graph::build_index()
{
	ppe2v.clear();
	for (int i = 0; i < v2el.size(); i++)
	{
		vector<int> &el = v2el[i];
		for (int j = 0; j < el.size(); j++)
		{
			int e = el[j];
			if (ppe2v.find(e) == ppe2v.end())
			{
				set<int> s;
				s.insert(i);
				ppe2v.insert(PISI(e, s));
			}
			else ppe2v[e].insert(i);
		}
	}
	return 0;
}

int overlap_graph::add_edges()
{
	for (int i = 0; i < vv.size(); i++)
	{
		vector<int> &el = v2el[i];
		// find vertice link to this one
		if (ppe2v.find(el[0]) == ppe2v.end()) continue;
		set<int> &ss = ppe2v[el[0]];
		set<int>::iterator it;
		for (it = ss.begin(); it != ss.end(); ++it)
		{
			int s = *it;
			if (s == i) continue;
			PBI b = compatible(s, i);
			if (b.first)
			{
				edge_descriptor e = add_edge(s, i);
				set_edge_weight(e, b.second);
			}
		}
	}
	return 0;
}

PBI overlap_graph::compatible(int x, int y)
{
	// if x compatibly link to y
	assert(x != y);
	bool b = false;
	int l = 0;
	vector<int> &xel = v2el[x];
	vector<int> &yel = v2el[y];
	int i = 0, j = 0;
	for (; i < xel.size(); i++) { if (xel[i] == yel[0]) break; }
	assert(i < xel.size());
	while (i < xel.size() && j < yel.size())
	{
		if (xel[i] != yel[j]) break;
		l++;
		i++;
		j++;
	}
	if (i < xel.size() || j == yel.size()) b = false;
	else b = true;
	return PBI(b, l);
}


PEB overlap_graph::find_seed_edge(set<edge_descriptor> &used)
{
	PEEI pei = edges();
	if (pei.first == pei.second) return PEB(null_edge, false);
	double max = 0;
	edge_iterator it;
	for (it = pei.first; it != pei.second; ++it)
	{
		if (used.find(*it) != used.end()) continue;
		double ew = get_edge_weight(*it);
		assert(ew >= 1);
		if (ew <= max) continue;
		max = ew;
	}
	if (max <= 0) return PEB(null_edge, false);
	double max_min = 0;
	edge_descriptor e;
	for (it = pei.first; it != pei.second; ++it)
	{
		if (used.find(*it) != used.end()) continue;
		double ew = get_edge_weight(*it);
		assert(ew >= 1);
		if (ew != max) continue;
		int s = (*it)->source();
		int t = (*it)->target();
		double sw = get_vertex_weight(s);
		double tw = get_vertex_weight(t);
		double min_st = (sw <= tw) ? sw : tw;
		if (min_st <= max_min) continue;
		max_min = min_st;
		e = (*it);
	}
	assert(max_min > 0);
	assert(e != null_edge);
	return PEB(e, true);
}

int overlap_graph::left_extend(int v)
{
	PEEI pi = in_edges(v);
	if (pi.first == pi.second) return -1;
	edge_iterator it;
	double max_ew = 0;
	for (it = pi.first; it != pi.second; ++it)
	{
		double ew = get_edge_weight(*it);
		if (ew <= max_ew) continue;
		max_ew = ew;
	}
	assert(max_ew >= 1);

	double max_vw = 0;
	int lv = -1;
	for (it = pi.first; it != pi.second; ++it)
	{
		double ew = get_edge_weight(*it);
		if (ew != max_ew) continue;
		int s = (*it)->source();
		assert((*it)->target() == v);
		double vw = get_vertex_weight(s);
		if (vw <= max_vw) continue;
		max_vw = vw;
		lv = s;
	}
	assert(max_vw >= min_cov);
	assert(lv != -1);
	return lv;
}

int overlap_graph::right_extend(int v)
{
	PEEI po = out_edges(v);
	if (po.first == po.second) return -1;
	edge_iterator it;
	double max_ew = 0;
	for (it = po.first; it != po.second; ++it)
	{
		double ew = get_edge_weight(*it);
		if (ew <= max_ew) continue;
		max_ew = ew;
	}
	assert(max_ew >= 1);

	double max_vw = 0;
	int rv = -1;
	for (it = po.first; it != po.second; ++it)
	{
		double ew = get_edge_weight(*it);
		if (ew != max_ew) continue;
		assert((*it)->source() == v);
		int t = (*it)->target();
		double vw = get_vertex_weight(t);
		if (vw <= max_vw) continue;
		max_vw = vw;
		rv = t;
	}
	assert(max_vw >= min_cov);
	assert(rv != -1);
	return rv;
}

double overlap_graph::get_path_abundance(const vector<int> &v)
{
	assert(v.size() >= 2);
	double min = get_vertex_weight(v[0]);
	for (int i = 0; i < v.size(); i++)
	{
		double vw = get_vertex_weight(v[i]);
		if (vw < min) min = vw;
	}
	return min;
}

int overlap_graph::update(vector<int> &v, double w, set<edge_descriptor> &used)
{
	assert(w > 0);
	for (int i = 0; i < v.size(); i++)
	{
		int vi = v[i];
		double vw = get_vertex_weight(vi);
		vw -= w;
		set_vertex_weight(vi, vw);
		if (vw >= min_cov)
		{
			if (i > 0)
			{
				PEB e = edge(v[i - 1], v[i]);
				if (e.second) used.insert(e.first);
			}
			continue;
		}
		PEEI pi = in_edges(vi);
		PEEI po = out_edges(vi);
		edge_iterator it;
		for (it = pi.first; it != pi.second; ++it)
		{
			edge_iterator it1 = used.find(*it);
			if (it1 != used.end()) used.erase(it1);
		}
		for (it = po.first; it != po.second; ++it)
		{
			edge_iterator it1 = used.find(*it);
			if (it1 != used.end()) used.erase(it1);
		}
		clear_vertex(vi);
	}
	return 0;
}


int overlap_graph::contract_from(int x)
{
	int t = x;
	int od = out_degree(x);
	assert(od == 1);
	int id = 1;
	vector<int> v;
	while (od == 1 && id == 1) 
	{
		v.push_back(t);
		PEEI po = out_edges(t);
		edge_iterator it = po.first;
		assert(++it == po.second);
		it = po.first;
		t = (*it)->target();
		od = out_degree(t);
		id = in_degree(t);
	}
	if (id == 1) v.push_back(t);
	if (v.size() <= 1) return 0;
	if (verbose >= 5)
	{
		cout << "contract vertices: ";
		writev(cout, v);
		cout << endl;
	}
	// link x->(v.back()->)
	PEEI po = out_edges(v.back());
	edge_iterator it;
	for (it = po.first; it != po.second; ++it)
	{
		t = (*it)->target();
		double ew = get_edge_weight(*it);
		edge_descriptor e = add_edge(x, t);
		set_edge_weight(e, ew);
	}
	double xw = get_path_abundance(v);
	set_vertex_weight(x, xw);
	//set<int> s;
	vector<int> es;
	// clear inner vertices
	for (int i = 0; i < v.size(); i++)
	{
		int vi = v[i];
		vector<int> &el = v2el[vi];
		//s.insert(el.begin(), el.end());
		if (es.empty()) es.insert(es.end(), el.begin(), el.end());
		else
		{
			int j;
			for (j = 0; j < el.size(); j++) { if (el[j] == es.back()) break; }
			for (j = j + 1; j < el.size(); j++) { es.push_back(el[j]); }
		}
		if (i == 0) continue;
		clear_vertex(vi);
	}
	//vector<int> el(s.begin(), s.end());
	//set_vertex_edge_list(x, el);
	set_vertex_edge_list(x, es);
	if (verbose >= 6)
	{
		write(cout);
		int in;
		cout << "please input an integer:" << endl;
		cin >> in;
	}
	return 0;
}


int overlap_graph::write(ostream &out)
{
	int vn = num_vertices();
	int en = num_edges();
	out << "overlap_graph: " << vn << " vertices, " << en << " edges.\n";
	if (verbose <= 2) return 0;
	out << vn << " vertices(v2el):\n";
	for (int i = 0; i < v2el.size(); i++)
	{
		out << i << "=(";
		writev(out, v2el[i]);
		double vw = get_vertex_weight(i);
		out << ")=" << vw << "\n";
	}
	out << num_edges() << " edges:";
	PEEI pei = edges();
	edge_iterator it1 = pei.first, it2 = pei.second;
	int pre_s = -1;
	for (; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double ew = get_edge_weight(*it1);
		if (s != pre_s) out << "\n";
		pre_s = s;
		out << s << "->" << t << ":" << ew << " ";
	}
	out << "\n";
	return 0;
}
