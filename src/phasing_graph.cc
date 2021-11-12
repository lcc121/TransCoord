#include "phasing_graph.h"

phasing_graph::phasing_graph()
{}

phasing_graph::phasing_graph(splice_graph &gr, VE &i2e, MEI &e2i, vector<path> &pp, hyper_set &hs)
{
	add_normal_nodes(gr, i2e, e2i, pp, hs);
	add_contracted_nodes(pp);
	if (verbose >= 2) write(cout);
}

phasing_graph::~phasing_graph()
{}

int phasing_graph::extract_paths(splice_graph &gr, VE &i2e)
{
	paths.clear();
	set<int> c_used, n_used;
	int iter = 0;
	while (true)
	{
		iter += 1;
		int seed = find_seed(c_used, n_used);
		if (seed == -1) break;
		if (verbose >= 5) cout << "phasing_graph seed=" << seed << endl;
		bool invalid = false;// found path not from ss~>tt
		int s = seed;
		vector<int> v(1, s);
		// left extend
		while (true)
		{
			s = left_extend(s);
			if (s == -2) invalid = true;
			if (s == -1 || s == -2) break;
			v.insert(v.begin(), s);
		}
		if (invalid)
		{
			if (seed < pos) n_used.insert(seed);
			else c_used.insert(seed);
			if (verbose >= 5)
			{
				cout << "invalid seed " << seed << ": ";
				writev(cout, v);
				cout << endl;
			}
			continue;
		}
		// right extend
		s = seed;
		while (true)
		{
			s = right_extend(s);
			if (s == -2) invalid = true;
			if (s == -1 || s == -2) break;
			v.push_back(s);
		}
		if (invalid)
		{
			if (seed < pos) n_used.insert(seed);
			else c_used.insert(seed);
			if (verbose >= 5)
			{
				cout << "invalid seed " << seed << ": ";
				writev(cout, v);
				cout << endl;
			}
			continue;
		}
		// extract path and update graph
		//set<int> es;
		vector<int> es;
		for (int i = 0; i < v.size(); i++)
		{
			vector<int> &el = v2el[v[i]];
			//es.insert(el.begin(), el.end());
			for (int j = 0; j < el.size(); j++)
			{
				if (j == 0 && !es.empty() && el[j] == es.back()) continue;
				es.push_back(el[j]);
			}
		}
		path p;
		p.v.assign(es.begin(), es.end());
		p.abd = get_path_abundance(v);
		if (verbose >= 5)
		{
			cout << "seed_lr_extend(vertices)=";
			writev(cout, v);
			cout << "=" << p.abd << endl;
			cout << "seed_lr_extend(el)=";
			writev(cout, p.v);
			cout << "=" << p.abd << endl;
		}
		paths.push_back(p);
		update(v, p.abd, c_used, n_used, i2e, gr.num_vertices());
		if (verbose >= 6)
		{
			write(cout);
			int in;
			cout << "please input an integer:" << endl;
			cin >> in;
		}
	}
	if (verbose >= 2) cout << "phasing_graph iterated " << iter << "times." << endl;
	if (verbose >= 5)
	{
		cout << "\n" << paths.size() << " paths:\n";
		for (int i = 0; i < paths.size(); i++) { paths[i].print(i); }
	}
	return 0;
}


int phasing_graph::set_vertex_edge_list(int v, vector<int> &el)
{
	assert(v >= 0 && v < vv.size());
	if (v2el.size() != vv.size()) v2el.resize(vv.size());
	v2el[v] = el;
	return 0;
}

int phasing_graph::set_vertex_weight(int v, double w)
{
	assert(v >= 0 && v < vv.size());
	if (vwrt.size() != vv.size()) vwrt.resize(vv.size());
	vwrt[v] = w;
	return 0;
}

double phasing_graph::get_vertex_weight(int v) const
{
	assert(v >= 0 && v < vwrt.size());
	return vwrt[v];
}

int phasing_graph::set_edge_weight(edge_descriptor e, int w)
{
	if (ewrt.find(e) != ewrt.end()) ewrt[e] = w;
	else ewrt.insert(PEI(e, w));
	return 0;
}

int phasing_graph::get_edge_weight(edge_descriptor e) const
{
	MEI::const_iterator it = ewrt.find(e);
	assert(it != ewrt.end());
	return it->second;
}


int phasing_graph::add_normal_nodes(splice_graph &gr, VE &i2e, MEI &e2i, vector<path> &pp, hyper_set &hs)
{
	pos = i2e.size();
	int vn = gr.num_vertices();
	// add vertices
	for (int i = 0; i < i2e.size(); i++)
	{
		add_vertex();
		double w = gr.get_edge_weight(i2e[i]);
		set_vertex_weight(i, w);
		vector<int> v(1, i);
		set_vertex_edge_list(i, v);
	}
	// add edges
	for (int i = 0; i < vn; i++)
	{
		vector<int> in, out;          // in/out edges of each vertex
		int max_in = -1, max_out = -1;// max in/out edge index of each vertex
		double max1 = 0, max2 = 0;    // max in/out edge weight
		PEEI pei = gr.in_edges(i);
		edge_iterator it;
		for (it = pei.first; it != pei.second; ++it)
		{
			//int ei = subg.get_edge_index(*it);
			assert(e2i.find(*it) != e2i.end());
			int ei = e2i.find(*it)->second;
			assert(ei >= 0 && ei < pos);
			in.push_back(ei);
			double ew = gr.get_edge_weight(*it);
			if (ew <= max1) continue;
			max1 = ew;
			max_in = ei;
		}
		pei = gr.out_edges(i);
		for (it = pei.first; it != pei.second; ++it)
		{
			//int ei = subg.get_edge_index(*it);
			assert(e2i.find(*it) != e2i.end());
			int ei = e2i.find(*it)->second;
			assert(ei >= 0 && ei < pos);
			out.push_back(ei);
			double ew = gr.get_edge_weight(*it);
			if (ew <= max2) continue;
			max2 = ew;
			max_out = ei;
		}
		if (in.size() == 0 || out.size() == 0) continue;
		for (int j = 0; j < in.size(); j++)
		{
			int s = in[j];
			int ss = i2e[s]->source();
			for (int k = 0; k < out.size(); k++)
			{
				int t = out[k];
				edge_descriptor e = add_edge(s, t);
				int tt = i2e[t]->target();
				if (ss == 0 || tt == vn - 1) set_edge_weight(e, 1);
				else set_edge_weight(e, 0);
			}
		}
		// heaviest in/out-edge of each vertex =1
		PEB e = edge(max_in, max_out);
		assert(e.second);
		set_edge_weight(e.first, 1);
	}
	// update edge weight
	for (int i = 0; i < pp.size(); i++)
	{
		path &p = pp[i];
		vector<int> &v = p.v;
		for (int j = 1; j < v.size(); j++)
		{
			int s = v[j - 1];
			int t = v[j];
			PEB e = edge(s, t);
			assert(e.second);
			set_edge_weight(e.first, 1);
		}
	}
	VVI &edges = hs.edges;
	vector<int> &eucnts = hs.eucnts;
	for (int i = 0; i < eucnts.size(); i++)
	{
		if (eucnts[i] <= 0) continue;
		vector<int> &el = edges[i];
		for (int j = 1; j < el.size(); j++)
		{
			int s = el[j - 1];
			int t = el[j];
			PEB e = edge(s, t);
			assert(e.second);
			set_edge_weight(e.first, 1);
		}
	}

	return 0;
}

int phasing_graph::add_contracted_nodes(vector<path> &pp)
{
	// add vertices and update vertex weight
	int index = pos;
	for (int i = 0; i < pp.size(); i++)
	{
		path &p = pp[i];
		vector<int> &v = p.v;
		//if (p.v.size() <= 2) continue;
		add_vertex();
		double w = p.abd;
		set_vertex_weight(index, w);
		set_vertex_edge_list(index, v);
		for (int j = 1; j < (int)v.size() - 1; j++)// no change of first and last
		{
			int vi = v[j];
			double vw = get_vertex_weight(vi);
			vw -= w;
			set_vertex_weight(vi, vw);
		}
		index += 1;
	}
	// add edges 
	for (int i = pos; i < vv.size(); i++)
	{
		vector<int> &el = v2el[i];
		int s = el[0];
		int t = el.back();
		assert(s >= 0 && s < pos);
		assert(t >= 0 && t < pos);
		edge_descriptor e1 = add_edge(s, i);
		set_edge_weight(e1, 1);
		edge_descriptor e2 = add_edge(i, t);
		set_edge_weight(e2, 1);
	}
	return 0;
}


int phasing_graph::find_seed(set<int> &c_used, set<int> &n_used)
{
	int max = 0;
	int max_vi = -1;
	int vn = num_vertices();
	int nvn = pos;
	int cvn = vn - pos;
	if (c_used.size() < cvn)
	{
		for (int i = pos; i < vn; i++)
		{
			if (c_used.find(i) != c_used.end()) continue;
			double vw = get_vertex_weight(i);
			if (vw <= max) continue;
			max = vw;
			max_vi = i;
		}
	}
	else if (n_used.size() < nvn)
	{
		for (int i = 0; i < pos; i++)
		{
			if (n_used.find(i) != n_used.end()) continue;
			double vw = get_vertex_weight(i);
			if (vw <= max) continue;
			max = vw;
			max_vi = i;
		}
	}
	
	return max_vi;
}

int phasing_graph::left_extend(int v)
{
	double max0 = 0.1, max1 = 0.1;
	int v0 = -2, v1 = -2;
	PEEI pei = in_edges(v);
	if (pei.first == pei.second) return -1;// reach terminal s
	edge_iterator it;
	for (it = pei.first; it != pei.second; ++it)
	{
		int ew = get_edge_weight(*it);
		assert((*it)->target() == v);
		int s = (*it)->source();
		double sw = get_vertex_weight(s); 
		if (ew == 0)
		{
			if (sw > max0) { max0 = sw; v0 = s; }
		}
		else
		{
			if (sw > max1) { max1 = sw; v1 = s; }
		}
	}
	if (v1 != -2) return v1;
	if (v0 != -2) return v0;
	return -2;
}

int phasing_graph::right_extend(int v)
{
	double max0 = 0.1, max1 = 0.1;
	int v0 = -2, v1 = -2;
	PEEI pei = out_edges(v);
	if (pei.first == pei.second) return -1;
	edge_iterator it;
	for (it = pei.first; it != pei.second; ++it)
	{
		int ew = get_edge_weight(*it);
		assert((*it)->source() == v);
		int t = (*it)->target();
		double tw = get_vertex_weight(t);
		if (ew == 0)
		{
			if (tw > max0) { max0 = tw; v0 = t; }
		}
		else
		{
			if (tw > max1) { max1 = tw; v1 = t; }
		}
	}
	if (v1 != -2) return v1;
	if (v0 != -2) return v0;
	return -2;
}

double phasing_graph::get_path_abundance(const vector<int> &v)
{
	assert(v.size() > 0);
	double min = get_vertex_weight(v[0]);
	for (int i = 0; i < v.size(); i++)
	{
		double vw = get_vertex_weight(v[i]);
		if (vw < min) min = vw;
	}
	return min;
}

int phasing_graph::update(const vector<int> &v, double w, set<int> &c_used, set<int> &n_used, const VE &i2e, int vn)
{
	bool cleared = false;
	for (int i = 0; i < v.size(); i++)
	{
		int vi = v[i];
		if (vi < pos) n_used.insert(vi);
		else c_used.insert(vi);
		double vw = get_vertex_weight(vi);
		vw -= w;
		set_vertex_weight(vi, vw);
		if (vw >= min_cov) continue;
		clear_vertex(vi);
		cleared = true;
	}
	// clear vertices unreachable from ss/tt
	while (cleared)
	{
		cleared = false;
		for (int i = 0; i < vv.size(); i++)
		{
			PEEI pei = in_edges(i);
			PEEI peo = out_edges(i);
			if (pei.first != pei.second && peo.first != peo.second) continue;
			if (pei.first == pei.second && peo.first == peo.second)
			{
				//already cleared or new isolated sstt-edge
				if (i < pos) n_used.insert(i);
				else c_used.insert(i);
				continue;
			}
			if (pei.first == pei.second)
			{
				if (i < pos)
				{
					int s = i2e[i]->source();
					if (s == 0) continue;
				}
			}
			else
			{
				if (i < pos)
				{
					int t = i2e[i]->target();
					if (t == vn - 1) continue;
				}
			}
			clear_vertex(i);// can not used as seed
			if (i < pos) n_used.insert(i);
			else c_used.insert(i);
			cleared = true;
		}
	}
	return 0;
}

int phasing_graph::write(ostream &out)
{
	int vn = num_vertices();
	int en = num_edges();
	out << "phasing_graph: " << vn << " vertices, " << en << " edges.\n";
	if (verbose <= 2) return 0;
	out << vn << " vertices(v2el):\n";
	for (int i = 0; i < vn; i++)
	{
		out << i << "=(";
		writev(out, v2el[i]);
		double vw = get_vertex_weight(i);
		out << ")=" << vw << "\n";
	}
	out << num_edges() << " edges:";
	PEEI pei = edges();
	edge_iterator it;
	int pre_s = -1;
	for (it = pei.first; it != pei.second; ++it)
	{
		int s = (*it)->source();
		int t = (*it)->target();
		int ew = get_edge_weight(*it);
		if (s != pre_s) cout << "\n";
		pre_s = s;
		out << s << "->" << t << ":" << ew << " ";
	}
	out << "\n";
	return 0;
}

