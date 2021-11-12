/*
This file was modified from a file named graph_base.cc in Scallop
The original license info:
BSD 3-Clause License

Copyright (c) 2017, Mingfu Shao, Carl Kingsford, and Carnegie Mellon University
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include "graph_base.h"

#include <climits>
#include <cstdio>
#include <cassert>
#include <algorithm>

using namespace std;

graph_base::graph_base()
{}

graph_base::~graph_base()
{
	clear();
}

graph_base::graph_base(const graph_base &gr)
{
	//copy(gr); !!!
}

int graph_base::copy(const graph_base &gr)
{
	clear();
	for(int i = 0; i < gr.num_vertices(); i++) add_vertex();

	PEEI p = gr.edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		add_edge((*it)->source(), (*it)->target());
	}
	return 0;
}

int graph_base::add_vertex()
{
	vertex_base *v = new vertex_base();
	vv.push_back(v);
	return 0;
}

int graph_base::clear_vertex(int x)
{
	vector<edge_base*> v;
	PEEI pi = vv[x]->in_edges();
	PEEI po = vv[x]->out_edges();
	for(edge_iterator it = pi.first; it != pi.second; it++)
	{
		v.push_back(*it);
	}
	for(edge_iterator it = po.first; it != po.second; it++)
	{
		v.push_back(*it);
	}
	for(int i = 0; i < v.size(); i++)
	{
		remove_edge(v[i]);
	}
	return 0;
}

int graph_base::clear()
{
	for(int i = 0; i < vv.size(); i++) delete vv[i];
	for(edge_iterator it = se.begin(); it != se.end(); it++)
	{
		delete (*it);
	}
	vv.clear();
	se.clear();
	return 0;
}

int graph_base::degree(int v) const
{
	return vv[v]->degree();
}

PEB graph_base::edge(int s, int t) 
{
	assert(s >= 0 && s < vv.size());
	assert(t >= 0 && t < vv.size());
	PEEI p = vv[s]->out_edges();
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		int x = (*it)->target();
		if(x != t) continue;
		return PEB(*it, true);
	}
	return PEB(null_edge, false);
}

vector<edge_descriptor> graph_base::edges(int s, int t)
{
	vector<edge_descriptor> v;
	PEEI p = out_edges(s);
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		if((*it)->target() != t) continue;
		v.push_back(*it);
	}
	return v;
}

PEEI graph_base::edges() const
{
	return PEEI(se.begin(), se.end());
}

set<int> graph_base::adjacent_vertices(int s)
{
	set<int> ss;
	PEEI p = out_edges(s);
	for(edge_iterator it = p.first; it != p.second; it++)
	{
		assert((*it)->source() == s);
		int t = (*it)->target();
		if(ss.find(t) == ss.end()) ss.insert(t);
	}
	return ss;
}

size_t graph_base::support_size() const
{
	int s = 0;
	for(int i = 0; i < num_vertices(); i++)
	{
		if(degree(i) == 0) continue;
		s++;
	}
	return s;
}

size_t graph_base::num_vertices() const
{
	return vv.size();
}

size_t graph_base::num_edges() const
{
	return se.size();
}

size_t graph_base::num_inner_edges() const
{
	PEEI pei = edges();
	edge_iterator it1 = pei.first, it2 = pei.second;
	size_t num = 0;
	for (; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		if (s == 0 || t == num_vertices() - 1) continue;
		num += 1;
	}
	return num;
}

int graph_base::get_edge_indices(VE &i2e, MEI &e2i)
{
	i2e.clear();
	e2i.clear();
	int index = 0;
	PEEI pei = edges();
	edge_iterator it1 = pei.first, it2 = pei.second;
	for(; it1 != it2; it1++)
	{
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	return 0;
}

int graph_base::get_inner_edge_indices(VE &i2e, MEI &e2i)
{
	// give index to all edges except ss and tt
	i2e.clear();
	e2i.clear();
	int index = 0;
	PEEI pei = edges();
	edge_iterator it1 = pei.first, it2 = pei.second;
	for (; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		if (s == 0 || t == num_vertices() - 1) continue;
		e2i.insert(PEI(*it1, index));
		i2e.push_back(*it1);
		index++;
	}
	int nie = num_inner_edges();
	assert(i2e.size() == nie);
	return 0;
}

bool graph_base::bfs(const vector<int> &vs, int t, const set<edge_descriptor> &fb)
{
	set<int> closed(vs.begin(), vs.end());
	vector<int> open = vs;
	int p = 0;
	while(p < open.size())
	{
		int x = open[p];
		if(x == t) return true;
		p++;
		PEEI pei = out_edges(x);
		edge_iterator it1 = pei.first, it2 = pei.second;
		for(; it1 != it2; it1++)
		{
			if(fb.find(*it1) != fb.end()) continue;
			int y = (*it1)->target();
			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);
		}
	}
	return false;
}

int graph_base::bfs(int s, vector<int> &v, vector<int> &b)
{
	v.clear();
	b.clear();
	vector<int> open;
	open.push_back(s);
	set<int> closed;
	closed.insert(s);
	v.push_back(s);
	b.push_back(-1);
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		p++;
		PEEI pei = out_edges(x);
		edge_iterator it1 = pei.first, it2 = pei.second;
		for(; it1 != it2; it1++)
		{
			int y = (*it1)->target();
			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);
			v.push_back(y);
			b.push_back(x);
		}
	}
	return 0;
}

int graph_base::bfs(int s, set<edge_descriptor> &ss)
{
	ss.clear();
	set<int> closed;
	vector<int> open;
	open.push_back(s);
	closed.insert(s);
	int p = 0;

	while(p < open.size())
	{
		int x = open[p];
		p++;
		PEEI pei = out_edges(x);
		edge_iterator it1 = pei.first, it2 = pei.second;
		for(; it1 != it2; it1++)
		{
			int y = (*it1)->target();
			ss.insert(*it1);
			if(closed.find(y) != closed.end()) continue;
			closed.insert(y);
			open.push_back(y);
		}
	}
	return 0;
}

int graph_base::bfs(int t, vector<int> &v) 
{
	vector<int> b;
	return bfs(t, v, b);
}

bool graph_base::check_path(int s, int t) 
{
	set<edge_descriptor> fb;
	vector<int> ss;
	ss.push_back(s);
	return bfs(ss, t, fb);
}

bool graph_base::compute_shortest_path(int s, int t, vector<int> &p) 
{
	p.clear();
	vector<int> v;
	vector<int> b;
	bfs(s, v, b);

	int n = v.size() - 1;
	int x = t;
	while(true)
	{
		p.push_back(x);
		if(x == s) break;

		int k;
		for(k = n; k >= 0; k--)
		{
			if(v[k] == x) break;	
		}
		if(k <= -1) return false;

		x = b[k];
		n = k - 1;
		if(n < 0) break;
	}
	reverse(p.begin(), p.end());
	return true;
}

bool graph_base::check_nested() 
{
	PEEI p = edges();
	for(edge_iterator i = p.first; i != p.second; i++)
	{
		edge_iterator j = i;
		j++;
		for(; j != p.second; j++)
		{
			bool b = intersect(*i, *j);
			if(b == true) return false;
		}
	}
	return true;
}

int graph_base::print() const
{
	printf("total %lu vertices, %lu edges\n", vv.size(), se.size());
	for(int i = 0; i < vv.size(); i++)
	{
		printf("vertex %d: ", i);
		vv[i]->print();
	}

	for(edge_iterator it = se.begin(); it != se.end(); it++)
	{
		(*it)->print();
	}
	return 0;
}

