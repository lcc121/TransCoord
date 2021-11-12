/*
This file was modified from a file named directed_graph.h in Scallop
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
#ifndef __DIRECTED_GRAPH_H__
#define __DIRECTED_GRAPH_H__

#include <vector>
#include <map>

#include "graph_base.h"

using namespace std;

class directed_graph : public graph_base
{
public:
	directed_graph();
	directed_graph(const directed_graph &gr);
	directed_graph& operator=(const directed_graph &gr);
	virtual ~directed_graph();

public:
	// modify the graph
	virtual edge_descriptor add_edge(int s, int t);
	virtual int remove_edge(edge_descriptor e);
	virtual int remove_edge(int s, int t);
	virtual int move_edge(edge_base *e, int x, int y);
	virtual int exchange(int x, int y, int z);
	virtual int rotate(int x, int y);

	// access functions
	virtual int in_degree(int v) const;
	virtual int out_degree(int v) const;
	virtual PEEI in_edges(int v);
	virtual PEEI out_edges(int v);

	// algorithms
	virtual int bfs_reverse(int t, vector<int> &v);
	virtual int bfs_reverse(int t, vector<int> &v, vector<int> &b);
	virtual int bfs_reverse(int t, set<edge_descriptor> &ss);
	virtual bool bfs_reverse(const vector<int> &vt, int s, const set<edge_descriptor> &fb);
	virtual bool compute_shortest_path(int x, int y, vector<int> &p);
	virtual bool compute_shortest_path(edge_descriptor ex, edge_descriptor ey, vector<int> &p);
	virtual bool check_path(int x, int y);
	virtual bool check_path(edge_descriptor ex, edge_descriptor ey);
	virtual bool intersect(edge_descriptor ex, edge_descriptor ey);
	virtual vector<int> topological_sort();
	virtual vector<int> topological_sort_reverse();
	virtual vector<int> topological_sort0();
	virtual int compute_in_partner(int x);
	virtual int compute_out_partner(int x);
	virtual int compute_in_equivalent_vertex(int x);
	virtual int compute_out_equivalent_vertex(int x);
	virtual int check_nest(int x, int r, set<edge_descriptor> &vv);
	virtual int check_nest(int x, int r, set<edge_descriptor> &vv, const vector<int> &tpo);
	virtual int check_nest(int x, int r, const vector<int> &tpo);

	// draw
	int draw(const string &file, const MIS &mis, const MES &mes, double len);
	int draw_v(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &colors);
	int draw_vwl(const string &file, const MIS &mis, const MES &mes, double len, const vector<int> &colors);
};

#endif
