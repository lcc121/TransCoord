/*
This file was modified from a file named scallop.h in Scallop
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

#ifndef __SCALLOP3_H__
#define __SCALLOP3_H__

#include "splice_graph.h"
#include "hyper_set.h"
#include "equation.h"
#include "router.h"
#include "path.h"

typedef map< edge_descriptor, vector<int> > MEV;
typedef pair< edge_descriptor, vector<int> > PEV;
typedef pair< vector<int>, vector<int> > PVV;
typedef pair<PEE, int> PPEEI;
typedef map<PEE, int> MPEEI;
typedef pair<int, int> PI;
typedef map<int, int> MI;

// for noisy splice graph
class scallop
{
public:
	scallop();
	scallop(const splice_graph &gr, const hyper_set &hs);
	virtual ~scallop();

public:
	int assemble();

public:
	splice_graph gr;					// splice graph
	MEI e2i;							// edge map, from edge to index
	VE i2e;								// edge map, from index to edge
	MEV mev;							// super edges
	vector<int> v2v;					// vertex map
	hyper_set hs;						// hyper edges
	int round;							// iteration
	set<int> nonzeroset;				// vertices with degree >= 1
	vector<path> paths;					// predicted paths

private:
	// init
	int classify();
	int init_vertex_map();
	int init_super_edges();
	int init_inner_weights();
	int init_nonzeroset();
	int add_pseudo_hyper_edges();
	int refine_splice_graph();

	// resolve iteratively
	bool resolve_trivial_vertex(int type, double jump_ratio);
	bool resolve_trivial_vertex_fast(double jump_ratio);
	bool resolve_single_trivial_vertex_fast(int i, double jump_ratio);
	bool resolve_smallest_edges(double max_ratio);
	bool resolve_negligible_edges(bool extend, double max_ratio);
	bool resolve_splittable_vertex(int type, int degree, double max_ratio);
	bool resolve_unsplittable_vertex(int type, int degree, double max_ratio);
	bool resolve_hyper_edge(int fsize);

	// smooth vertex
	int balance_vertex(int x);
	double compute_balance_ratio(int x);

	// decomposing subroutines
	int compute_smallest_edge(int x, double &ratio);
	int decompose_trivial_vertex(int v);
	int decompose_vertex_extend(int v, MPID &pe2w);
	int decompose_vertex_replace(int v, MPID &pe2w);
	int classify_trivial_vertex(int v, bool fast);
	int exchange_sink(int old_sink, int new_sink);
	int split_vertex(int x, const vector<int> &xe, const vector<int> &ye);
	int split_edge(int exi, double w);
	int merge_adjacent_edges(int x, int y);
	int merge_adjacent_edges(int x, int y, double ww);
	int merge_adjacent_equal_edges(int x, int y);
	int remove_edge(int e);
	int split_merge_path(const VE &p, double w);
	int split_merge_path(const vector<int> &p, double w);
	int collect_existing_st_paths();
	int collect_path(int e);
	int compute_length(const path &p);
	int greedy_decompose();

	// stats, print, and draw
	int print();
	int stats();
	int summarize_vertices();
	int draw_splice_graph(const string &file);
	vector<int> topological_sort();
};

#endif

