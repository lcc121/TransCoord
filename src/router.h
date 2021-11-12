/*
This file was modified from a file named router.h in Scallop
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

#ifndef __ROUTER_H__
#define __ROUTER_H__

#include <vector>
#include "util.h"
#include "splice_graph.h"
#include "equation.h"
#include "undirected_graph.h"
#include "hyper_set.h"

typedef pair<int, double> PID;
typedef map<int, double> MID;
typedef pair<PI, double> PPID;
typedef map<PI, double> MPID;

using namespace std;

class router
{
public:
	router(int r, splice_graph &g, MEI &ei, VE &ie);
	router(int r, splice_graph &g, MEI &ei, VE &ie, const MPII &mpi);
	router& operator=(const router &rt);

public:
	int root;					// central vertex
	splice_graph &gr;			// reference splice graph
	MEI &e2i;					// reference map of edge to index
	VE &i2e;					// reference map of index to edge
	vector<PI> routes;			// pairs of connections
	vector<int> counts;			// counts for routes

	MI e2u;						// edge to index
	vector<int> u2e;			// index to edge
	MED u2w;					// weights of edges of ug
	undirected_graph ug;		// bipartite graph

	int type;					// trivial, splitable, single, or multiple 
	int degree;					// level
	double ratio;				// worst ratio
	vector<equation> eqns;		// split results
	MPID pe2w;					// decompose results (for pairs of edges)

#ifdef USECLP
	MID se2w;
#endif

public:
	int classify();												// compute status
	int build();												// give solution

	// init
	int build_indices();										// build u2e and e2u
	int build_bipartite_graph();								// build bipartite graph
	vector<double> compute_balanced_weights();					// balanced weights

	// decompose splitable vertex
	int split();												// for splitable vertices

	// decompose unsplitable vertex with greedy algorithm
	int thread();												// for unsplitable vertices
	int thread_isolate1(int k, vector<double> &vw);
	int thread_isolate2(int k, vector<double> &vw);
	bool thread_leaf(vector<double> &vw);
	bool thread_turn(vector<double> &vw);

#ifdef USECLP
	// decompose unsplitable vertex with LP 
	int lpsolve();
	int extend_bipartite_graph_max();							// extended graph
	int extend_bipartite_graph_all();							// extended graph
	int build_maximum_spanning_tree();							// make ug a (maximum) spanning tree
	int decompose0_clp();										// solve LP with CLP
	int decompose1_clp();										// solve LP with CLP
	int decompose2_clp();										// solve LP with CLP
#endif

	// print and stats
	int print();
	int stats();
};

bool compare_edge_weight(const PED &x, const PED &y);

#endif
