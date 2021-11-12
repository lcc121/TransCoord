/*
This file was modified from a file named hyper_set.h in Scallop
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
#ifndef __HYPER_SET_H__
#define __HYPER_SET_H__

#include <map>
#include <set>
#include <vector>

#include "util.h"
#include "directed_graph.h"
#include "splice_graph.h"

using namespace std;

typedef pair<int, int> PII;
typedef pair<vector<int>, pair<int, int>> PVIPII;
typedef map<vector<int>, pair<int, int>> MVIPII;
struct node_info
{
	int count;  // read count
	int ucnt;   // unique-read count
	int uppcnt; // paired-path count
};
typedef pair<vector<int>, node_info> PVINI;
typedef map<vector<int>, node_info> MVINI;
typedef map< int, set<int> > MISI;
typedef pair< int, set<int> > PISI;
typedef vector< vector<int> > VVI;
typedef map< pair<int, int>, int> MPII;
typedef pair< pair<int, int>, int> PPII;

class hyper_set
{
public:
	//MVIPII nodes;	    // hyper-edges using list-of-nodes
	MVINI nodes;

	VVI edges;			// hyper-edges using list-of-edges
	vector<int> ecnts;	// counts for edges
	vector<int> eucnts; // unique map counts for edges
	vector<int> uppcnts;// unique paired-path count
	MISI e2s;			// index: from edge to hyper-edges

public:
	int clear();
	int add_node_list(const set<int> &s, int uc, int upe);
	int add_node_list(const set<int> &s, int c, int uc, int upe);
	int add_node_list(const vector<int> &s, int c, int uc, int upe);
	int build(directed_graph &gr, MEI &e2i);
	int build_edges(directed_graph &gr, MEI &e2i);
	int build_index();
	int update_index();
	set<int> get_intersection(const vector<int> &v);
	MI get_successors(int e);
	MI get_predecessors(int e);
	int get_heaviest_successor(const set<int> &s, int e);
	int get_heaviest_predecessor(const set<int> &s, int e);
	MI get_compatible_predecessors(vector<int> &e);
	MI get_compatible_successors(vector<int> &e);
	int get_heaviest_compatible_predecessor(const set<int> &s, vector<int> &e);
	int get_heaviest_compatible_successor(const set<int> &s, vector<int> &e);
	int get_heaviest_exist_compatible(MI &m, const set<int> &s, vector<int> &e);
	MPII get_routes(int x, directed_graph &gr, MEI &e2i);
	MEI get_single_edge_counts(splice_graph &gr, VE &i2e);
	int print();
	int write(ostream & out);
	int write_edges(ostream & out);

public:
	int replace(int x, int e);
	int replace(int x, int y, int e);
	int replace(const vector<int> &x, int e);
	int remove(int e);
	int remove(const vector<int> &x);
	int remove(const set<int> &x);
	int remove_pair(int x, int y);
	int insert_between(int x, int y, int e);
	bool useful(const vector<int> &v, int k1, int k2);
	bool extend(int e);
	bool left_extend(int e);
	bool left_extend(const vector<int> &s);
	bool right_extend(int e);
	bool right_extend(const vector<int> &s);
	bool left_dominate(int e);
	bool right_dominate(int e);
};

#endif
