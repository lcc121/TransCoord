/*
This file was modified from a file named hyper_set.cc in Scallop
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
#include "hyper_set.h"
#include "config.h"
#include <algorithm>
#include <cstdio>

int hyper_set::clear()
{
	nodes.clear();
	edges.clear();
	e2s.clear();
	ecnts.clear();
	eucnts.clear();
	uppcnts.clear();
	return 0;
}

int hyper_set::add_node_list(const set<int> &s, int uc, int upe)
{
	return add_node_list(s, 1, uc, upe);
}

int hyper_set::add_node_list(const set<int> &s, int c, int uc, int upe)
{
	vector<int> v(s.begin(), s.end());
	return add_node_list(v, c, uc, upe);
}

int hyper_set::add_node_list(const vector<int> &s, int c, int uc, int upe)
{
	vector<int> v = s;
	sort(v.begin(), v.end());
	for (int i = 0; i < v.size(); i++) v[i]++;
	node_info ni;
	ni.count = c;
	ni.ucnt = uc;
	ni.uppcnt = upe;
	if (nodes.find(v) == nodes.end()) nodes.insert(PVINI(v, ni));
	else
	{
		nodes[v].count += c;
		nodes[v].ucnt += uc;
		nodes[v].uppcnt += upe;
	}
	return 0;
}

int hyper_set::build(directed_graph &gr, MEI& e2i)
{
	build_edges(gr, e2i);
	build_index();
	return 0;
}

int hyper_set::build_edges(directed_graph &gr, MEI& e2i)
{
	edges.clear();
	ecnts.clear();
	eucnts.clear();
	uppcnts.clear();
	for (MVINI::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		int c = it->second.count;
		int uc = it->second.ucnt;
		int uppc = it->second.uppcnt;
		if (c < min_router_count) continue;//1

		const vector<int> &vv = it->first;
		vector<int> ve;
		bool b = true;
		for (int k = 0; k < vv.size() - 1; k++)
		{
			assert(vv[k] < vv[k + 1]);
			PEB p = gr.edge(vv[k], vv[k + 1]);
			if (p.second == false) b = false;
			if (p.second == false) ve.push_back(-1);
			else ve.push_back(e2i[p.first]);
		}

		if (b == true && ve.size() >= 2)
		{
			edges.push_back(ve);
			ecnts.push_back(c);
			eucnts.push_back(uc);
			uppcnts.push_back(uppc);
		}
		continue;

		vector<int> v;
		for (int k = 0; k < ve.size(); k++)
		{
			if (ve[k] == -1)
			{
				if (v.size() >= 2)
				{
					edges.push_back(v);
					ecnts.push_back(c);
					eucnts.push_back(uc);
					uppcnts.push_back(uppc);
				}
				v.clear();
			}
			else
			{
				v.push_back(ve[k]);
			}
		}
		if (v.size() >= 2)
		{
			edges.push_back(v);
			ecnts.push_back(c);
			eucnts.push_back(uc);
			uppcnts.push_back(uppc);
		}
	}
	assert(edges.size() == ecnts.size());
	assert(edges.size() == eucnts.size());
	assert(edges.size() == uppcnts.size());
	return 0;
}

int hyper_set::build_index()
{
	e2s.clear();
	for (int i = 0; i < edges.size(); i++)
	{
		vector<int> &v = edges[i];
		for (int j = 0; j < v.size(); j++)
		{
			int e = v[j];
			if (e == -1) continue;
			if (e2s.find(e) == e2s.end())
			{
				set<int> s;
				s.insert(i);
				e2s.insert(PISI(e, s));
			}
			else
			{
				e2s[e].insert(i);
			}
		}
	}
	return 0;
}

int hyper_set::update_index()
{
	vector<int> fb1;
	for (MISI::iterator p = e2s.begin(); p != e2s.end(); p++)
	{
		int e = p->first;
		set<int> &ss = p->second;
		vector<int> fb2;
		for (set<int>::iterator it = ss.begin(); it != ss.end(); it++)
		{
			vector<int> &v = edges[*it];
			for (int i = 0; i < v.size(); i++)
			{
				if (v[i] != e) continue;
				bool b1 = false, b2 = false;
				if (i == 0 || v[i - 1] == -1) b1 = true;
				if (i == v.size() - 1 || v[i + 1] == -1) b2 = true;
				if (b1 == true && b2 == true) fb2.push_back(*it);
				break;
			}
		}
		for (int i = 0; i < fb2.size(); i++) ss.erase(fb2[i]);
		if (ss.size() == 0) fb1.push_back(e);
	}
	for (int i = 0; i < fb1.size(); i++) e2s.erase(fb1[i]);
	return 0;
}

set<int> hyper_set::get_intersection(const vector<int> &v)
{
	set<int> ss;
	if (v.size() == 0) return ss;
	assert(v[0] >= 0);
	if (e2s.find(v[0]) == e2s.end()) return ss;
	ss = e2s[v[0]];
	vector<int> vv(ss.size());
	for (int i = 1; i < v.size(); i++)
	{
		assert(v[i] >= 0);
		set<int> s;
		if (e2s.find(v[i]) == e2s.end()) return s;
		s = e2s[v[i]];
		vector<int>::iterator it = set_intersection(ss.begin(), ss.end(), s.begin(), s.end(), vv.begin());//algorithm, 求两个有序集合的交集存入vv，返回最后一个元素下一个
		ss = set<int>(vv.begin(), it);
	}
	return ss;
}

MI hyper_set::get_successors(int e)
{
	MI s;
	if (e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for (set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		int c = ecnts[*it];
		for (int i = 0; i < v.size(); i++)
		{
			if (v[i] != e) continue;
			if (i >= v.size() - 1) continue;
			int k = v[i + 1];
			if (k == -1) continue;
			if (s.find(k) == s.end()) s.insert(PI(k, c));
			else s[k] += c;
		}
	}
	return s;
}
MI hyper_set::get_predecessors(int e)
{
	MI s;
	if (e2s.find(e) == e2s.end()) return s;
	set<int> &ss = e2s[e];
	for (set<int>::iterator it = ss.begin(); it != ss.end(); it++)
	{
		vector<int> &v = edges[*it];
		int c = ecnts[*it];
		for (int i = 0; i < v.size(); i++)
		{
			if (v[i] != e) continue;
			if (i == 0) continue;
			int k = v[i - 1];
			if (k == -1) continue;
			if (s.find(k) == s.end()) s.insert(PI(k, c));
			else s[k] += c;
		}
	}
	return s;
}

int hyper_set::get_heaviest_successor(const set<int> &s, int e)
{
	int ei = -1;
	int max = 0;
	MI sc = get_successors(e);
	MI::iterator it;
	for (it = sc.begin(); it != sc.end(); ++it)
	{
		if (s.find(it->first) == s.end()) continue;
		if (it->second <= max) continue;
		max = it->second;
		ei = it->first;
	}
	if (verbose >= 6)
	{
		cout << "find hs successors: ";
		for (it = sc.begin(); it != sc.end(); ++it) cout << it->first << "(counts=" << it->second << ") ";
		cout << "\nchoose " << ei << endl;
	}
	return ei;
}
int hyper_set::get_heaviest_predecessor(const set<int> &s, int e)
{
	int ei = -1;
	int max = 0;
	MI pc = get_predecessors(e);
	MI::iterator it;
	for (it = pc.begin(); it != pc.end(); ++it)
	{
		if (s.find(it->first) == s.end()) continue;
		if (it->second <= max) continue;
		max = it->second;
		ei = it->first;
	}
	if (verbose >= 6)
	{
		cout << "find hs predecessors: ";
		for (it = pc.begin(); it != pc.end(); ++it) cout << it->first << "(counts=" << it->second << ") ";
		cout << "\nchoose " << ei << endl;
	}
	return ei;
}

MI hyper_set::get_compatible_predecessors(vector<int> &e)
{
	MI m;
	if (e.empty()) return m;
	if (e2s.find(e[0]) == e2s.end()) return m;
	set<int> &ss = e2s[e[0]];
	for (set<int>::iterator it = ss.begin(); it != ss.end(); ++it)
	{
		vector<int> &v = edges[*it];
		int i = 0, j = 0;
		for (; i < v.size(); i++) { if (v[i] == e[0]) break; }
		assert(i < v.size());
		if (i == 0) continue;
		int l = v[i - 1];
		if (l == -1) continue;
		while (i < v.size() && j < e.size())
		{
			if (v[i] == e[j]) { i += 1; j += 1; }
			else break;
		}
		if (i < v.size() && j < e.size()) continue;
		int c = ecnts[*it];
		if (m.find(l) == m.end()) m.insert(PI(l, c));
		else m[l] += c;
	}
	return m;
}
MI hyper_set::get_compatible_successors(vector<int> &e)
{
	MI m;
	if (e.empty()) return m;
	if (e2s.find(e.back()) == e2s.end()) return m;
	set<int> &ss = e2s[e.back()];
	for (set<int>::iterator it = ss.begin(); it != ss.end(); ++it)
	{
		vector<int> &v = edges[*it];
		int i = 0, j = (int)e.size() - 1;
		for (; i < v.size(); i++) { if (v[i] == e.back()) break; }
		assert(i < v.size());
		if (i == (int)v.size() - 1) continue;
		int r = v[i + 1];
		while (i >= 0 && j >= 0)
		{
			if (v[i] == e[j]) { i -= 1; j -= 1; }
			else break;
		}
		if (i >= 0 && j >= 0) continue;
		int c = ecnts[*it];
		if (m.find(r) == m.end()) m.insert(PI(r, c));
		else m[r] += c;
	}
	return m;
}
int hyper_set::get_heaviest_compatible_predecessor(const set<int> &s, vector<int> &e)
{
	MI pc = get_compatible_predecessors(e);
	return get_heaviest_exist_compatible(pc, s, e);
}
int hyper_set::get_heaviest_compatible_successor(const set<int> &s, vector<int> &e)
{
	MI sc = get_compatible_successors(e);
	return get_heaviest_exist_compatible(sc, s, e);
}
int hyper_set::get_heaviest_exist_compatible(MI &m, const set<int> &s, vector<int> &e)
{
	int ei = -1;
	int max = 0;
	MI::iterator it;
	for (it = m.begin(); it != m.end(); ++it)
	{
		if (s.find(it->first) == s.end()) continue;
		if (it->second <= max) continue;
		max = it->second;
		ei = it->first;
	}
	if (verbose >= 6)
	{
		cout << "find compatible: ";
		for (it = m.begin(); it != m.end(); ++it) { cout << it->first << "(counts=" << it->second << ") "; }
		cout << "\nchoose " << ei << endl;
	}
	return ei;
}

MPII hyper_set::get_routes(int x, directed_graph &gr, MEI &e2i)
{
	MPII mpi;
	edge_iterator it1, it2;
	PEEI pei;
	vector<PI> v;
	for (pei = gr.in_edges(x), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		assert(e2i.find(*it1) != e2i.end());
		int e = e2i[*it1];
		MI s = get_successors(e);
		for (MI::iterator it = s.begin(); it != s.end(); it++)
		{
			PI p(e, it->first);
			mpi.insert(PPII(p, it->second));
		}
	}
	return mpi;
}

MEI hyper_set::get_single_edge_counts(splice_graph &gr, VE &i2e)
{
	MEI secnts;
	MED::iterator it;
	for (it = gr.ewrt.begin(); it != gr.ewrt.end(); ++it) secnts.insert(PEI(it->first, (int)it->second));
	for (int i = 0; i < edges.size(); i++)
	{
		vector<int> &v = edges[i];
		for (int j = 0; j < v.size(); j++)
		{
			edge_descriptor e = i2e[v[j]];
			assert(secnts.find(e) != secnts.end());
			if (secnts[e] > ecnts[i]) secnts[e] -= ecnts[i];
			else secnts[e] = 0;
		}
	}

	return secnts;
}


int hyper_set::replace(int x, int e)
{
	vector<int> v;
	v.push_back(x);
	replace(v, e);
	return 0;
}

int hyper_set::replace(int x, int y, int e)
{
	vector<int> v;
	v.push_back(x);
	v.push_back(y);
	replace(v, e);
	return 0;
}

int hyper_set::replace(const vector<int> &v, int e)
{
	if (v.size() == 0) return 0;
	set<int> s = get_intersection(v);

	vector<int> fb;
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		vector<int> bv = consecutive_subset(vv, v);

		if (bv.size() <= 0) continue;
		assert(bv.size() == 1);

		int b = bv[0];
		vv[b] = e;

		bool b1 = useful(vv, 0, b);
		bool b2 = useful(vv, b + v.size() - 1, vv.size() - 1);

		if (b1 == false && b2 == false)
		{
			fb.push_back(k);
			continue;
		}

		vv.erase(vv.begin() + b + 1, vv.begin() + b + v.size());

		if (e2s.find(e) == e2s.end())
		{
			set<int> ss;
			ss.insert(k);
			e2s.insert(PISI(e, ss));
		}
		else
		{
			e2s[e].insert(k);
		}
	}

	for (int i = 0; i < v.size(); i++)
	{
		int u = v[i];
		if (e2s.find(u) == e2s.end()) continue;
		for (int k = 0; k < fb.size(); k++) e2s[u].erase(fb[k]);
		if (e2s[u].size() == 0) e2s.erase(u);
	}
	return 0;
}

int hyper_set::remove(const set<int> &s)
{
	return remove(vector<int>(s.begin(), s.end()));
}

int hyper_set::remove(const vector<int> &v)
{
	for (int i = 0; i < v.size(); i++) remove(v[i]);
	return 0;
}

int hyper_set::remove(int e)
{
	if (e2s.find(e) == e2s.end()) return 0;
	set<int> &s = e2s[e];
	vector<int> fb;
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 0; i < vv.size(); i++)
		{
			if (vv[i] != e) continue;

			vv[i] = -1;

			bool b1 = useful(vv, 0, i - 1);
			bool b2 = useful(vv, i + 1, vv.size() - 1);
			if (b1 == false && b2 == false) fb.push_back(k);

			break;
		}
	}

	for (int i = 0; i < fb.size(); i++) s.erase(fb[i]);
	e2s.erase(e);
	return 0;
}

int hyper_set::remove_pair(int x, int y)
{
	if (e2s.find(x) == e2s.end()) return 0;
	set<int> &s = e2s[x];
	vector<int> fb;
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 0; i < vv.size(); i++)
		{
			if (i == vv.size() - 1) continue;
			if (vv[i] != x) continue;
			if (vv[i + 1] != y) continue;

			bool b1 = useful(vv, 0, i);
			bool b2 = (b1 == true) ? true : useful(vv, i + 1, vv.size() - 1);

			if (b1 == false && b2 == false) fb.push_back(k);
			else vv.insert(vv.begin() + i + 1, -1);

			break;
		}
	}

	for (int i = 0; i < fb.size(); i++) s.erase(fb[i]);
	if (s.size() == 0) e2s.erase(x);

	return 0;
}

bool hyper_set::useful(const vector<int> &v, int k1, int k2)
{
	for (int i = k1; i < k2; i++)
	{
		if (v[i] >= 0 && v[i + 1] >= 0) return true;
	}
	return false;
}

int hyper_set::insert_between(int x, int y, int e)
{
	if (e2s.find(x) == e2s.end()) return 0;
	set<int> s = e2s[x];
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 0; i < vv.size(); i++)
		{
			if (i == vv.size() - 1) continue;
			if (vv[i] != x) continue;
			if (vv[i + 1] != y) continue;
			vv.insert(vv.begin() + i + 1, e);

			if (e2s.find(e) == e2s.end())
			{
				set<int> ss;
				ss.insert(k);
				e2s.insert(PISI(e, ss));
			}
			else
			{
				e2s[e].insert(k);
			}

			//printf("line %d: insert %d between (%d, %d) = (%d, %d, %d)\n", k, e, x, y, vv[i], vv[i + 1], vv[i + 2]);

			break;
		}
	}
	return 0;
}

bool hyper_set::extend(int e)
{
	return (left_extend(e) || right_extend(e));
}

bool hyper_set::left_extend(int e)
{
	if (e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 1; i < vv.size(); i++)
		{
			if (vv[i] == e && vv[i - 1] != -1) return true;
		}
	}
	return false;
}

bool hyper_set::right_extend(int e)
{
	if (e2s.find(e) == e2s.end()) return false;
	set<int> s = e2s[e];
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 0; i < vv.size() - 1; i++)
		{
			if (vv[i] == e && vv[i + 1] != -1) return true;
		}
	}
	return false;
}

bool hyper_set::left_extend(const vector<int> &s)
{
	for (int i = 0; i < s.size(); i++)
	{
		if (left_extend(s[i]) == true) return true;
	}
	return false;
}

bool hyper_set::right_extend(const vector<int> &s)
{
	for (int i = 0; i < s.size(); i++)
	{
		if (right_extend(s[i]) == true) return true;
	}
	return false;
}

bool hyper_set::left_dominate(int e)
{
	// for each appearance of e
	// if right is not empty then left is also not empty
	if (e2s.find(e) == e2s.end()) return true;

	set<PI> x1;
	set<PI> x2;
	set<int> s = e2s[e];
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);

		for (int i = 0; i < vv.size() - 1; i++)
		{
			if (vv[i] != e) continue;
			if (vv[i + 1] == -1) continue;

			if (i == 0 || vv[i - 1] == -1)
			{
				if (i + 2 < vv.size()) x1.insert(PI(vv[i + 1], vv[i + 2]));
				else x1.insert(PI(vv[i + 1], -1));
			}
			else
			{
				x2.insert(PI(vv[i + 1], -1));
				if (i + 2 < vv.size()) x2.insert(PI(vv[i + 1], vv[i + 2]));
			}
		}
	}

	for (set<PI>::iterator it = x1.begin(); it != x1.end(); it++)
	{
		PI p = (*it);
		if (x2.find(p) == x2.end()) return false;
	}
	return true;
}

bool hyper_set::right_dominate(int e)
{
	// for each appearance of e
	// if left is not empty then right is also not empty
	if (e2s.find(e) == e2s.end()) return true;
	set<PI> x1;
	set<PI> x2;
	set<int> s = e2s[e];
	for (set<int>::iterator it = s.begin(); it != s.end(); it++)
	{
		int k = (*it);
		vector<int> &vv = edges[k];
		assert(vv.size() >= 1);
		for (int i = 1; i < vv.size(); i++)
		{
			if (vv[i] != e) continue;
			if (vv[i - 1] == -1) continue;

			if (i == vv.size() - 1 || vv[i + 1] == -1)
			{
				if (i - 2 >= 0) x1.insert(PI(vv[i - 1], vv[i - 2]));
				else x1.insert(PI(vv[i - 1], -1));
			}
			else
			{
				x2.insert(PI(vv[i - 1], -1));
				if (i - 2 >= 0) x2.insert(PI(vv[i - 1], vv[i - 2]));
			}
		}
	}

	for (set<PI>::iterator it = x1.begin(); it != x1.end(); it++)
	{
		PI p = (*it);
		if (x2.find(p) == x2.end()) return false;
	}

	return true;
}

int hyper_set::print()
{
	//printf("PRINT HYPER_SET\n");
	for (MVINI::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second.count;
		int uc = it->second.ucnt;
		int uppc = it->second.uppcnt;
		printf("hyper-edge (nodes), count = %d, ucnt = %d, uppcnt = %d, list = ( ", c, uc, uppc);
		writev(cout, v);
		printf(")\n");
	}

	/*
	for(int i = 0; i < edges.size(); i++)
	{
	printf("hyper-edge (edges) %d: ( ", i);
	printv(edges[i]);
	printf(")\n");
	}
	*/
	return 0;
}

int hyper_set::write(ostream &out)
{
	for (MVINI::iterator it = nodes.begin(); it != nodes.end(); it++)
	{
		const vector<int> &v = it->first;
		int c = it->second.count;
		int uc = it->second.ucnt;
		int uppc = it->second.uppcnt;
		out << "hyper-set(nodes), count=" << c << ", ucnt=" << uc << ", uppcnt=" << uppc << ", list={";
		writev(out, v);
		out << "}\n";
	}
	return 0;
}

int hyper_set::write_edges(ostream &out)
{
	for (int i = 0; i < edges.size(); i++)
	{
		out << "hyper-set(edges[" << i << "]), counts=" << ecnts[i] << ", unique_counts=" << eucnts[i] << ", uppcnts=" << uppcnts[i] << ", list{";
		writev(out, edges[i]);
		out << "}\n";
	}
	return 0;
}
