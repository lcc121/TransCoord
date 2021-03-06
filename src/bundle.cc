/*
This file was modified from a file named bundle.cc in Scallop
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
#include <cassert>
#include <cstdio>
#include <map>
#include <iomanip>
#include <fstream>

#include "bundle.h"
#include "region.h"
#include "config.h"
#include "util.h"
#include "undirected_graph.h"

bundle::bundle()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';
}

bundle::~bundle()
{}

int bundle::add_hit(const hit &ht)
{
	// store new hit
	hits.push_back(ht);

	// calcuate the boundaries on reference
	if (ht.pos < lpos) lpos = ht.pos;
	if (ht.rpos > rpos) rpos = ht.rpos;

	// set tid
	if (tid == -1) tid = ht.tid;
	assert(tid == ht.tid);

	// set strand
	if (hits.size() <= 1) strand = ht.strand;
	assert(strand == ht.strand);

	// DEBUG
	/*
	if(strand != ht.strand)
	{
	printf("strand = %c, ht.strand = %c, ht.xs = %c,\n", strand, ht.strand, ht.xs);
	}
	*/

	for (int k = 0; k < ht.itvm.size(); k++)
	{
		int32_t s = high32(ht.itvm[k]);
		int32_t t = low32(ht.itvm[k]);
		//printf(" add interval %d-%d\n", s, t);
		mmap += make_pair(ROI(s, t), 1);
	}

	for (int k = 0; k < ht.itvi.size(); k++)
	{
		int32_t s = high32(ht.itvi[k]);
		int32_t t = low32(ht.itvi[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	for (int k = 0; k < ht.itvd.size(); k++)
	{
		int32_t s = high32(ht.itvd[k]);
		int32_t t = low32(ht.itvd[k]);
		imap += make_pair(ROI(s, t), 1);
	}

	return 0;
}

int bundle::build()
{
	compute_strand();        

	check_left_ascending();  

	build_junctions();       
							 

	build_regions();         
	build_partial_exons();   

	build_partial_exon_map();
	link_partial_exons();    
	build_splice_graph();    

	revise_splice_graph();   

	build_hyper_edges2();    

	return 0;
}

int bundle::clear()
{
	tid = -1;
	chrm = "";
	lpos = 1 << 30;
	rpos = 0;
	strand = '.';

	hits.clear();
	mmap.clear();
	imap.clear();

	junctions.clear();
	regions.clear();
	pexons.clear();
	regional.clear();
	pmap.clear();
	gr.clear();
	hs.clear();
	return 0;
}

int bundle::clear_memory()
{
	// left only gr&hs
	hits.clear();
	hits.shrink_to_fit();
	mmap.clear();
	imap.clear();

	junctions.clear();
	junctions.shrink_to_fit();
	regions.clear();
	regions.shrink_to_fit();
	pexons.clear();
	pexons.shrink_to_fit();
	regional.clear();
	regional.shrink_to_fit();
	pmap.clear();
	return 0;
}


int bundle::compute_strand()
{
	if (library_type != UNSTRANDED) assert(strand != '.');
	if (library_type != UNSTRANDED) return 0;

	int n0 = 0, np = 0, nq = 0;
	for (int i = 0; i < hits.size(); i++)
	{
		if (hits[i].xs == '.') n0++;
		if (hits[i].xs == '+') np++;
		if (hits[i].xs == '-') nq++;
	}

	if (np > nq) strand = '+';
	else if (np < nq) strand = '-';
	else strand = '.';

	return 0;
}

int bundle::check_left_ascending()
{
	for (int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].pos;
		int32_t p2 = hits[i].pos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::check_right_ascending()
{
	for (int i = 1; i < hits.size(); i++)
	{
		int32_t p1 = hits[i - 1].rpos;
		int32_t p2 = hits[i].rpos;
		assert(p1 <= p2);
	}
	return 0;
}

int bundle::build_junctions()
{
	int min_max_boundary_quality = min_mapping_quality;
	map< int64_t, vector<int> > m;
	for (int i = 0; i < hits.size(); i++)
	{
		vector<int64_t> v = hits[i].spos;
		if (v.size() == 0) continue;

		//hits[i].print();
		for (int k = 0; k < v.size(); k++)
		{
			int64_t p = v[k];

			// DEBUG
			/*
			int32_t x1 = low32(p);
			int32_t x2 = high32(p);
			if(fabs(x1 - 9364768) <= 2 || fabs(x2 - 9364768) <=2)
			{
			printf("HIT ");
			hits[i].print();
			}
			*/
			if (m.find(p) == m.end())
			{
				vector<int> hv;
				hv.push_back(i);
				m.insert(pair< int64_t, vector<int> >(p, hv));
			}
			else
			{
				m[p].push_back(i);
			}
		}
	}

	map< int64_t, vector<int> >::iterator it;
	junctions.reserve(m.size());
	for (it = m.begin(); it != m.end(); it++)
	{
		vector<int> &v = it->second;
		if (v.size() < min_splice_boundary_hits) continue;

		int32_t p1 = high32(it->first);
		int32_t p2 = low32(it->first);

		int s0 = 0;
		int s1 = 0;
		int s2 = 0;
		int nm = 0;
		for (int k = 0; k < v.size(); k++)
		{
			hit &h = hits[v[k]];
			nm += h.nm;
			if (h.xs == '.') s0++;
			if (h.xs == '+') s1++;
			if (h.xs == '-') s2++;
		}

		//printf("junction: %s:%d-%d (%d, %d, %d) %d\n", chrm.c_str(), p1, p2, s0, s1, s2, s1 < s2 ? s1 : s2);

		junction jc(it->first, v.size());
		jc.nm = nm;
		if (s1 == 0 && s2 == 0) jc.strand = '.';
		else if (s1 >= 1 && s2 >= 1) jc.strand = '.';
		else if (s1 > s2) jc.strand = '+';
		else jc.strand = '-';
		junctions.push_back(jc);

		/*
		uint32_t max_qual = 0;
		for(int k = 0; k < v.size(); k++)
		{
		hit &h = hits[v[k]];
		if(h.qual > max_qual) max_qual = h.qual;
		}
		assert(max_qual >= min_max_boundary_quality);
		*/
	}
	junctions.shrink_to_fit();
	return 0;
}

int bundle::correct_junctions()
{
	if (junctions.size() == 0) return 0;
	sort(junctions.begin(), junctions.end(), junction_cmp_length);
	set<int> fb;
	for (int k = 1; k < junctions.size(); k++)
	{
		if (fb.find(k - 1) != fb.end()) continue;
		if (fb.find(k - 0) != fb.end()) continue;

		junction &j1 = junctions[k - 1];
		junction &j2 = junctions[k - 0];

		if (fabs(j1.lpos - j2.lpos) >= 10) continue;
		if (fabs(j1.rpos - j2.rpos) >= 10) continue;
		if (fabs((j1.rpos - j1.lpos) - (j2.rpos - j2.lpos)) >= 10) continue;

		double nm1 = j1.nm * 1.0 / j1.count;
		double nm2 = j2.nm * 1.0 / j2.count;
		if (nm1 < nm2 - 0.8)
		{
			// correct nm2 to nm1
			fb.insert(k);

			if (j1.lpos < j2.lpos) mmap += make_pair(ROI(j1.lpos + 1, j2.lpos + 1), -1);
			else if (j1.lpos > j2.lpos) mmap += make_pair(ROI(j2.lpos + 1, j1.lpos + 1), 1);

			if (j1.rpos < j2.rpos) mmap += make_pair(ROI(j1.rpos, j2.rpos), 1);
			else if (j1.rpos > j2.rpos) mmap += make_pair(ROI(j2.rpos, j1.rpos), -1);

			if (verbose >= 2)
			{
				j1.print(chrm, k - 1);
				j2.print(chrm, k - 0);
			}
		}
		else if (nm2 < nm1 - 0.8)
		{
			// correct nm1 to nm2
			fb.insert(k - 1);
			if (j2.lpos < j1.lpos) mmap += make_pair(ROI(j2.lpos + 1, j1.lpos + 1), -1);
			else if (j2.lpos > j1.lpos) mmap += make_pair(ROI(j1.lpos + 1, j2.lpos + 1), 1);

			if (j2.rpos < j1.rpos) mmap += make_pair(ROI(j2.rpos, j1.rpos), 1);
			else if (j2.rpos > j1.rpos) mmap += make_pair(ROI(j1.rpos, j2.rpos), -1);

			if (verbose >= 2)
			{
				j1.print(chrm, k - 1);
				j2.print(chrm, k - 0);
			}
		}
	}

	vector<junction> v;
	for (int i = 0; i < junctions.size(); i++)
	{
		if (fb.find(i) != fb.end()) continue;
		v.push_back(junctions[i]);
	}
	junctions = v;
	return 0;
}

int bundle::build_regions()
{
	MPI s;
	s.insert(PI(lpos, START_BOUNDARY));
	s.insert(PI(rpos, END_BOUNDARY));  
	for (int i = 0; i < junctions.size(); i++)
	{
		junction &jc = junctions[i];

		double ave, dev;
		evaluate_rectangle(mmap, jc.lpos, jc.rpos, ave, dev);

		int32_t l = jc.lpos;
		int32_t r = jc.rpos;

		if (s.find(l) == s.end()) s.insert(PI(l, LEFT_SPLICE));
		else if (s[l] == RIGHT_SPLICE) s[l] = LEFT_RIGHT_SPLICE;

		if (s.find(r) == s.end()) s.insert(PI(r, RIGHT_SPLICE));
		else if (s[r] == LEFT_SPLICE) s[r] = LEFT_RIGHT_SPLICE;
	}

	for (int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		if (s.find(p.lpos) != s.end()) s.insert(PI(p.lpos, p.ltype));
		if (s.find(p.rpos) != s.end()) s.insert(PI(p.rpos, p.rtype));
	}

	vector<PPI> v(s.begin(), s.end());
	sort(v.begin(), v.end());

	regions.clear();
	regions.reserve((int)v.size() - 1);
	for (int k = 0; k < (int)v.size() - 1; k++)
	{
		int32_t l = v[k].first;
		int32_t r = v[k + 1].first;
		int ltype = v[k].second;
		int rtype = v[k + 1].second;

		if (ltype == LEFT_RIGHT_SPLICE) ltype = RIGHT_SPLICE;
		if (rtype == LEFT_RIGHT_SPLICE) rtype = LEFT_SPLICE;

		regions.push_back(region(l, r, ltype, rtype, &mmap, &imap));
	}

	return 0;
}

int bundle::build_partial_exons()
{
	pexons.clear();
	regional.clear();
	for (int i = 0; i < regions.size(); i++)
	{
		region &r = regions[i];
		for (int k = 0; k < r.pexons.size(); k++)
		{
			partial_exon &pe = r.pexons[k];
			pexons.push_back(pe);
			if ((pe.lpos != lpos || pe.rpos != rpos) && pe.ltype == START_BOUNDARY && pe.rtype == END_BOUNDARY) regional.push_back(true);
			else regional.push_back(false);
			//printf("regional = %s, ", regional.back() ? "TRUE" : "FALSE"); pe.print(k);
		}
		//pexons.insert(pexons.end(), r.pexons.begin(), r.pexons.end());
	}
	return 0;
}

int bundle::build_partial_exon_map()
{
	pmap.clear();
	for (int i = 0; i < pexons.size(); i++)
	{
		partial_exon &p = pexons[i];
		pmap += make_pair(ROI(p.lpos, p.rpos), i + 1);
	}
	return 0;
}

int bundle::locate_left_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x, x + 1));
	if (it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p2 >= x);
	assert(p1 <= x);

	if (x - p1 > min_flank_length && p2 - x < min_flank_length) k++;

	if (k >= pexons.size()) return -1;
	return k;
}

int bundle::locate_right_partial_exon(int32_t x)
{
	SIMI it = pmap.find(ROI(x - 1, x));
	if (it == pmap.end()) return -1;
	assert(it->second >= 1);
	assert(it->second <= pexons.size());

	int k = it->second - 1;
	int32_t p1 = lower(it->first);
	int32_t p2 = upper(it->first);
	assert(p1 < x);
	assert(p2 >= x);

	if (p2 - x > min_flank_length && x - p1 <= min_flank_length) k--;
	return k;
}

int bundle::build_hyper_edges2()
{
	//sort(hits.begin(), hits.end(), hit_compare_by_name);
	sort(hits.begin(), hits.end());

	/*
	printf("----------------------\n");
	for(int k = 9; k < hits.size(); k++) hits[k].print();
	printf("======================\n");
	*/

	hs.clear();

	string qname;
	int hi = -2;
	bool pp = false;
	vector<int> sp1;
					
	for (int i = 0; i < hits.size(); i++)
	{
		hit &h = hits[i];

		/*
		printf("sp1 = ( ");
		printv(sp1);
		printf(")\n");
		h.print();*/

		if (h.qname != qname || h.hi != hi)
		{
			set<int> s(sp1.begin(), sp1.end());
			if (s.size() >= 2)
			{
				int ucount = (h.nh == 1) ? 1 : 0;
				int upe = pp ? 1 : 0;
				hs.add_node_list(s, ucount, upe);
			}
			sp1.clear();
			pp = false;
		}

		qname = h.qname;
		hi = h.hi;

		if ((h.flag & 0x4) >= 1) continue;

		vector<int> sp2;
		for (int k = 0; k < h.itvm.size(); k++)
		{
			int32_t p1 = high32(h.itvm[k]);
			int32_t p2 = low32(h.itvm[k]);

			int k1 = locate_left_partial_exon(p1);
			int k2 = locate_right_partial_exon(p2);
			if (k1 < 0 || k2 < 0) continue;

			for (int j = k1; j <= k2; j++) sp2.push_back(j);
		}

		if (sp1.size() <= 0 || sp2.size() <= 0)
		{
			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
			continue;
		}

		/*
		printf("sp2 = ( ");
		printv(sp2);
		printf(")\n");
		*/
		//paired-end 
		int x1 = -1, x2 = -1;
		if (h.isize < 0)  
		{
			x1 = sp1[max_element(sp1)];
			x2 = sp2[min_element(sp2)];
		}
		else
		{
			x1 = sp2[max_element(sp2)];
			x2 = sp1[min_element(sp1)];
		}

		vector<int> sp3;
		bool c = bridge_read(x1, x2, sp3);
		if (c == false)
		{
			set<int> s(sp1.begin(), sp1.end());
			if (s.size() >= 2)
			{
				int ucount = (h.nh == 1) ? 1 : 0;
				hs.add_node_list(s, ucount, 0);
			}
			sp1 = sp2;
		}
		else
		{
			sp1.insert(sp1.end(), sp2.begin(), sp2.end());
			sp1.insert(sp1.end(), sp3.begin(), sp3.end());
			if (sp1.size() >= 3) pp = true;
		}
	}

	return 0;
}

bool bundle::bridge_read(int x, int y, vector<int> &v)
{
	v.clear();
	if (x >= y) return true;

	PEB e = gr.edge(x + 1, y + 1);
	if (e.second == true) return true;
	//else return false;

	if (y - x >= 6) return false;

	long max = 9999999999;
	vector<long> table;
	vector<int> trace;
	int n = y - x + 1;
	table.resize(n, 0);
	trace.resize(n, -1);
	table[0] = 1;
	trace[0] = -1;
	for (int i = x + 1; i <= y; i++)
	{
		edge_iterator it1, it2;
		PEEI pei;
		for (pei = gr.in_edges(i + 1), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			int s = (*it1)->source() - 1;
			int t = (*it1)->target() - 1;
			assert(t == i);
			if (s < x) continue;
			if (table[s - x] <= 0) continue;
			table[t - x] += table[s - x];
			trace[t - x] = s - x;
			if (table[t - x] >= max) return false;
		}
	}

	//printf("x = %d, y = %d, num-paths = %ld\n", x, y, table[n - 1]);
	if (table[n - 1] != 1) return false;

	//printf("path = ");

	v.clear();
	int p = n - 1;
	while (p >= 0)
	{
		p = trace[p];
		if (p <= 0) break;
		v.push_back(p + x);
		//printf("%d ", p + x);
	}
	//printf("\n");
	//assert(v.size() >= 1);

	return true;
}

int bundle::link_partial_exons()
{
	if (pexons.size() == 0) return 0;

	MPI lm;
	MPI rm;
	for (int i = 0; i < pexons.size(); i++)
	{
		int32_t l = pexons[i].lpos;
		int32_t r = pexons[i].rpos;

		assert(lm.find(l) == lm.end());
		assert(rm.find(r) == rm.end());
		lm.insert(PPI(l, i));
		rm.insert(PPI(r, i));
	}

	for (int i = 0; i < junctions.size(); i++)
	{
		junction &b = junctions[i];
		MPI::iterator li = rm.find(b.lpos);
		MPI::iterator ri = lm.find(b.rpos);

		assert(li != rm.end());
		assert(ri != lm.end());

		if (li != rm.end() && ri != lm.end())
		{
			b.lexon = li->second;
			b.rexon = ri->second;
		}
		else
		{
			b.lexon = b.rexon = -1;
		}
	}
	return 0;
}

int bundle::build_splice_graph()
{
	gr.clear();

	// vertices: start, each region, end
	gr.add_vertex();
	vertex_info vi0;
	vi0.lpos = lpos;
	vi0.rpos = lpos;
	gr.set_vertex_weight(0, 0);
	gr.set_vertex_info(0, vi0);
	for (int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];
		int length = r.rpos - r.lpos;
		assert(length >= 1);
		gr.add_vertex();
		gr.set_vertex_weight(i + 1, r.ave < 1.0 ? 1.0 : r.ave);
		vertex_info vi;
		vi.lpos = r.lpos;
		vi.rpos = r.rpos;
		vi.length = length;
		vi.stddev = r.dev;// < 1.0 ? 1.0 : r.dev;
		vi.regional = regional[i];
		gr.set_vertex_info(i + 1, vi);
	}

	gr.add_vertex();
	vertex_info vin;
	vin.lpos = rpos;
	vin.rpos = rpos;
	gr.set_vertex_weight(pexons.size() + 1, 0);
	gr.set_vertex_info(pexons.size() + 1, vin);

	// edges: each junction => and e2w
	for (int i = 0; i < junctions.size(); i++)
	{
		const junction &b = junctions[i];

		if (b.lexon < 0 || b.rexon < 0) continue;

		//const partial_exon &x = pexons[b.lexon];
		//const partial_exon &y = pexons[b.rexon];
		edge_descriptor p = gr.add_edge(b.lexon + 1, b.rexon + 1);
		assert(b.count >= 1);
		edge_info ei;
		ei.weight = b.count;
		ei.strand = b.strand;
		gr.set_edge_info(p, ei);
		gr.set_edge_weight(p, b.count);
	}

	// edges: connecting start/end and pexons
	int ss = 0;
	int tt = pexons.size() + 1;
	for (int i = 0; i < pexons.size(); i++)
	{
		const partial_exon &r = pexons[i];

		if (r.ltype == START_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(ss, i + 1);
			double w = r.ave;
			if (i >= 1 && pexons[i - 1].rpos == r.lpos) w -= pexons[i - 1].ave;
			if (w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}

		if (r.rtype == END_BOUNDARY)
		{
			edge_descriptor p = gr.add_edge(i + 1, tt);
			double w = r.ave;
			if (i < pexons.size() - 1 && pexons[i + 1].lpos == r.rpos) w -= pexons[i + 1].ave;
			if (w < 1.0) w = 1.0;
			gr.set_edge_weight(p, w);
			edge_info ei;
			ei.weight = w;
			gr.set_edge_info(p, ei);
		}
	}

	// edges: connecting adjacent pexons => e2w
	for (int i = 0; i < (int)(pexons.size()) - 1; i++)
	{
		const partial_exon &x = pexons[i];
		const partial_exon &y = pexons[i + 1];

		if (x.rpos != y.lpos) continue;

		assert(x.rpos == y.lpos);

		int xd = gr.out_degree(i + 1);
		int yd = gr.in_degree(i + 2);
		double wt = (xd < yd) ? x.ave : y.ave;
		//int32_t xr = compute_overlap(mmap, x.rpos - 1);
		//int32_t yl = compute_overlap(mmap, y.lpos);
		//double wt = xr < yl ? xr : yl;

		edge_descriptor p = gr.add_edge(i + 1, i + 2);
		double w = (wt < 1.0) ? 1.0 : wt;
		gr.set_edge_weight(p, w);
		edge_info ei;
		ei.weight = w;
		gr.set_edge_info(p, ei);
	}

	gr.strand = strand;
	gr.chrm = chrm;
	return 0;
}

int bundle::revise_splice_graph()
{
	while (true)
	{
		bool b = false;

		b = extend_boundaries();
		if (b == true) continue;

		b = remove_inner_boundaries();
		if (b == true) continue;

		b = remove_small_exons();
		if (b == true) refine_splice_graph();
		if (b == true) continue;

		b = remove_small_junctions();
		if (b == true) refine_splice_graph();
		if (b == true) continue;

		b = keep_surviving_edges();
		if (b == true) refine_splice_graph();
		if (b == true) continue;

		b = remove_intron_contamination();
		if (b == true) continue;

		break;
	}

	refine_splice_graph();

	return 0;
}

bool bundle::extend_boundaries()
{
	edge_iterator it1, it2;
	PEEI pei;
	for (pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		int s = e->source();
		int t = e->target();
		int32_t p = gr.get_vertex_info(t).lpos - gr.get_vertex_info(s).rpos;
		double we = gr.get_edge_weight(e);
		double ws = gr.get_vertex_weight(s);
		double wt = gr.get_vertex_weight(t);

		if (p <= 0) continue;
		if (s == 0) continue;
		if (t == gr.num_vertices() - 1) continue;

		bool b = false;
		if (gr.out_degree(s) == 1 && ws >= 10.0 * we * we + 10.0) b = true;
		if (gr.in_degree(t) == 1 && wt >= 10.0 * we * we + 10.0) b = true;

		if (b == false) continue;

		if (gr.out_degree(s) == 1)
		{
			edge_descriptor ee = gr.add_edge(s, gr.num_vertices() - 1);
			gr.set_edge_weight(ee, ws);
			gr.set_edge_info(ee, edge_info());
		}
		if (gr.in_degree(t) == 1)
		{
			edge_descriptor ee = gr.add_edge(0, t);
			gr.set_edge_weight(ee, wt);
			gr.set_edge_info(ee, edge_info());
		}
		if (fixed_bundle_idx == -1 && verbose >= 4) printf("extend_boundaries: remove edge (%d, %d), weight = %.2lf\n", e->source(), e->target(), gr.get_edge_weight(e));
		gr.remove_edge(e);

		return true;
	}

	return false;
}

bool bundle::remove_inner_boundaries()
{
	// s->i->t
	bool flag = false;
	int n = gr.num_vertices() - 1;
	for (int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if (gr.in_degree(i) != 1) continue;
		if (gr.out_degree(i) != 1) continue;

		PEEI pei = gr.in_edges(i);
		edge_iterator it1 = pei.first, it2 = pei.second;
		edge_descriptor e1 = (*it1);

		pei = gr.out_edges(i);
		it1 = pei.first;
		it2 = pei.second;
		edge_descriptor e2 = (*it1);
		vertex_info vi = gr.get_vertex_info(i);
		int s = e1->source();
		int t = e2->target();

		if (s != 0 && t != n) continue;
		if (s != 0 && gr.out_degree(s) == 1) continue;
		if (t != n && gr.in_degree(t) == 1) continue;

		if (vi.stddev >= 0.01) continue;

		if (verbose >= 4 && fixed_bundle_idx == -1)
		{
			printf("remove_inner_boundaries: ");
			write_clear_vertex(cout, i);
		}
		gr.clear_vertex(i);
		flag = true;
	}
	return flag;
}

bool bundle::remove_small_exons()
{
	bool flag = false;
	for (int i = 1; i < gr.num_vertices() - 1; i++)
	{
		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;

		if (p2 - p1 >= min_exon_length) continue;
		if (gr.degree(i) <= 0) continue;

		for (pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			//if(gr.out_degree(s) <= 1) b = false;
			if (s != 0 && gr.get_vertex_info(s).rpos == p1) b = false;
			if (b == false) break;
		}
		for (pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			//if(gr.in_degree(t) <= 1) b = false;
			if (t != gr.num_vertices() - 1 && gr.get_vertex_info(t).lpos == p2) b = false;
			if (b == false) break;
		}

		if (b == false) continue;

		// only consider boundary small exons
		if (gr.edge(0, i).second == false && gr.edge(i, gr.num_vertices() - 1).second == false) continue;

		if (verbose >= 4 && fixed_bundle_idx == -1)
		{
			printf("remove_small_exons: ");
			write_clear_vertex(cout, i);
		}
		gr.clear_vertex(i);
		flag = true;
	}
	return flag;
}

bool bundle::remove_small_junctions()
{
	SE se;
	for (int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if (gr.degree(i) <= 0) continue;

		bool b = true;
		edge_iterator it1, it2;
		PEEI pei;
		int32_t p1 = gr.get_vertex_info(i).lpos;
		int32_t p2 = gr.get_vertex_info(i).rpos;
		double wi = gr.get_vertex_weight(i);

		// compute max in-adjacent edge
		double ws = 0;
		for (pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_vertex_weight(s);
			if (s == 0) continue;
			if (gr.get_vertex_info(s).rpos != p1) continue;
			if (w < ws) continue;
			ws = w;
		}

		// remove small in-junction
		for (pei = gr.in_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int s = e->source();
			double w = gr.get_edge_weight(e);
			if (s == 0) continue;
			if (gr.get_vertex_info(s).rpos == p1) continue;
			if (ws < 2.0 * w * w + 18.0) continue;
			if (wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

		// compute max out-adjacent edge
		double wt = 0;
		for (pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			int t = e->target();
			double w = gr.get_vertex_weight(t);
			if (t == gr.num_vertices() - 1) continue;
			if (gr.get_vertex_info(t).lpos != p2) continue;
			if (w < wt) continue;
			wt = w;
		}

		// remove small in-junction
		for (pei = gr.out_edges(i), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
		{
			edge_descriptor e = (*it1);
			double w = gr.get_edge_weight(e);
			int t = e->target();
			if (t == gr.num_vertices() - 1) continue;
			if (gr.get_vertex_info(t).lpos == p2) continue;
			if (ws < 2.0 * w * w + 18.0) continue;
			if (wi < 2.0 * w * w + 18.0) continue;

			se.insert(e);
		}

	}

	if (se.size() <= 0) return false;

	for (SE::iterator it = se.begin(); it != se.end(); it++)
	{
		edge_descriptor e = (*it);
		if (verbose >= 4 && fixed_bundle_idx == -1)
			printf("remove_samll_junctions: remove edge (%d, %d), weight = %.2lf\n", e->source(), e->target(), gr.get_edge_weight(e));
		gr.remove_edge(e);
	}

	return true;
}

VE bundle::compute_maximal_edges()
{
	typedef pair<double, edge_descriptor> PDE;
	vector<PDE> ve;

	undirected_graph ug;
	edge_iterator it1, it2;
	PEEI pei;
	for (int i = 0; i < gr.num_vertices(); i++) ug.add_vertex();
	for (pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		edge_descriptor e = (*it1);
		double w = gr.get_edge_weight(e);
		int s = e->source();
		int t = e->target();
		if (s == 0) continue;
		if (t == gr.num_vertices() - 1) continue;
		ug.add_edge(s, t);
		ve.push_back(PDE(w, e));
	}

	vector<int> vv = ug.assign_connected_components();

	sort(ve.begin(), ve.end());

	for (int i = 1; i < ve.size(); i++) assert(ve[i - 1].first <= ve[i].first);

	VE x;
	set<int> sc;
	for (int i = ve.size() - 1; i >= 0; i--)
	{
		edge_descriptor e = ve[i].second;
		double w = gr.get_edge_weight(e);
		if (w < 1.5) break;
		int s = e->source();
		int t = e->target();
		if (s == 0) continue;
		if (t == gr.num_vertices()) continue;
		int c1 = vv[s];
		int c2 = vv[t];
		assert(c1 == c2);
		if (sc.find(c1) != sc.end()) continue;
		x.push_back(e);
		sc.insert(c1);
	}
	return x;
}

bool bundle::keep_surviving_edges()
{
	set<int> sv1;
	set<int> sv2;
	SE se;
	edge_iterator it1, it2;
	PEEI pei;
	for (pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		int s = (*it1)->source();
		int t = (*it1)->target();
		double w = gr.get_edge_weight(*it1);
		int32_t p1 = gr.get_vertex_info(s).rpos;
		int32_t p2 = gr.get_vertex_info(t).lpos;
		if (w < min_surviving_edge_weight) continue;
		se.insert(*it1);
		sv1.insert(t);
		sv2.insert(s);
	}

	VE me = compute_maximal_edges();
	for (int i = 0; i < me.size(); i++)
	{
		edge_descriptor ee = me[i];
		se.insert(ee);
		sv1.insert(ee->target());
		sv2.insert(ee->source());
	}

	while (true)
	{
		bool b = false;
		for (SE::iterator it = se.begin(); it != se.end(); it++)
		{
			edge_descriptor e = (*it);
			int s = e->source();
			int t = e->target();
			if (sv1.find(s) == sv1.end() && s != 0)
			{
				edge_descriptor ee = gr.max_in_edge(s);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(s);
				sv2.insert(ee->source());
				b = true;
			}
			if (sv2.find(t) == sv2.end() && t != gr.num_vertices() - 1)
			{
				edge_descriptor ee = gr.max_out_edge(t);
				assert(ee != null_edge);
				assert(se.find(ee) == se.end());
				se.insert(ee);
				sv1.insert(ee->target());
				sv2.insert(t);
				b = true;
			}
			if (b == true) break;
		}
		if (b == false) break;
	}

	VE ve;
	for (pei = gr.edges(), it1 = pei.first, it2 = pei.second; it1 != it2; it1++)
	{
		if (se.find(*it1) != se.end()) continue;
		ve.push_back(*it1);
	}

	for (int i = 0; i < ve.size(); i++)
	{
		if (verbose >= 4 && fixed_bundle_idx == -1)
			printf("keep_surviving_edges: remove edge (%d, %d), weight = %.2lf\n", ve[i]->source(), ve[i]->target(), gr.get_edge_weight(ve[i]));
		gr.remove_edge(ve[i]);
	}

	if (ve.size() >= 1) return true;
	else return false;
}

bool bundle::remove_intron_contamination()
{
	bool flag = false;
	for (int i = 1; i < gr.num_vertices(); i++)
	{
		if (gr.in_degree(i) != 1) continue;
		if (gr.out_degree(i) != 1) continue;

		edge_iterator it1, it2;
		PEEI pei = gr.in_edges(i);
		it1 = pei.first;
		edge_descriptor e1 = (*it1);
		pei = gr.out_edges(i);
		it1 = pei.first;
		edge_descriptor e2 = (*it1);
		int s = e1->source();
		int t = e2->target();
		double wv = gr.get_vertex_weight(i);
		vertex_info vi = gr.get_vertex_info(i);

		if (s == 0) continue;
		if (t == gr.num_vertices() - 1) continue;
		if (gr.get_vertex_info(s).rpos != vi.lpos) continue;
		if (gr.get_vertex_info(t).lpos != vi.rpos) continue;

		PEB p = gr.edge(s, t);
		if (p.second == false) continue;

		edge_descriptor ee = p.first;
		double we = gr.get_edge_weight(ee);

		if (wv > we) continue;
		if (wv > max_intron_contamination_coverage) continue;

		if (verbose >= 4 && fixed_bundle_idx == -1)
		{
			printf("remove_intron_contamination: ");
			write_clear_vertex(cout, i);
		}
		gr.clear_vertex(i);
		flag = true;
	}
	return flag;
}

int bundle::refine_splice_graph()
{
	while (true)
	{
		bool b = false;
		for (int i = 1; i < gr.num_vertices() - 1; i++)
		{
			if (gr.degree(i) == 0) continue;
			if (gr.in_degree(i) >= 1 && gr.out_degree(i) >= 1) continue;
			if (verbose >= 4 && fixed_bundle_idx == -1)
			{
				printf("refine_splice_graph: ");
				write_clear_vertex(cout, i);
			}
			gr.clear_vertex(i);
			b = true;
		}
		if (b == false) break;
	}
	return 0;
}


int bundle::count_junctions() const
{
	int x = 0;
	for (int i = 0; i < junctions.size(); i++)
	{
		x += junctions[i].count;
	}
	return x;
}

int bundle::write(ostream &out, int index, int mode)
{
	out << "Bundle " << index << ": ";
	// statistic xs
	int n0 = 0, np = 0, nq = 0;
	for (int i = 0; i < hits.size(); i++)
	{
		if (hits[i].xs == '.') n0++;
		if (hits[i].xs == '+') np++;
		if (hits[i].xs == '-') nq++;
	}

	out << "tid = " << tid
		<< ", #hits = " << hits.size()
		<< ", #partial-exons = " << pexons.size()
		<< ", range = " << chrm << ":" << lpos << "-" << rpos
		<< ", orient = " << strand << " (" << n0 << ", " << np << ", " << nq << ")\n";

	if (mode == 0 && verbose <= 3) return 0;
	if (fixed_bundle_idx != -1) return 0;
	// verbose>=4 || verbose==-1;
	// print hits
	//for (int i = 0; i < hits.size(); i++) hits[i].write(out);

	// print regions
	for (int i = 0; i < regions.size(); i++) regions[i].write(out, i);

	// print junctions 
	for (int i = 0; i < junctions.size(); i++) junctions[i].write(out, chrm, i);

	// print partial exons
	for (int i = 0; i < pexons.size(); i++) pexons[i].write(out, i);

	// print hyper-edges
	hs.write(out);

	out << "\n";
	return 0;
}

int bundle::write_clear_vertex(ostream &out, int i)
{
	vertex_info vi = gr.get_vertex_info(i);
	out << "clear vertex " << i << ", weight=" << gr.get_vertex_weight(i) << ", length=" << vi.length << ", pos=[" << vi.lpos << ", " << vi.rpos << "), stddev=" << vi.stddev << "\n";
	out << "remove in-edges: ";
	PEEI pei = gr.in_edges(i);
	edge_iterator it1, it2;
	edge_descriptor e;
	for (it1 = pei.first, it2 = pei.second; it1 != it2; ++it1)
	{
		e = (*it1);
		out << "(" << e->source() << ", " << e->target() << ")=" << gr.get_edge_weight(e) << ", ";
	}
	out << "\nremove out-edges: ";
	pei = gr.out_edges(i);
	for (it1 = pei.first, it2 = pei.second; it1 != it2; ++it1)
	{
		e = (*it1);
		out << "(" << e->source() << ", " << e->target() << ")=" << gr.get_edge_weight(e) << ", ";
	}
	out << "\n";
	return 0;
}
/*
int bundle::output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const
{
for(int i = 0; i < p.size(); i++)
{
string tid = gid + "." + tostring(i);
output_transcript(fout, p[i], gid, tid);
}
return 0;
}

int bundle::output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const
{
fout.precision(2);
fout<<fixed;

const vector<int> &v = p.v;
double coverage = p.abd;		// number of molecular

assert(v[0] == 0);
assert(v[v.size() - 1] == pexons.size() + 1);
if(v.size() < 2) return 0;

int ss = v[1];
int tt = v[v.size() - 2];
int32_t ll = pexons[ss - 1].lpos;
int32_t rr = pexons[tt - 1].rpos;

fout<<chrm.c_str()<<"\t";		// chromosome name
fout<<algo.c_str()<<"\t";		// source
fout<<"transcript\t";			// feature
fout<<ll + 1<<"\t";				// left position
fout<<rr<<"\t";					// right position
fout<<1000<<"\t";				// score, now as abundance
fout<<strand<<"\t";				// strand
fout<<".\t";					// frame
fout<<"gene_id \""<<gid.c_str()<<"\"; ";
fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
fout<<"coverage \""<<coverage<<"\";"<<endl;

join_interval_map jmap;
for(int k = 1; k < v.size() - 1; k++)
{
const partial_exon &r = pexons[v[k] - 1];
jmap += make_pair(ROI(r.lpos, r.rpos), 1);
}

int cnt = 0;
for(JIMI it = jmap.begin(); it != jmap.end(); it++)
{
fout<<chrm.c_str()<<"\t";			// chromosome name
fout<<algo.c_str()<<"\t";			// source
fout<<"exon\t";						// feature
fout<<lower(it->first) + 1<<"\t";	// left position
fout<<upper(it->first)<<"\t";		// right position
fout<<1000<<"\t";					// score
fout<<strand<<"\t";					// strand
fout<<".\t";						// frame
fout<<"gene_id \""<<gid.c_str()<<"\"; ";
fout<<"transcript_id \""<<tid.c_str()<<"\"; ";
fout<<"exon_number \""<<++cnt<<"\"; ";
fout<<"coverage \""<<coverage<<"\";"<<endl;
}
return 0;
}

int bundle::output_transcripts(vector<transcript> &trsts, const vector<path> &p, const string &gid) const
{
trsts.clear();
for(int i = 0; i < p.size(); i++)
{
string tid = gid + "." + tostring(i);
transcript trst;
output_transcript(trst, p[i], gid, tid);
trsts.push_back(trst);
}
return 0;
}

int bundle::output_transcripts(gene &gn, const vector<path> &p, const string &gid) const
{
for(int i = 0; i < p.size(); i++)
{
string tid = gid + "." + tostring(i);
transcript trst;
output_transcript(trst, p[i], gid, tid);
gn.add_transcript(trst);
}
return 0;
}

int bundle::output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const
{
trst.seqname = chrm;
trst.source = algo;
trst.gene_id = gid;
trst.transcript_id = tid;
trst.coverage = p.abd;
trst.strand = strand;

const vector<int> &v = p.v;
join_interval_map jmap;
for(int k = 1; k < v.size() - 1; k++)
{
const partial_exon &r = pexons[v[k] - 1];
jmap += make_pair(ROI(r.lpos, r.rpos), 1);
}
for(JIMI it = jmap.begin(); it != jmap.end(); it++)
{
trst.add_exon(lower(it->first), upper(it->first));
}
return 0;
}*/
