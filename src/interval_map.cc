/*
This file was modified from a file named interval_map.cc in Scallop
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

#include "interval_map.h"

int create_split(split_interval_map &imap, int32_t p)
{
	SIMI it = imap.find(p);
	if(it == imap.end()) return 0;
	int32_t l = lower(it->first);
	int32_t r = upper(it->first);
	int32_t w = it->second;
	assert(l <= p);
	assert(r >= p);
	if(l == p || r == p) return 0;
	imap -= make_pair(ROI(l, r), w);
	imap += make_pair(ROI(l, p), w);
	imap += make_pair(ROI(p, r), w);
	return 0;
}

int compute_overlap(const split_interval_map &imap, int32_t p)
{
	SIMI it = imap.find(p);
	if(it == imap.end()) return 0;
	return it->second;
}

SIMI locate_right_iterator(const split_interval_map &imap, int32_t x)
{
	return imap.upper_bound(ROI(x - 1, x));
}

SIMI locate_left_iterator(const split_interval_map &imap, int32_t x)
{
	SIMI it = imap.lower_bound(ROI(x - 1, x));
	if(it == imap.end() && it == imap.begin()) return it;
	if(it == imap.end()) it--;

	while(upper(it->first) > x)
	{
		if(it == imap.begin()) return imap.end();
		it--;
	}
	return it;
}

PSIMI locate_boundary_iterators(const split_interval_map &imap, int32_t x, int32_t y)
{
	SIMI lit, rit;//????iterator
	lit = locate_right_iterator(imap, x);
	if(lit == imap.end() || upper(lit->first) > y) lit = imap.end();

	rit = locate_left_iterator(imap, y);
	if(rit == imap.end() || lower(rit->first) < x) rit = imap.end();

	if(lit == imap.end()) assert(rit == imap.end());
	if(rit == imap.end() && lit != imap.end()) 
	{
		//printf("x = %d, y = %d, lit = [%d, %d)\n", x, y, lower(lit->first), upper(lit->first));
		assert(lit == imap.end());
	}

	return PSIMI(lit, rit); 
}

int32_t compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int32_t x = it->second;
		if(x > s) s = x;
	}

	if(q != imap.end())
	{
		int32_t x = q->second;
		if(x > s) s = x;
	}

	return s;
}

int32_t compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		int l = lower(it->first);
		int u = upper(it->first);
		assert(u > l);

		//printf(" AA add [%d, %d) : %d\n", lower(it->first), upper(it->first), it->second);

		s += (u - l) * it->second;
	}
	if(q != imap.end())
	{
		//printf(" BB add [%d, %d) : %d\n", lower(q->first), upper(q->first), q->second);
		s += (upper(q->first) - lower(q->first)) * q->second;
	}
	return s;
}

int32_t compute_coverage(const split_interval_map &imap, SIMI &p, SIMI &q)
{
	if(p == imap.end()) return 0;

	int32_t s = 0;
	for(SIMI it = p; it != q; it++)
	{
		s += upper(it->first) - lower(it->first);
	}

	if(q != imap.end()) s += upper(q->first) - lower(q->first);

	return s;
}

int evaluate_rectangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev)
{
	// split_interval_map(interval_map.h)
	//typedef icl::split_interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map;
	//typedef split_interval_map::const_iterator SIMI;
	//typedef pair<SIMI, SIMI> PSIMI;
	ave = 0;
	dev = 1.0;

	PSIMI pei = locate_boundary_iterators(imap, ll, rr);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == imap.end()) return 0;
	if(rit == imap.end()) return 0;

	ave = 1.0 * compute_sum_overlap(imap, lit, rit) / (rr - ll);
	//printf("compute average %d-%d = %.2lf\n", ll, rr, ave);

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		var += (it->second - ave) * (it->second - ave) * (upper(it->first) - lower(it->first));

		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	//if(dev < 1.0) dev = 1.0;

	return 0;
}

int evaluate_triangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev)
{
	ave = 0;
	dev = 1.0;

	PSIMI pei = locate_boundary_iterators(imap, ll, rr);
	SIMI lit = pei.first, rit = pei.second;

	if(lit == imap.end()) return 0;
	if(rit == imap.end()) return 0;

	vector<double> xv;
	vector<double> yv;
	double xm = 0;
	double ym = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (lower(it->first) + upper(it->first)) / 2.0;
		double yi = it->second;
		xv.push_back(xi);
		yv.push_back(yi);
		xm += xi;
		ym += yi;
		if(it == rit) break;
	}

	xm /= xv.size();
	ym /= yv.size();

	double f1 = 0;
	double f2 = 0;
	for(int i = 0; i < xv.size(); i++)
	{
		f1 += (xv[i] - xm) * (yv[i] - ym);
		f2 += (xv[i] - xm) * (xv[i] - xm);
	}

	double b1 = f1 / f2;
	double b0 = ym - b1 * xm;

	double a1 = b1 * rr + b0;
	double a0 = b1 * ll + b0;
	ave = (a1 > a0) ? a1 : a0;

	double var = 0;
	for(SIMI it = lit; ; it++)
	{
		assert(upper(it->first) > lower(it->first));
		double xi = (upper(it->first) + lower(it->first)) / 2.0;
		double yi = b1 * xi + b0;
		var += (it->second - yi) * (it->second - yi) * (upper(it->first) - lower(it->first));
		if(it == rit) break;
	}

	dev = sqrt(var / (rr - ll));
	if(dev < 1.0) dev = 1.0;

	return 0;
}

int test_split_interval_map()
{
	split_interval_map imap;

	imap += make_pair(ROI(6, 7), 3);
	imap += make_pair(ROI(1, 3), 3);
	imap += make_pair(ROI(1, 2), 1);
	imap += make_pair(ROI(2, 5), 2);

	create_split(imap, 4);

	SIMI it;
	
	for(it = imap.begin(); it != imap.end(); it++)
	{
		printf("interval: [%d,%d) -> %d\n", lower(it->first), upper(it->first), it->second);
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.find(i);
		if(it == imap.end())
		{
			printf("find %d: does not exist\n", i);
		}
		else
		{
			printf("find %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.lower_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("lower bound %d: does not exist\n", i);
		}
		else
		{
			printf("lower bound %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 1; i <= 7; i++)
	{
		it = imap.upper_bound(ROI(i, i + 1));

		if(it == imap.end())
		{
			printf("upper bound %d: does not exist\n", i);
		}
		else
		{
			printf("upper bound %d: [%d,%d) -> %d\n", i, lower(it->first), upper(it->first), it->second);
		}
	}

	for(int i = 0; i <= 8; i++)
	{
		for(int j = i; j <= 8; j++)
		{
			pair<SIMI, SIMI> p = locate_boundary_iterators(imap, i, j);
			int s = compute_coverage(imap, p.first, p.second);
			printf("coverage [%d,%d) = %d\n", i, j, s);
		}
	}

	return 0;
}


