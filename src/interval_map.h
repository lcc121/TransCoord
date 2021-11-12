/*
This file was modified from a file named interval_map.h in Scallop
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

#ifndef __INTERVAL_MAP_H__
#define __INTERVAL_MAP_H__

// boost::interval map
#include "boost/icl/interval_map.hpp"
#include "boost/icl/split_interval_map.hpp"

#include <vector>

using namespace boost;
using namespace std;

typedef icl::right_open_interval<int32_t> ROI;

// join interval map
typedef icl::interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> join_interval_map;
typedef join_interval_map::const_iterator JIMI;
typedef pair<JIMI, JIMI> PJIMI;

// split interval map
typedef icl::split_interval_map<int32_t, int32_t, icl::partial_absorber, less, icl::inplace_plus, icl::inter_section, ROI> split_interval_map;
typedef split_interval_map::const_iterator SIMI;
typedef pair<SIMI, SIMI> PSIMI;

// if p is inside an interval, split this interval into 2
int create_split(split_interval_map &imap, int32_t p);

// return the overlap at position p
int compute_overlap(const split_interval_map &imap, int32_t p);

// find the leftmost iterator whose upper posistion <= x
SIMI locate_right_iterator(const split_interval_map &imap, int32_t x);

// find the rightmost interval whose lower position >= x
SIMI locate_left_iterator(const split_interval_map &imap, int32_t x);

// locate boundary iterators
PSIMI locate_boundary_iterators(const split_interval_map &imap, int32_t x, int32_t y);

// return the sum of the lengths of intervals from p to q (include q)
int compute_coverage(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the maximum overlap of the intervals from p to q (include q)
int compute_max_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

// return the sum of the overlap of the intervals from p to q (include q)
int compute_sum_overlap(const split_interval_map &imap, SIMI &p, SIMI &q);

// evaluate a region
int evaluate_rectangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev);
int evaluate_triangle(const split_interval_map &imap, int ll, int rr, double &ave, double &dev);

// testing
int test_split_interval_map();

#endif
