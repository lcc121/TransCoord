/*
This file was modified from a file named equation.cc in Scallop
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

#include "equation.h"
#include "util.h"
#include <cstdio>

equation::equation()
{
	e = 0;
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(double _e)
	:e(_e)
{
	f = 0;
	a = 0;
	d = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t)
	: s(_s), t(_t)
{
	e = 0;
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

equation::equation(const vector<int> &_s, const vector<int> &_t, double _e)
	: s(_s), t(_t), e(_e)
{
	f = 0;
	d = 0;
	a = 0;
	w = 0;
}

int equation::clear()
{
	s.clear();
	t.clear();
	e = 0;
	f = 0;
	a = 0;
	d = 0;
	w = 0;
	return 0;
}

int equation::print(int index) const
{
	printf("equation %3d: (%2lu, %2lu) edges, error = %.2lf, f = %d, w = %d, adjacent = %2d, distant = %2d. ", 
			index, s.size(), t.size(), e, f, w, a, d);

	printf("S = ( ");
	//printv(s);
	writev(cout, s);
	printf("), T = ( ");
	//printv(t);
	writev(cout, s);
	printf(")\n");
	
	return 0;
}

bool equation_cmp1(const equation &x, const equation &y) 
{
	if(x.e < y.e - 0.00001) return true;
	else if(x.e > y.e + 0.00001) return false;
	else if(x.s.size() + x.t.size() < y.s.size() + y.t.size()) return true;
	else return false;
}

bool equation_cmp2(const equation &x, const equation &y) 
{
	if(x.f > y.f) return true;
	if(x.f < y.f) return false;
	if(x.d < y.d) return true;
	if(x.d > y.d) return false;
	if(x.s.size() + x.t.size() < y.s.size() + y.t.size()) return true;
	if(x.s.size() + x.t.size() > y.s.size() + y.t.size()) return false;
	if(x.w < y.w) return true;
	return false;
}
