/*
This file was modified from a file named subsetsum.h in Scallop
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

#ifndef __SUBSETSUM4_H__
#define __SUBSETSUM4_H__

#include <vector>
#include "equation.h"

using namespace std;

typedef pair<int, int> PI;

// partition s and t into s1/s2 and t1/t2
// such that sum(s1) is close to sum(t1)
// AND sum(s2) is close to sum(t2)
class subsetsum
{
public:
	subsetsum(const vector<PI> &s, const vector<PI> &t);

private:
	vector<PI> source;					// given input
	vector<PI> target;					// given target numbers
	int ubound1;						// ubound for source
	int ubound2;						// ubound for target
	vector< vector<int> > table1;		// dp table2
	vector< vector<int> > table2;		// dp table2

public:
	equation eqn;

public:
	int solve();
	int print();
	static int test();

private:
	int rescale();
	int init(const vector<PI> &vv, vector< vector<int> > &table, int ubound);
	int fill(const vector<PI> &vv, vector< vector<int> > &table, int ubound);
	int backtrace(int vi, const vector<PI> &vv, const vector< vector<int> > &table, vector<int> &ss);
	int optimize();
};

#endif
