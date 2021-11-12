/*
This file was modified from a file named path.cc in Scallop
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

#include "path.h"

#include <cassert>
#include <cstdio>

path::path()
{
	v.clear();
	abd = 0;
	reads = 0;
	//length = 0;
}

path::~path()
{}

int path::clear()
{
	v.clear();
	abd = 0;
	reads = 0;
	//length = 0;
	return 0;
}

int path::print(int index) const
{
	if(v.size() == 0) return 0;
	//printf("path %d: abundance = %.2lf, length = %d, reads = %.2lf, vertices = ", index, abd, length, reads);
	printf("path %d: abundance = %f, length = %d, reads = %f, vertices = ", index, abd, length, reads);
	for(int i = 0; i < v.size() - 1; i++)
	{
		printf("%d, ", v[i]);
	}
	printf("%d\n", v[v.size() - 1]);
	return 0;
}

vector<int> path::index(int n) const
{
	vector<int> vv;
	vv.resize(n, -1);
	for(int i = 1; i < v.size(); i++)
	{
		int s = v[i - 1];
		int t = v[i];
		assert(s >= 0 && s < n);
		assert(t >= 0 && t < n);
		vv[s] = t;
	}
	return vv;
}
