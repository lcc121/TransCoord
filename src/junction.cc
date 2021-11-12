/*
This file was modified from a file named junction.cc in Scallop
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

#include <cstdio>
#include "junction.h"
#include "config.h"

junction::junction()
{}

junction::junction(int64_t _p)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = 0;
	strand = '.';
	lexon = -1;
	rexon = -1;
	nm = 0;
}

junction::junction(int64_t _p, int _c)
{
	lpos = high32(_p);
	rpos = low32(_p);
	count = _c;
	strand = '.';
	lexon = -1;
	rexon = -1;
	nm = 0;
}

junction::junction(const junction &sp)
{
	lpos = sp.lpos;
	rpos = sp.rpos;
	count = sp.count;
	lexon = sp.lexon;
	rexon = sp.rexon;
	strand = sp.strand;
	nm = sp.nm;
}

bool junction::operator<(const junction &x) const
{
	if(lpos <= x.lpos) return true;
	else return false;
}

int junction::print(const string &chrm, int index) const
{
	printf("junction %d: region = %s:%d-%d, %d -> %d, length = %d, count = %d, strand = %c, nm = %d\n", 
			index, chrm.c_str(), lpos, rpos, lexon, rexon, rpos - lpos, count, strand, nm);
	return 0;
}

int junction::write(ostream &out, const string &chrm, int index) const
{
	out << "junction " << index << ": region=" << chrm << ":" << lpos << "-" << rpos << ", " << lexon << "->" << rexon << ", length=" << rpos - lpos << ", count=" << count << ", strand=" << strand << ", nm=" << nm << "\n";
	return 0;
}

bool junction_cmp_length(const junction &x, const junction &y)
{
	int32_t p1 = x.rpos - x.lpos;
	int32_t p2 = y.rpos - y.lpos;
	if(p1 < p2) return true;
	else return false;
}
