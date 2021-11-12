/*
This file was modified from a file named partial_exon.cc in Scallop
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

#include "partial_exon.h"
#include "util.h"
#include <cstdio>

partial_exon::partial_exon(int32_t _lpos, int32_t _rpos, int _ltype, int _rtype)
	: lpos(_lpos), rpos(_rpos), ltype(_ltype), rtype(_rtype)
{
}

string partial_exon::label() const
{
	string l = tostring((lpos + 1) % 100000);
	string r = tostring(rpos % 100000);
	return (l + "-" + r);
}

int partial_exon::print(int index) const
{
	printf("partial_exon %d: [%d-%d), type = (%d, %d), length = %d, ave-abd = %.1lf, std-abd = %.1lf\n",
			index, lpos, rpos, ltype, rtype, rpos - lpos, ave, dev);
	return 0;
}

int partial_exon::write(ostream &out, int index)
{
	out << "partial_exon " << index << ": [" << lpos << ", " << rpos << "), type=(" << ltype << ", " << rtype << "), length=" << rpos - lpos << ", ave-abd=" << ave << ", std-abd=" << dev << "\n";
	return 0;
}