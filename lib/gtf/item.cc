/*
This file was modified from a file named item.cc in Scallop
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

#include "item.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <cassert>
#include <cstdio>
#include <cmath>

item::item(const string &s)
{
	parse(s);
}

int item::parse(const string &s)
{
	char buf[10240];
	stringstream sstr(s);
	sstr>>buf;
	seqname.assign(buf);
	sstr>>buf;
	source.assign(buf);
	sstr>>buf;
	feature.assign(buf);
	sstr>>start>>end;
	start--;			// TODO gtf: (from 1, both inclusive)
	sstr>>buf;
	if(buf[0] == '.') score = -1;
	else score = atof(buf);
	sstr>>buf;
	strand = buf[0];
	sstr>>buf;
	frame = buf[0];

	char buf2[10240];
	coverage = 0;
	RPKM = 0;
	TPM = 0;
	while(sstr.eof() == false)//end of file
	{
		sstr>>buf;
		/*
		sstr.getline(buf2, 10240, '"');
		sstr.getline(buf2, 10240, '"');
		*/
		sstr.getline(buf2, 10240, ';');
		string v(buf2);
		int k1 = v.find('"');
		int k2 = v.rfind('"');
		if(k1 >= 0 && k1 < k2 && k2 < v.size()) v = v.substr(k1 + 1, k2 - k1 - 1);

		if(string(buf) == "" || v == "") break;

		//printf("|%s|%s|\n", buf, v.c_str());

		if(string(buf) == "transcript_id") transcript_id = v;
		else if(string(buf) == "transcript_type") transcript_type = v;
		else if(string(buf) == "gene_type") gene_type = v;
		else if(string(buf) == "gene_id") gene_id = v;
		else if(string(buf) == "cov") coverage = atof(v.c_str());
		else if(string(buf) == "coverage") coverage = atof(v.c_str());
		else if(string(buf) == "expression") coverage = atof(v.c_str());
		else if(string(buf) == "expr") coverage = atof(v.c_str());
		else if(string(buf) == "TPM") TPM = atof(v.c_str());
		else if(string(buf) == "RPKM") RPKM = atof(v.c_str());
		else if(string(buf) == "FPKM") FPKM = atof(v.c_str());

		//sstr.getline(buf2, 10240, ';');
	}

	return 0;
}

int item::print() const
{
	printf("%s\t%s\t%s\t%d\t%d\t%.1lf\t%c\t%c\ttranscript_id \"%s\"; gene_id \"%s\"; coverage \"%.2f\"; RPKM \"%.2lf\"\n",
			seqname.c_str(), source.c_str(), feature.c_str(), start, end, score, strand, frame,
			transcript_id.c_str(), gene_id.c_str(), coverage, RPKM);
	return 0;
}

bool item::operator<(const item &ge) const
{
	if(start < ge.start) return true;
	else return false;
}

int item::length() const
{
	assert(end > start);
	return end - start;
}
