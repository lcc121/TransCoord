/*
This file was modified from a file named previewer.cc in Scallop
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
#include <cassert>
#include <sstream>

#include "previewer.h"
#include "config.h"

previewer::previewer()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
}

previewer::~previewer()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int previewer::preview()
{
	int total = 0;
	int single = 0;
	int paired = 0;

	int first = 0;
	int second = 0;
	vector<int> sp1;
	vector<int> sp2;

    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(total >= max_preview_reads) break;
		if(sp1.size() >= max_preview_spliced_reads && sp2.size() >= max_preview_spliced_reads) break;

		bam1_core_t &p = b1t->core;

		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// qstrandary alignment
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality
		if(p.n_cigar < 1) continue;												// should never happen

		total++;

		hit ht(b1t);
		ht.set_tags(b1t);

		if((ht.flag & 0x1) >= 1) paired ++;
		if((ht.flag & 0x1) <= 0) single ++;

		if(ht.xs == '.') continue;
		if(ht.xs == '+' && sp1.size() >= max_preview_spliced_reads) continue;
		if(ht.xs == '-' && sp2.size() >= max_preview_spliced_reads) continue;

		// predicted strand
		char xs = '.';

		// for paired read
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '-';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) >= 1 && (ht.flag & 0x20) <= 0 && (ht.flag & 0x40) >= 1 && (ht.flag & 0x80) <= 0) xs = '+';
		if((ht.flag & 0x1) >= 1 && (ht.flag & 0x10) <= 0 && (ht.flag & 0x20) >= 1 && (ht.flag & 0x40) <= 0 && (ht.flag & 0x80) >= 1) xs = '+';

		// for single read
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) <= 0) xs = '-';
		if((ht.flag & 0x1) <= 0 && (ht.flag & 0x10) >= 1) xs = '+';

		if(xs == '+' && xs == ht.xs) sp1.push_back(1);
		if(xs == '-' && xs == ht.xs) sp2.push_back(1);
		if(xs == '+' && xs != ht.xs) sp1.push_back(2);
		if(xs == '-' && xs != ht.xs) sp2.push_back(2);
	}

	int sp = sp1.size() < sp2.size() ? sp1.size() : sp2.size();

	for(int k = 0; k < sp; k++)
	{
		if(sp1[k] == 1) first++;
		if(sp2[k] == 1) first++;
		if(sp1[k] == 2) second++;
		if(sp2[k] == 2) second++;
	}

	vector<string> vv;
	vv.push_back("empty");
	vv.push_back("unstranded");
	vv.push_back("first");
	vv.push_back("second");

	int s1 = UNSTRANDED;
	if(sp >= min_preview_spliced_reads && first > preview_infer_ratio * 2.0 * sp) s1 = FR_FIRST;
	if(sp >= min_preview_spliced_reads && second > preview_infer_ratio * 2.0 * sp) s1 = FR_SECOND;

	if(verbose >= 1)
	{
		printf("preview: reads = %d, single = %d, paired = %d, spliced reads = %d, first = %d, second = %d, inferred library_type = %s, given library_type = %s\n",
			total, single, paired, sp, first, second, vv[s1 + 1].c_str(), vv[library_type + 1].c_str());
	}

	if(library_type == EMPTY) library_type = s1;

	return 0;
}
