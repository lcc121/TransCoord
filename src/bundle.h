/*
This file was modified from a file named bundle.h in Scallop
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
#ifndef __BUNDLE_H__
#define __BUNDLE_H__

#include "hit.h"
#include "interval_map.h"
#include "junction.h"
#include "region.h"
#include "partial_exon.h"
#include "splice_graph.h"
#include "hyper_set.h"
#include "transcript.h"

using namespace std;

class bundle
{
public:
	bundle();
	~bundle();

public:
	int32_t tid;					// chromosome ID
	string chrm;					// chromosome name
	int32_t lpos;					// the leftmost boundary on reference
	int32_t rpos;					// the rightmost boundary on reference
	char strand;					// strandness
	vector<hit> hits;				// hits
	split_interval_map mmap;		// matched interval map
	split_interval_map imap;		// indel interval map

	vector<junction> junctions;		// splice junctions
	vector<region> regions;			// regions
	vector<partial_exon> pexons;	// partial exons
	vector<bool> regional;			// if a pe is regional
	split_interval_map pmap;		// partial exon map                                                                                                                                                                                                             
	splice_graph gr;				// splice graph
	hyper_set hs;					// hyper edges

public:
	int add_hit(const hit & ht);
	int build();
	int clear();
	int clear_memory();
	//int output_transcripts(ofstream &fout, const vector<path> &p, const string &gid) const;	
	//int output_transcripts(gene &gn, const vector<path> &p, const string &gid) const;	
	//int output_transcripts(vector<transcript> &trsts, const vector<path> &p, const string &gid) const;	
	//int output_transcript(ofstream &fout, const path &p, const string &gid, const string &tid) const;	
	//int output_transcript(transcript &trst, const path &p, const string &gid, const string &tid) const;	
	int count_junctions() const;
	int write(ostream & out, int index, int mode);
	int write_clear_vertex(ostream & out, int i);

private:
	// check and init
	int check_left_ascending();
	int check_right_ascending();
	int compute_strand();

	// splice graph
	int build_junctions();
	int correct_junctions();
	int build_regions();
	int build_partial_exons();
	int link_partial_exons();
	int build_splice_graph();
	int build_partial_exon_map();
	int locate_left_partial_exon(int32_t x);
	int locate_right_partial_exon(int32_t x);

	// revise splice graph
	VE compute_maximal_edges();
	int revise_splice_graph();
	int refine_splice_graph();
	bool keep_surviving_edges();
	bool extend_boundaries();
	bool remove_small_junctions();
	bool remove_small_exons();
	bool remove_inner_boundaries();
	bool remove_intron_contamination();

	// super edges
	int build_hyper_edges2();			// paired end
	bool bridge_read(int x, int y, vector<int> &s);

};

#endif
