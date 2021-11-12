/*
This file was modified from a file named gene.cc in Scallop
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

#include <map>
#include <algorithm>
#include <cassert>
#include "gene.h"

int gene::clear()
{
	transcripts.clear();
	t2i.clear();
	return 0;
}

int gene::assign(const vector<transcript> &v)
{
	clear();
	for(int i = 0; i < v.size(); i++)
	{
		add_transcript(v[i]);
	}
	return 0;
}

int gene::add_transcript(const item&e)
{
	assert(e.feature == "transcript");
	if(t2i.find(e.transcript_id) == t2i.end())
	{
		t2i.insert(pair<string, int>(e.transcript_id, t2i.size()));
		transcripts.push_back(transcript(e));
	}
	else
	{
		int k = t2i[e.transcript_id];
		transcripts[k].assign(e);
	}
	return 0;
}

int gene::add_transcript(const transcript &t)
{
	assert(t2i.find(t.transcript_id) == t2i.end());
	t2i.insert(pair<string, int>(t.transcript_id, t2i.size()));
	transcripts.push_back(t);
	return 0;
}

int gene::add_exon(const item&e)
{
	assert(e.feature == "exon");
	if(t2i.find(e.transcript_id) != t2i.end())
	{
		int k = t2i[e.transcript_id];
		transcripts[k].add_exon(e);
	}
	else
	{
		transcript t;
		t.assign(e);
		t.add_exon(e);
		add_transcript(t);
	}
	return 0;
}

int gene::set_gene_id(const string &id) 
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].gene_id = id;
	}
	return 0;
}

int gene::filter_single_exon_transcripts()
{
	vector<transcript> vv;
	t2i.clear();
	for(int i = 0; i < transcripts.size(); i++)
	{
		if(transcripts[i].exons.size() <= 1) continue;
		t2i.insert(pair<string, int>(transcripts[i].transcript_id, vv.size()));
		vv.push_back(transcripts[i]);
	}
	transcripts = vv;
	return 0;
}

int gene::filter_low_coverage_transcripts(double min_coverage)
{
	vector<transcript> vv;
	t2i.clear();
	for(int i = 0; i < transcripts.size(); i++)
	{
		if(transcripts[i].coverage < min_coverage) continue;
		t2i.insert(pair<string, int>(transcripts[i].transcript_id, vv.size()));
		vv.push_back(transcripts[i]);
	}
	transcripts = vv;
	return 0;
}

int gene::sort()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].sort();
	}
	return 0;
}

int gene::shrink()
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].shrink();
	}
	return 0;
}

int gene::assign_RPKM(double factor)
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].assign_RPKM(factor);
	}
	return 0;
}

set<int32_t> gene::get_exon_boundaries() const
{
	set<int32_t> s;
	for(int k = 0; k < transcripts.size(); k++)
	{
		const vector<PI32> &v = transcripts[k].exons;
		for(int i = 0; i < v.size(); i++)
		{
			int32_t l = v[i].first;
			int32_t r = v[i].second;
			s.insert(l);
			s.insert(r);
		}
	}
	return s;
}

PI32 gene::get_bounds() const
{
	PI32 pp(-1, -1);
	for(int k = 0; k < transcripts.size(); k++)
	{
		PI32 p = transcripts[k].get_bounds();
		if(pp.first == -1 || pp.first > p.first) pp.first = p.first;
		if(pp.second == -1 || pp.second < p.second) pp.second = p.second;
	}
	return pp;
}

string gene::get_seqname() const
{
	if(transcripts.size() == 0) return "";
	else return transcripts[0].seqname;
}

string gene::get_gene_id() const
{
	if(transcripts.size() == 0) return "";
	else return transcripts[0].gene_id;
}

int gene::write(ofstream &fout) const
{
	for(int i = 0; i < transcripts.size(); i++)
	{
		transcripts[i].write(fout);
	}
	return 0;
}
