/*
This file was modified from a file named assembler.cc in Scallop
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
#include <cstring> 

#include "config.h"
#include "assembler.h"
#include "super_graph.h"
#include "filter.h"

assembler::assembler()
{
    sfn = sam_open(input_file.c_str(), "r");
    hdr = sam_hdr_read(sfn);
    b1t = bam_init1();
	index = 0;
	terminate = false;
	qlen = 0;
	qcnt = 0;
}

assembler::~assembler()
{
    bam_destroy1(b1t);
    bam_hdr_destroy(hdr);
    sam_close(sfn);
}

int assembler::assemble()
{
	assert(output_file != "");
    while(sam_read1(sfn, hdr, b1t) >= 0)
	{
		if(terminate == true) return 0;

		bam1_core_t &p = b1t->core;

		if(p.tid < 0) continue;
		if((p.flag & 0x4) >= 1) continue;										// read is not mapped
		if((p.flag & 0x100) >= 1 && use_second_alignment == false) continue;	// secondary alignment(default=false)
		if(p.n_cigar > max_num_cigar) continue;									// ignore hits with more than max-num-cigar types
		if(p.qual < min_mapping_quality) continue;								// ignore hits with small quality(default=1)
		if(p.n_cigar < 1) continue;												// should never happen

		hit ht(b1t);
		ht.set_tags(b1t);
		ht.set_strand();

		//if(ht.nh >= 2 && p.qual < min_mapping_quality) continue;
		//if(ht.nm > max_edit_distance) continue;

		qlen += ht.qlen;
		qcnt += 1;

		// truncate
		if(ht.tid != b1.tid || ht.pos > b1.rpos + min_bundle_gap)
		{
			process(b1);
			b1.clear();
		}
		if(ht.tid != b2.tid || ht.pos > b2.rpos + min_bundle_gap)
		{
			process(b2);
			b2.clear();
		}

		//printf("read strand = %c, xs = %c, ts = %c\n", ht.strand, ht.xs, ht.ts);

		// add hit
		if(uniquely_mapped_only == true && ht.nh != 1) continue;
		if(library_type != UNSTRANDED && ht.strand == '+' && ht.xs == '-') continue;
		if(library_type != UNSTRANDED && ht.strand == '-' && ht.xs == '+') continue;
		if(library_type != UNSTRANDED && ht.strand == '.' && ht.xs != '.') ht.strand = ht.xs;
		if(library_type != UNSTRANDED && ht.strand == '+') b1.add_hit(ht);
		if(library_type != UNSTRANDED && ht.strand == '-') b2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') b1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '.') b2.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '+') b1.add_hit(ht);
		if(library_type == UNSTRANDED && ht.xs == '-') b2.add_hit(ht);
	}

	process(b1);
	process(b2);
	b1.clear();
	b2.clear();
	
	assign_RPKM();

	filter ft(trsts);
	if (filter_tc.find('c') != string::npos) ft.merge_single_exon_transcripts();
	trsts = ft.trs;
	if (!output_all_paths) write();

	/*write_lp_stats(cout);
	ofstream fout(stats_file.c_str(), ios_base::app);
	if (fout.fail()) exit(3);
	write_lp_stats(fout);
	fout.close();*/
	
	return 0;
}

int assembler::process(bundle &b)
{
	if (b.hits.size() < min_num_hits_in_bundle) return 0;
	if (b.tid < 0) return 0;

	// build
	char buf[1024];
	strcpy(buf, hdr->target_name[b.tid]);
	b.chrm = string(buf);
	b.build();
	b.write(cout, index, 0);
	
	// assemble
	if (fixed_bundle_idx != -1 && index > fixed_bundle_idx) exit(0);
	if (fixed_bundle_idx != -1 && fixed_bundle_idx != index)
	{
		index++;
		return 0;
	}
	if (verbose == -1)
	{
		ofstream fout(log_file.c_str(), ios_base::trunc);
		if (fout.fail()) exit(2);
		b.write(fout, index, 1);
		fout.close();
	}
	b.clear_memory();
	assemble(b.gr, b.hs);
	index++;

	return 0;
}

int assembler::assemble(const splice_graph &gr0, const hyper_set &hs0)
{
	super_graph sg(gr0, hs0);
	sg.build();
	
	ofstream fout;
	if (verbose == -1)
	{
		fout.open(log_file.c_str(), ios_base::app);
		if (fout.fail()) exit(2);
	}

	vector<transcript> gv;
	for(int k = 0; k < sg.subs.size(); k++)
	{
		string gid = "gene." + tostring(index) + "." + tostring(k);
		if(fixed_gene_name != "" && gid != fixed_gene_name) continue;
		if(verbose >= 2 && (k == 0 || fixed_gene_name != "")) sg.print();
		if (verbose == -1 && (k == 0 || fixed_gene_name != "")) sg.write(fout);

		splice_graph &gr = sg.subs[k];
		hyper_set &hs = sg.hss[k];
		if (verbose == -1)
		{
			gr.write(fout, k);
			hs.write(fout);
		}
		if (verbose >= 2)
		{
			gr.write(cout, k);
			hs.write(cout);
		}

		if(determine_regional_graph(gr) == true) continue;
		if(gr.num_edges() <= 0) continue;

		gr.gid = gid;
		if (verbose == -1 || verbose >= 1) printf("\nprocess super_graph.subs[%d]: %s\n", k, gid.c_str());
		if (true)
		{
			subgraph subg(gr, hs);
			if (verbose >= 1) subg.write(cout, k);
			if (verbose == -1) subg.write(fout, k);
			subg.assemble();
			if (verbose >= 2)
			{
				printf("transcripts:\n");
				for (int i = 0; i < subg.trsts.size(); i++) subg.trsts[i].write(cout);
			}

			filter ft(subg.trsts);
			if (filter_tc.find('a') != string::npos) ft.join_single_exon_transcripts();
			ft.filter_length_coverage();
			if (ft.trs.size() >= 1) gv.insert(gv.end(), ft.trs.begin(), ft.trs.end());
			if (verbose >= 2)
			{
				printf("transcripts after filtering:\n");
				for (int i = 0; i < ft.trs.size(); i++) ft.trs[i].write(cout);
			}
		}

		if(fixed_gene_name != "" && gid == fixed_gene_name) terminate = true;
		if(terminate == true) return 0;
	}
	filter ft(gv);
	if (filter_tc.find('b') != string::npos) ft.remove_nested_transcripts();
	if (ft.trs.size() >= 1) trsts.insert(trsts.end(), ft.trs.begin(), ft.trs.end());
	fout.close();
	return 0;
}

bool assembler::determine_regional_graph(splice_graph &gr)
{
	bool all_regional = true;
	for(int i = 1; i < gr.num_vertices() - 1; i++)
	{
		if(gr.get_vertex_info(i).regional == false) all_regional = false;
		if(all_regional == false) break;
	}
	return all_regional;
}

int assembler::assign_RPKM()
{
	double factor = 1e9 / qlen;
	for(int i = 0; i < trsts.size(); i++)
	{
		trsts[i].assign_RPKM(factor);
	}
	return 0;
}

int assembler::write()
{
	ofstream fout(output_file.c_str());
	if (fout.fail()) exit(1);
	for(int i = 0; i < trsts.size(); i++)
	{
		transcript &t = trsts[i];
		t.write(fout);
	}
	fout.close();
	return 0;
}

int assembler::write_lp_stats(ostream &out)
{
	out << "lp_stats\t\tlp1\t\tlp2\n";
	out << "lp_cnt\t\t" << lp1_stats[0] << "\t\t" << lp2_stats[0] << "\n";
	out << "lp_merr\t\t" << lp1_stats[1] << "\t\t" << lp2_stats[1] << "\n";
	out << "lp_opt\t\t" << lp1_stats[2] << "\t\t" << lp2_stats[2] << "\n";
	out << "lp_tl_sol\t\t" << lp1_stats[3] << "\t\t" << lp2_stats[3] << "\n";
	out << "lp_tl_nsol\t\t" << lp1_stats[4] << "\t\t" << lp2_stats[4] << "\n";
	out << "lp_weird\t\t" << lp1_stats[5] << "\t\t" << lp2_stats[5] << "\n";
	out << "lp1merr_gids: " << lp1merr_gids << "\n";
	out << "lp1nsol_gids: " << lp1nsol_gids << "\n";
	out << "lp2merr_gids: " << lp2merr_gids << "\n";
	out << "lp2nsol_gids: " << lp2nsol_gids << "\n";
}

/*
int assembler::compare(splice_graph &gr, const string &file, const string &texfile)
{
	if(file == "") return 0;

	genome g(file);
	if(g.genes.size() <= 0) return 0;

	gtf gg(g.genes[0]);

	splice_graph gt;
	gg.build_splice_graph(gt);

	sgraph_compare sgc(gt, gr);
	sgc.compare(texfile);

	return 0;
}*/
