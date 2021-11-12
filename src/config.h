/*
This file was modified from a file named config.h in Scallop
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
#ifndef __CONFIG_H__
#define __CONFIG_H__

#include "util.h"
#include <stdint.h>
#include <map>
#include <sstream>

using namespace std;

//// constants
#define START_BOUNDARY 1
#define END_BOUNDARY 2
#define LEFT_SPLICE 3 
#define RIGHT_SPLICE 4
#define LEFT_RIGHT_SPLICE 5
#define MIDDLE_CUT 6

#define TRIVIAL 0
#define NORMAL 1

// five types for decomposition
#define SMALLEST_EDGE 0
#define NEGLIGIBLE_EDGE 1
#define SPLITTABLE_SIMPLE 2
#define SPLITTABLE_HYPER 3
#define UNSPLITTABLE_SINGLE 4
#define UNSPLITTABLE_MULTIPLE 5
#define TRIVIAL_VERTEX 6

#define EMPTY -1
#define UNSTRANDED 0
#define FR_FIRST 1
#define FR_SECOND 2

//// parameters
// for bam file and reads
extern int min_flank_length;
extern int max_num_cigar;
extern int max_edit_distance;
extern int32_t min_bundle_gap;
extern int min_num_hits_in_bundle;
extern uint32_t min_mapping_quality;
extern int32_t min_splice_boundary_hits;
extern bool uniquely_mapped_only;
extern bool use_second_alignment;

// for preview
extern bool preview_only;
extern int max_preview_reads;
extern int max_preview_spliced_reads;
extern int min_preview_spliced_reads;
extern double preview_infer_ratio;

// for identifying subgraphs
extern int32_t min_subregion_gap;
extern double min_subregion_overlap;
extern int32_t min_subregion_length;
extern int min_subregion_ladders;

// for subsetsum and router
extern int max_dp_table_size;
extern int min_router_count;

// for splice graph
extern double max_intron_contamination_coverage;
extern double min_surviving_edge_weight;
extern double max_decompose_error_ratio[7];
extern double min_transcript_numreads;
extern double min_transcript_coverage;
extern double min_single_exon_coverage;
extern double min_single_node_coverage;
extern double min_transcript_coverage_ratio; 
extern int min_transcript_length_base;
extern int min_transcript_length_increase;
extern int min_exon_length;
extern int max_num_exons;

// for simulation
extern int simulation_num_vertices;
extern int simulation_num_edges;
extern int simulation_max_edge_weight;

// input and output
extern string algo;
extern string input_file;
extern string ref_file;
extern string ref_file1;
extern string ref_file2;
extern string output_file;
extern string output_prefix;
extern string log_file;
extern string stats_file;

// for controling
extern bool output_tex_files;
extern string fixed_gene_name;
extern int fixed_bundle_idx;
extern int max_num_bundles;
extern int library_type;
extern int min_gtf_transcripts_num;
//extern int batch_bundle_size;
extern int verbose;
extern string version;

// for subgraph finding paths
extern string path_options;
extern string filter_path;
extern int small_path_mode;
extern int small_path_num;
extern int simple_path_num;
extern bool output_all_paths;
extern bool use_backup_paths;
extern double min_cov;

// for lp
extern int lp_option;
extern int abs_tol_option;
extern int max_AB_element;
extern int thread_count;
extern double tlb;
extern string glog_file1;
extern string glog_file2;
extern string glp_file1;
extern string glp_file2;
extern string lp1merr_gids;
extern string lp2merr_gids;
extern string lp1nsol_gids;
extern string lp2nsol_gids;
extern int lp1_stats[6];
extern int lp2_stats[6];

extern string filter_tc;

// parse arguments
int print_command_line(int argc, const char ** argv);
int parse_arguments(int argc, const char ** argv);
int print_parameters();
int print_copyright();
int print_help();

#endif
