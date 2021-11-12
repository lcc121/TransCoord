/*
This file was modified from a file named config.cc in Scallop
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
#include "config.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <cstring>
#include <unistd.h>

using namespace std;

//// parameters
// for bam file and reads
int min_flank_length = 3;
int max_num_cigar = 7;
int max_edit_distance = 10;
int32_t min_bundle_gap = 50;
int min_num_hits_in_bundle = 20;
uint32_t min_mapping_quality = 1;
int32_t min_splice_boundary_hits = 1;
bool use_second_alignment = false;
bool uniquely_mapped_only = false;
int library_type = EMPTY;

// for preview
int max_preview_reads = 2000000;
int max_preview_spliced_reads = 50000;
int min_preview_spliced_reads = 10000;
double preview_infer_ratio = 0.95;
bool preview_only = false;

// for identifying subgraphs
int32_t min_subregion_gap = 3;
double min_subregion_overlap = 1.5;
int32_t min_subregion_length = 15;

// for revising/decomposing splice graph
double max_intron_contamination_coverage = 2.0;
double min_surviving_edge_weight = 1.5;
double max_decompose_error_ratio[7] = {0.33, 0.05, 0.0, 0.25, 0.30, 0.0, 1.1};

// for selecting paths
double min_transcript_coverage = 0;
double min_transcript_coverage_ratio = 0.005;
double min_single_exon_coverage = 100;
double min_single_node_coverage = 90;
double min_transcript_numreads = 20;
int min_transcript_length_base = 150;
int min_transcript_length_increase = 50;
int min_exon_length = 20;
int max_num_exons = 1000;

// for subsetsum and router
int max_dp_table_size = 10000;
int min_router_count = 1;

// for simulation
int simulation_num_vertices = 0;
int simulation_num_edges = 0;
int simulation_max_edge_weight = 0;

// input and output
string algo = "transcoord";
string input_file;
string ref_file;
string ref_file1;
string ref_file2;
string output_file;
string output_prefix;
string log_file;     
string stats_file;   

// for controling
bool output_tex_files = false;
string fixed_gene_name = "";
int fixed_bundle_idx = -1;
//int batch_bundle_size = 100;
int verbose = 0;              
string version = "v0.0.0";

// for subgraph finding paths
string path_options = "de";          
string filter_path = "b";           
int small_path_mode = 1;            
int small_path_num = 12;
int simple_path_num = 12;
bool output_all_paths = false;    
bool use_backup_paths = false; 
double min_cov = 0.1;

// for lp
int lp_option = 4;                  
int abs_tol_option = 1;             
int max_AB_element = 70000000;
int thread_count = 1;
double tlb = 0.0025;                
string glog_file1;
string glog_file2;
string glp_file1;
string glp_file2;
string lp1merr_gids;                
string lp2merr_gids;
string lp1nsol_gids;                
string lp2nsol_gids;
int lp1_stats[6] = { 0,0,0,0,0,0 }; 
int lp2_stats[6] = { 0,0,0,0,0,0 };

// for transcript
string filter_tc = "aceh";      

int parse_arguments(int argc, const char ** argv)
{
	for(int i = 1; i < argc; i++)
	{
		// necessary ones
		if(string(argv[i]) == "-i")
		{
			input_file = string(argv[i + 1]);
			if (verbose >= 3) cout << "input_file=" << input_file << endl;
			i++;
		}
		else if(string(argv[i]) == "-o")
		{
			output_file = string(argv[i + 1]);
			if (verbose >= 3) cout << "output_file=" << output_file << endl;

			const size_t last_slash_idx = output_file.rfind('/');
			if (std::string::npos != last_slash_idx) {
				output_prefix = output_file.substr(0, last_slash_idx);
				log_file = output_prefix + "/bundle.txt";
				glog_file1 = output_prefix + "/gurobi_mip1.log";
				glog_file2 = output_prefix + "/gurobi_mip2.log";
				glp_file1 = output_prefix + "/gurobi_mip1.lp";
				glp_file2 = output_prefix + "/gurobi_mip2.lp";
				stats_file = output_prefix + "/stats.txt";
			}
			if (access(output_prefix.c_str(), 0) == -1) {
				//mkdir(prefix.c_str());
				printf("error: directory %s not exit!", output_prefix.c_str());
				exit(1);
			}
			i++;
		}

		// internal use
		else if(string(argv[i]) == "-a")
		{
			algo = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r")
		{
			ref_file = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r1")
		{
			ref_file1 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-r2")
		{
			ref_file2 = string(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-g")
		{
			fixed_gene_name = string(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--fixed_bundle_idx")
		{
			fixed_bundle_idx = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "-t")
		{
			output_tex_files = true;
		}

		// user specified
		else if(string(argv[i]) == "--version")
		{
			printf("%s\n", version.c_str());
			exit(0);
		}
		else if(string(argv[i]) == "--help")
		{
			print_copyright();
			print_help();
			printf("\n");
			//print_logo();
			exit(0);
		}
		else if(string(argv[i]) == "--min_flank_length")
		{
			min_flank_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_cigar")
		{
			max_num_cigar = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_edit_distance")
		{
			max_edit_distance = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_bundle_gap")
		{
			min_bundle_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_num_hits_in_bundle")
		{
			min_num_hits_in_bundle = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_mapping_quality")
		{
			min_mapping_quality = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_splice_boundary_hits")
		{
			min_splice_boundary_hits = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_preview_spliced_reads")
		{
			max_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_preview_spliced_reads")
		{
			min_preview_spliced_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview")
		{
			preview_only = true;
		}
		else if(string(argv[i]) == "--max_preview_reads")
		{
			max_preview_reads = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--preview_infer_ratio")
		{
			preview_infer_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_gap")
		{
			min_subregion_gap = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_length")
		{
			min_subregion_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_subregion_overlap")
		{
			min_subregion_overlap = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_surviving_edge_weight")
		{
			min_surviving_edge_weight = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_intron_contamination_coverage")
		{
			max_intron_contamination_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_coverage")
		{
			min_transcript_coverage = atof(argv[i + 1]);
			i++;
			//if(fabs(min_transcript_coverage - 1.0) < 0.01) min_transcript_coverage = 1.01;
		}
		else if(string(argv[i]) == "--min_transcript_coverage_ratio")
		{
			min_transcript_coverage_ratio = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_single_exon_coverage")
		{
			min_single_exon_coverage = atof(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--min_single_node_coverage")
		{
			min_single_node_coverage = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_numreads")
		{
			min_transcript_numreads = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_base")
		{
			min_transcript_length_base = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_transcript_length_increase")
		{
			min_transcript_length_increase = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_exon_length")
		{
			min_exon_length = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_num_exons")
		{
			max_num_exons = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_dp_table_size")
		{
			max_dp_table_size = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--min_router_count")
		{
			min_router_count = atoi(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio0")
		{
			max_decompose_error_ratio[0] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio1")
		{
			max_decompose_error_ratio[1] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio2")
		{
			max_decompose_error_ratio[2] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio3")
		{
			max_decompose_error_ratio[3] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio4")
		{
			max_decompose_error_ratio[4] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio5")
		{
			max_decompose_error_ratio[5] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--max_decompose_error_ratio6")
		{
			max_decompose_error_ratio[6] = atof(argv[i + 1]);
			i++;
		}
		else if(string(argv[i]) == "--library_type")
		{
			string s(argv[i + 1]);
			if(s == "empty") library_type = EMPTY;
			if(s == "unstranded") library_type = UNSTRANDED;
			if(s == "first") library_type = FR_FIRST;
			if(s == "second") library_type = FR_SECOND;
			i++;
		}
		else if(string(argv[i]) == "--use_second_alignment")
		{
			string s(argv[i + 1]);
			if(s == "true") use_second_alignment = true;
			else use_second_alignment = false;
			i++;
		}
		else if(string(argv[i]) == "--uniquely_mapped_only")
		{
			string s(argv[i + 1]);
			if(s == "true") uniquely_mapped_only = true;
			else uniquely_mapped_only = false;
			i++;
		}
		else if(string(argv[i]) == "--verbose")
		{
			verbose = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--max_remain_paths")
		{
			max_AB_element = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--path_options")
		{
			path_options = string(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--lp_option")
		{
			lp_option = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--abs_tol_option")
		{
			abs_tol_option = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--tlb")
		{
			tlb = atof(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--small_path_mode")
		{
			small_path_mode = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--small_path_num")
		{
			small_path_num = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--simple_path_num")
		{
			simple_path_num = atoi(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--output_all_paths")
		{
			string s(argv[i + 1]);
			if (s == "true") output_all_paths = true;
			else output_all_paths = false;
			i++;
		}
		else if (string(argv[i]) == "--use_backup_paths")
		{
			string s(argv[i + 1]);
			if (s == "true") use_backup_paths = true;
			else use_backup_paths = false;
			i++;
		}
		else if (string(argv[i]) == "--filter_tc")
		{
			filter_tc = string(argv[i + 1]);
			i++;
		}
		else if (string(argv[i]) == "--filter_path")
		{
			filter_path = string(argv[i + 1]);
			i++;
		}
	}

	if(min_surviving_edge_weight < 0.1 + min_transcript_coverage) 
	{
		min_surviving_edge_weight = 0.1 + min_transcript_coverage;
	}

	// verify arguments
	if(input_file == "")
	{
		printf("error: input-file is missing.\n");
		exit(0);
	}

	if(output_file == "" && preview_only == false)
	{
		printf("error: output-file is missing.\n");
		exit(0);
	}

	return 0;
}

int print_parameters()
{
	printf("parameters:\n");

	// for bam file and reads
	printf("min_flank_length = %d\n", min_flank_length);
	printf("max_num_cigar = %d\n", max_num_cigar);
	printf("max_edit_distance = %d\n", max_edit_distance);
	printf("min_bundle_gap = %d\n", min_bundle_gap);
	printf("min_num_hits_in_bundle = %d\n", min_num_hits_in_bundle);
	printf("min_mapping_quality = %d\n", min_mapping_quality);
	printf("min_splice_boundary_hits = %d\n", min_splice_boundary_hits);

	// for preview
	printf("preview_only = %c\n", preview_only ? 'T' : 'F');
	printf("max_preview_reads = %d\n", max_preview_reads);
	printf("max_preview_spliced_reads = %d\n", max_preview_spliced_reads);
	printf("min_preview_spliced_reads = %d\n", min_preview_spliced_reads);
	printf("preview_infer_ratio = %.3lf\n", preview_infer_ratio);

	// for identifying subgraphs
	printf("min_subregion_gap = %d\n", min_subregion_gap);
	printf("min_subregion_length = %d\n", min_subregion_length);
	printf("min_subregion_overlap = %.2lf\n", min_subregion_overlap);

	// for splice graph
	printf("max_intron_contamination_coverage = %.2lf\n", max_intron_contamination_coverage);
	printf("min_surviving_edge_weight = %.2lf\n", min_surviving_edge_weight);
	printf("min_transcript_coverage = %.2lf\n", min_transcript_coverage);
	printf("min_transcript_coverage_ratio = %.2lf\n", min_transcript_coverage_ratio);
	printf("min_single_exon_coverage = %.2lf\n", min_single_exon_coverage);
	printf("min_transcript_numreads = %.2lf\n", min_transcript_numreads);
	printf("min_transcript_length_base = %d\n", min_transcript_length_base);
	printf("min_transcript_length_increase = %d\n", min_transcript_length_increase);
	printf("max_num_exons = %d\n", max_num_exons);

	// for input and output
	printf("algo = %s\n", algo.c_str());
	printf("input_file = %s\n", input_file.c_str());
	printf("ref_file = %s\n", ref_file.c_str());
	printf("ref_file1 = %s\n", ref_file1.c_str());
	printf("ref_file2 = %s\n", ref_file2.c_str());
	printf("output_file = %s\n", output_file.c_str());

	// for controling
	printf("library_type = %d\n", library_type);
	printf("output_tex_files = %c\n", output_tex_files ? 'T' : 'F');
	printf("fixed_gene_name = %s\n", fixed_gene_name.c_str());
	printf("use_second_alignment = %c\n", use_second_alignment ? 'T' : 'F');
	printf("uniquely_mapped_only = %c\n", uniquely_mapped_only ? 'T' : 'F');
	printf("verbose = %d\n", verbose);

	printf("\n");

	return 0;
}

int print_command_line(int argc, const char ** argv)
{
	printf("command line: ");
	for(int i = 0; i < argc; i++)
	{
		printf("%s ", argv[i]);
	}
	printf("\n");
	return 0;
}

int print_help()
{
	printf("\n");
	printf("Usage: transcoord -i <bam-file> -o <gtf-file> [options]\n");
	printf("\n");
	printf("Options:\n");
	printf(" %-42s  %s\n", "--help",  "print usage of TransCoord and exit");
	printf(" %-42s  %s\n", "--version",  "print current version of TransCoord and exit");
	printf(" %-42s  %s\n", "--verbose <0, 1, 2>",  "0: quiet; 1: one line for each graph; 2: with details, default: 1");
	printf(" %-42s  %s\n", "--library_type <first, second, unstranded>",  "library type of the sample, default: unstranded");
	printf(" %-42s  %s\n", "--min_transcript_coverage <float>",  "minimum coverage required for a multi-exon transcript, default: 0");
	printf(" %-42s  %s\n", "--min_single_exon_coverage <float>",  "minimum coverage required for a single-exon transcript, default: 100");
	printf(" %-42s  %s\n", "--min_mapping_quality <integer>",  "ignore reads with mapping quality less than this value, default: 1");
	printf(" %-42s  %s\n", "--max_num_cigar <integer>",  "ignore reads with CIGAR size larger than this value, default: 7");
	printf(" %-42s  %s\n", "--min_bundle_gap <integer>",  "minimum distances required to start a new bundle, default: 50");
	printf(" %-42s  %s\n", "--min_num_hits_in_bundle <integer>",  "minimum number of reads required in a bundle, default: 20");
	printf(" %-42s  %s\n", "--min_flank_length <integer>",  "minimum match length in each side for a spliced read, default: 3");
	printf(" %-42s  %s\n", "--min_splice_bundary_hits <integer>",  "minimum number of spliced reads required for a junction, default: 1");
	return 0;
}

int print_copyright()
{
	printf("TransCoord %s \n", version.c_str());
	return 0;
}
