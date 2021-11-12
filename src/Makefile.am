bin_PROGRAMS = transcoord


GTF_INCLUDE = $(top_srcdir)/lib/gtf
UTIL_INCLUDE = $(top_srcdir)/lib/util
GRAPH_INCLUDE = $(top_srcdir)/lib/graph

GTF_LIB = $(top_builddir)/lib/gtf
UTIL_LIB = $(top_builddir)/lib/util
GRAPH_LIB = $(top_builddir)/lib/graph

#transcoord_CPPFLAGS = -std=c++11 -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
transcoord_CPPFLAGS = -I$(GTF_INCLUDE) -I$(GRAPH_INCLUDE) -I$(UTIL_INCLUDE)
transcoord_LDFLAGS = -pthread -L$(GTF_LIB) -L$(GRAPH_LIB) -L$(UTIL_LIB)
transcoord_LDADD = -lgtf -lgraph -lutil

transcoord_SOURCES = splice_graph.h splice_graph.cc \
				     super_graph.h super_graph.cc \
				     vertex_info.h vertex_info.cc \
				     edge_info.h edge_info.cc \
				     interval_map.h interval_map.cc \
				     config.h config.cc \
				     hit.h hit.cc \
				     partial_exon.h partial_exon.cc \
				     hyper_set.h hyper_set.cc \
				     region.h region.cc \
				     junction.h junction.cc \
				     bundle.h bundle.cc \
				     previewer.h previewer.cc \
				     assembler.h assembler.cc \
				     filter.h filter.cc \
				     main.cc \
				     subgraph.h subgraph.cc \
                     path.h path.cc \
                     gurobi_mip.h gurobi_mip.cc \
				     scallop.h scallop.cc \
                     equation.h equation.cc \
                     router.h router.cc \
                     subsetsum.h subsetsum.cc \
                     overlap_graph.h overlap_graph.cc \
                     phasing_graph.h phasing_graph.cc
