/*
This file was modified from a file named vertex_base.cc in Scallop
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

#include "vertex_base.h"
#include "edge_base.h"

#include <cstdio>
#include <cassert>
#include <cstdio>

using namespace std;

vertex_base::vertex_base()
{}

vertex_base::~vertex_base()
{}

int vertex_base::add_in_edge(edge_base *e)
{
	assert(si.find(e) == si.end());
	si.insert(e);
	return 0;
}

int vertex_base::add_out_edge(edge_base *e)
{
	assert(so.find(e) == so.end());
	so.insert(e);
	return 0;
}

int vertex_base::remove_in_edge(edge_base *e)
{
	assert(si.find(e) != si.end());
	si.erase(e);
	return 0;
}

int vertex_base::remove_out_edge(edge_base *e)
{
	assert(so.find(e) != so.end());
	so.erase(e);
	return 0;
}

int vertex_base::degree() const
{
	return in_degree() + out_degree();
}

int vertex_base::in_degree() const
{
	return (int)(si.size());
}

int vertex_base::out_degree() const
{
	return (int)(so.size());
}

PEEI vertex_base::in_edges() const
{
	return PEEI(si.begin(), si.end());
}

PEEI vertex_base::out_edges() const
{
	return PEEI(so.begin(), so.end());
}

int vertex_base::print() const
{
	printf("in-edges = ( ");
	for(edge_iterator it = si.begin(); it != si.end(); it++)
	{
		printf("[%d, %d] ", (*it)->source(), (*it)->target());
	}
	printf("), out-edges = ( ");
	for(edge_iterator it = so.begin(); it != so.end(); it++)
	{
		printf("[%d, %d] ", (*it)->source(), (*it)->target());
	}
	printf(")\n");
	return 0;
}
