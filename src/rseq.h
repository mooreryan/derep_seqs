/*
Copyright 2017 Ryan Moore
Contact: moorer@udel.edu

This file is part of derep_seqs.

derep_seqs is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

derep_seqs is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef RSEQ_H
#define RSEQ_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>
#include "kseq_helper.h"
#include "vendor/tommyhashlin.h"
#include "vendor/tommyarray.h"

typedef struct rseq_t {
  char* head;
  char* seq;
  unsigned long len;
  tommy_hashlin* kmers;
  tommy_array* hashed_kmers;
  unsigned long unique_kmers;
} rseq_t;

rseq_t* rseq_init(kseq_t* kseq);

void
rseq_destroy(rseq_t* rseq, tommy_foreach_func* free_func);

void
rseq_print(FILE* fstream, rseq_t* rseq);

#endif
