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

#include "rseq.h"

rseq_t* rseq_init(kseq_t* kseq)
{
  char buf[1000];
  rseq_t* rseq = malloc(sizeof(rseq_t));
  assert(rseq != NULL);

  if (kseq->comment.l) {
    sprintf(buf,
            "%s %s",
            kseq->name.s,
            kseq->comment.s);
    rseq->head = strdup(buf);
  } else {
    rseq->head = strdup(kseq->name.s);
  }

  rseq->seq = strdup(kseq->seq.s);
  rseq->len = kseq->seq.l;

  rseq->kmers = NULL;
  rseq->hashed_kmers = NULL;
  rseq->unique_kmers = 0;

  return rseq;
}

void
rseq_destroy(rseq_t* rseq,
             tommy_foreach_func* rseq_free_func)
{
  free(rseq->head);
  free(rseq->seq);

  if (rseq->kmers != NULL) {
    tommy_hashlin_foreach(rseq->kmers, rseq_free_func);
    tommy_hashlin_done(rseq->kmers);
    free(rseq->kmers);
  }

  if (rseq->hashed_kmers != NULL) {
    tommy_array_done(rseq->hashed_kmers);
    free(rseq->hashed_kmers);
  }

  free(rseq);
}

void
rseq_print(FILE* fstream, rseq_t* rseq)
{
  fprintf(fstream,
          ">%s\n"
          "%s\n",
          rseq->head,
          rseq->seq);
}
