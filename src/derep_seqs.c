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

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <zlib.h>

#include "vendor/kseq.h"
#include "vendor/ssef.h"
#include "vendor/tommyarray.h"
#include "pointers.h"
#include "rseq.h"
#include "kseq_helper.h"

#define VERSION "0.1.0"

int
main(int argc, char *argv[])
{
  if (argc < 2) {
    fprintf(stderr,
            "Version: %s --- Usage: %s <contigs.fasta>\n",
            VERSION,
            argv[0]);

    return 1;
  }

  long l = 0;
  kseq_t* seq;

  rseq_t* rseq;
  rseq_t* search_term;
  rseq_t* text;

  gzFile fp;

  tommy_size_t num_seqs = 0;
  tommy_size_t i = 0;
  tommy_size_t j = 0;

  unsigned long seqs_too_short = 0;
  unsigned long seqs_not_printed = 0;
  unsigned long seqs_printed = 0;

  fp = gzopen(argv[1], "r");
  assert(fp);

  seq = kseq_init(fp);
  assert(seq);

  tommy_array* seqs = malloc(sizeof(tommy_array));
  assert(seqs != NULL);
  tommy_array_init(seqs);

  int result = 0;

  fprintf(stderr, "LOG -- reading seqs into memory\n");
  while ((l = kseq_read(seq)) >= 0) {
    rseq = rseq_init(seq);

    tommy_array_insert(seqs, rseq);
  }
  num_seqs = tommy_array_size(seqs);

  fprintf(stderr, "LOG -- read %lu seqs\n", num_seqs);

  fprintf(stderr, "LOG -- dereplicating seqs\n");
  for (i = 0; i < num_seqs - 1; ++i) {
    if ((i % 10000) == 0) {
      fprintf(stderr,
              "LOG -- searching for seq: %lu of %lu (%.2f%%)\r",
              i+1,
              num_seqs,
              (i+1) / (double)num_seqs * 100);
    }
    result = -2;
    search_term = tommy_array_get(seqs, i);
    for (j = i + 1; j < num_seqs; ++j) {
      text = tommy_array_get(seqs, j);

      result = ssef_search((unsigned char*)search_term->seq,
                           search_term->len,
                           (unsigned char*)text->seq,
                           text->len);

      if (result == -1) {
        fprintf(stderr,
                "WARN -- The sequence '%s' to search was too short! "
                "Printing it anyway, but it may be a substring of another seq!\n",
                search_term->head);
        rseq_print(stdout, search_term);
        ++seqs_too_short;
        break;
      } else if (result > 0) {
        ++seqs_not_printed;
        break; /* query i has a match with target j, stop searching for i */
      }
    }

    if (result == 0) { /* seq i is unique */
      ++seqs_printed;
      rseq_print(stdout, search_term);
    }
  }

  /* print out the last seq (longest one can't be a substring!) */
  search_term = tommy_array_get(seqs, i);
  rseq_print(stdout, search_term);
  ++seqs_printed;

  fprintf(stderr,
          "LOG -- searching for seq: %lu of %lu (%.2f%%)\n",
          i+1,
          num_seqs,
          (i+1) / (double)num_seqs * 100);


  double total = seqs_printed + seqs_not_printed;

  fprintf(stderr,
          "LOG -- seqs printed: %lu (%.2f%%)\n",
          seqs_printed,
          seqs_printed / total * 100);
  fprintf(stderr,
          "LOG -- seqs not printed: %lu (%.2f%%)\n",
          seqs_not_printed,
          seqs_not_printed / total * 100);
  fprintf(stderr,
          "LOG -- seqs too short (and printed): %lu (%.2f%%)\n",
          seqs_too_short,
          seqs_too_short / total * 100);

  fprintf(stderr, "LOG -- freeing memory\n");

  kseq_destroy(seq);
  gzclose(fp);

  for (i = 0; i < num_seqs; ++i) {
    rseq_destroy(tommy_array_get(seqs, i));
  }
  tommy_array_done(seqs);
  free(seqs);

  fprintf(stderr, "LOG -- done\n");

  return 0;
}
