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
#include "vendor/hash3.h"
#include "vendor/tommyarray.h"
#include "vendor/tommyhash.h"
#include "pointers.h"
#include "rseq.h"
#include "kseq_helper.h"

#define VERSION "0.4.0"
#define WORD_SIZE 20

/* from http://mgronhol.github.io/fast-strcmp/ */
int fast_compare( const char *ptr0, const char *ptr1, int len ){
  int fast = len/sizeof(size_t) + 1;
  int offset = (fast-1)*sizeof(size_t);
  int current_block = 0;

  if( len <= sizeof(size_t)){ fast = 0; }


  size_t *lptr0 = (size_t*)ptr0;
  size_t *lptr1 = (size_t*)ptr1;

  while( current_block < fast ){
    if( (lptr0[current_block] ^ lptr1[current_block] )){
      int pos;
      for(pos = current_block*sizeof(size_t); pos < len ; ++pos ){
        if( (ptr0[pos] ^ ptr1[pos]) || (ptr0[pos] == 0) || (ptr1[pos] == 0) ){
          return  (int)((unsigned char)ptr0[pos] - (unsigned char)ptr1[pos]);
        }
      }
    }

    ++current_block;
  }

  while( len > offset ){
    if( (ptr0[offset] ^ ptr1[offset] )){
      return (int)((unsigned char)ptr0[offset] - (unsigned char)ptr1[offset]);
    }
    ++offset;
  }


  return 0;
}

struct str_t {
  char* str;
  tommy_node node;
};

struct str_t*
str_dup(struct str_t* str)
{
  struct str_t* new_str = malloc(sizeof(struct str_t));
  assert(new_str != NULL);
  new_str->str = strdup(str->str);

  return new_str;
}

void
str_destroy(struct str_t* str)
{
  free(str->str);
  free(str);
}

int
str_compare(const void* arg, const void* str)
{
  return strcmp((const char*)arg,
                ((const struct str_t*)str)->str);
}

struct hashed_kmer_t {
  tommy_uint32_t hashed_kmer;
  tommy_node node;
};

int
hashed_kmer_compare(const void* arg, const void* hk)
{
  return *(const tommy_uint32_t*)arg !=
    ((const struct hashed_kmer_t*)hk)->hashed_kmer;
}

struct hashed_kmer_t*
hashed_kmer_init(tommy_uint32_t hashed_val)
{
  struct hashed_kmer_t* hashed_kmer = malloc(sizeof(struct hashed_kmer_t));
  assert(hashed_kmer != NULL);

  hashed_kmer->hashed_kmer = hashed_val;

  return hashed_kmer;
}

void
hashed_kmer_destroy(struct hashed_kmer_t* hashed_kmer)
{
  free(hashed_kmer);
}

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

  unsigned long num_kmers = 0;
  unsigned int kmer_pos_i = 0;

  char* kmer = malloc((WORD_SIZE + 1) * sizeof(char));
  assert(kmer != NULL);

  tommy_uint32_t hash_val = 0;

  unsigned long kmer_i = 0;

  unsigned long strcmp_calls = 0;

  unsigned long passed_filter = 0;

  unsigned long idx = 0;
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

  struct hashed_kmer_t* hashed_kmer;
  struct hashed_kmer_t* hashed_kmer_ret;
  /* tommy_array* hashed_kmers; */
  tommy_hashlin* counting_hash = NULL;

  int result = 0;

  fprintf(stderr, "LOG -- reading seqs into memory\n");
  while ((l = kseq_read(seq)) >= 0) {
    if ((++idx % 1000) == 0) {
      fprintf(stderr, "%lu\r", idx);
    }
    rseq = rseq_init(seq);


    /* TODO print out too short seqs here */
    tommy_array_insert(seqs, rseq);
  }
  num_seqs = tommy_array_size(seqs);

  int unique_seqs[num_seqs];

  fprintf(stderr, "LOG -- read %lu seqs\n", num_seqs);

  fprintf(stderr, "LOG -- dereplicating seqs\n");

  /* first one is longest, always unique */
  rseq_print(stdout, tommy_array_get(seqs, 0));
  unique_seqs[0] = 1;
  ++seqs_printed;

  idx = 0;

  for (i = 1; i < num_seqs; ++i) {
    result = -2;
    search_term = tommy_array_get(seqs, i);

    if (i % 100 == 0) {
      fprintf(stderr,
              "comparing seq: %lu of %lu (%.2f%%)\r",
              i,
              num_seqs,
              (i) / (double)num_seqs * 100);
    }

    if (search_term->len < 33 ) { /* 32 is hard limit for ssef */
      rseq_print(stdout, search_term);
      ++seqs_too_short;
      ++seqs_printed;
      unique_seqs[i] = 1;

    } else {

      for (j = 0; j < i; ++j) {

        if (unique_seqs[j] == 1) { /* only check if j hasn't already
                                      been eliminated */
          text = tommy_array_get(seqs, j);

          /* TODO is strcmp or hashing faster than this when lengths are
             equal? */
          if (text->len == search_term->len) {
            ++strcmp_calls;
            result = fast_compare(text->seq, search_term->seq, text->len);
            if (result == 0) { /* match */
              result = 1;
            } else {
              result = 0;
            }
          } else {
            /* if the first 33 (32 is min search term for SSEF) bases in
               search_term don't have a match, the whole thing wont have
               a match */

            /* kmers = eigg_kmers_new(text->seq, text->len, WORD_SIZE); */

            counting_hash = malloc(sizeof(tommy_hashlin));
            assert(counting_hash != NULL);
            tommy_hashlin_init(counting_hash);

            /* hash the kmers of the text */
            num_kmers = text->len - WORD_SIZE + 1;
            for (kmer_i = 0; kmer_i < num_kmers; ++kmer_i) {
              /* unroll this loop for speed */
              for (kmer_pos_i = 0; kmer_pos_i < WORD_SIZE; kmer_pos_i += 10) {
                kmer[kmer_pos_i  ] = text->seq[kmer_i + kmer_pos_i    ];
                kmer[kmer_pos_i+1] = text->seq[kmer_i + kmer_pos_i + 1];
                kmer[kmer_pos_i+2] = text->seq[kmer_i + kmer_pos_i + 2];
                kmer[kmer_pos_i+3] = text->seq[kmer_i + kmer_pos_i + 3];
                kmer[kmer_pos_i+4] = text->seq[kmer_i + kmer_pos_i + 4];
                kmer[kmer_pos_i+5] = text->seq[kmer_i + kmer_pos_i + 5];
                kmer[kmer_pos_i+6] = text->seq[kmer_i + kmer_pos_i + 6];
                kmer[kmer_pos_i+7] = text->seq[kmer_i + kmer_pos_i + 7];
                kmer[kmer_pos_i+8] = text->seq[kmer_i + kmer_pos_i + 8];
                kmer[kmer_pos_i+9] = text->seq[kmer_i + kmer_pos_i + 9];
              }
              kmer[kmer_pos_i] = '\0';

              hash_val = tommy_hash_u32(0, kmer, WORD_SIZE);
              hashed_kmer_ret = tommy_hashlin_search(counting_hash,
                                                     hashed_kmer_compare,
                                                     &hash_val,
                                                     hash_val);

              if (!hashed_kmer_ret) {
                hashed_kmer = hashed_kmer_init(hash_val);
                tommy_hashlin_insert(counting_hash,
                                     &hashed_kmer->node,
                                     hashed_kmer,
                                     hashed_kmer->hashed_kmer);
              }
            }

            /* TODO should we save these hash vals */
            /* eigg_kmers_destroy(kmers); */

            /* check the kmers of the search_term against those in
               counting hash */
            /* kmers = eigg_kmers_new(search_term->seq, search_term->len, WORD_SIZE); */

            num_kmers = search_term->len - WORD_SIZE + 1;
            for (kmer_i = 0; kmer_i < num_kmers; ++kmer_i) {
              for (kmer_pos_i = 0;
                   kmer_pos_i < WORD_SIZE;
                   kmer_pos_i += 10) {
                kmer[kmer_pos_i  ] = search_term->seq[kmer_i + kmer_pos_i    ];
                kmer[kmer_pos_i+1] = search_term->seq[kmer_i + kmer_pos_i + 1];
                kmer[kmer_pos_i+2] = search_term->seq[kmer_i + kmer_pos_i + 2];
                kmer[kmer_pos_i+3] = search_term->seq[kmer_i + kmer_pos_i + 3];
                kmer[kmer_pos_i+4] = search_term->seq[kmer_i + kmer_pos_i + 4];
                kmer[kmer_pos_i+5] = search_term->seq[kmer_i + kmer_pos_i + 5];
                kmer[kmer_pos_i+6] = search_term->seq[kmer_i + kmer_pos_i + 6];
                kmer[kmer_pos_i+7] = search_term->seq[kmer_i + kmer_pos_i + 7];
                kmer[kmer_pos_i+8] = search_term->seq[kmer_i + kmer_pos_i + 8];
                kmer[kmer_pos_i+9] = search_term->seq[kmer_i + kmer_pos_i + 9];

              }
              kmer[kmer_pos_i] = '\0';

              hash_val = tommy_hash_u32(0, kmer, WORD_SIZE);
              /* check if this kmer is in the hash table */
              hashed_kmer_ret = tommy_hashlin_search(counting_hash,
                                                     hashed_kmer_compare,
                                                     &hash_val,
                                                     hash_val);

              /* TODO only screen some of the kmers? */
              if (hashed_kmer_ret) {
                /* fprintf(stderr, "match! %s, %u\n", kmer, hash_val); */
                result = 1;
              } else {
                result = 0;
                break;
              }
            }
            /* eigg_kmers_destroy(kmers); */
            /* if i checked actual kmers instead of hash vals and they
               all matched, I could avoid this check */
            if (result > 0) {
              ++passed_filter;
              /* passed the first filter, test the whole pattern */
              result = hash3_search((unsigned char*)search_term->seq,
                                    search_term->len,
                                    (unsigned char*)text->seq,
                                    text->len);
            }
            /* clear this counting hash */
            tommy_hashlin_foreach(counting_hash, free);
            tommy_hashlin_done(counting_hash);
            free(counting_hash);

          }
          /* assert(result != -1);  */

          if (result > 0) { /* there is a match */
            unique_seqs[i] = 0;
            ++seqs_not_printed;
            break; /* query i has a match with target j, stop searching for i */
          }

        } else {
          continue;
        }
      }

      if (result == 0) { /* seq i is unique */
        rseq_print(stdout, search_term);
        unique_seqs[i] = 1;
        ++seqs_printed;
      }
    }
  }
  fprintf(stderr,
          "                                               "
          "                                               \r");

  fprintf(stderr,
          "LOG -- seqs printed: %lu (%.2f%%)\n",
          seqs_printed,
          seqs_printed / (double)num_seqs * 100);
  fprintf(stderr,
          "LOG -- seqs not printed: %lu (%.2f%%)\n",
          seqs_not_printed,
          seqs_not_printed / (double)num_seqs * 100);
  fprintf(stderr,
          "LOG -- seqs too short (and printed): %lu (%.2f%%)\n",
          seqs_too_short,
          seqs_too_short / (double)num_seqs * 100);
  fprintf(stderr,
          "LOG -- strcmp calls: %lu\n",
          strcmp_calls);
  fprintf(stderr,
          "LOG -- passed filter (size %d): %lu\n",
          WORD_SIZE,
          passed_filter);

  fprintf(stderr, "LOG -- freeing memory\n");

  kseq_destroy(seq);
  gzclose(fp);

  for (i = 0; i < num_seqs; ++i) {
    rseq_destroy(tommy_array_get(seqs, i),
                 (tommy_foreach_func*)hashed_kmer_destroy);
  }
  tommy_array_done(seqs);
  free(seqs);

  fprintf(stderr, "LOG -- done\n");

  free(kmer);

  return 0;
}
