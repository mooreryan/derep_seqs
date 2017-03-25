#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <zlib.h>

#include "rabin_karp.h"
#include "vendor/hash3.h"
#include "vendor/kseq.h"
#include "vendor/tommyarray.h"
#include "vendor/tommyhashlin.h"

#define VERSION "0.7.0"
pthread_mutex_t mutex;
unsigned long passed_prefilter = 0;
unsigned long non_unique_seqs = 0;
int* non_unique_seq_index;
tommy_array* patterns = NULL;

KSEQ_INIT(gzFile, gzread)

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

void apple( long cycles)
{
  for ( long i = 0; i < cycles; ++i) {
    ;
  }
}

typedef struct rseq_t {
  char* head;
  char* seq;
  unsigned long len;
} rseq_t;

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

  return rseq;
}

void
rseq_destroy(rseq_t* rseq)
{
  free(rseq->head);
  free(rseq->seq);
  free(rseq);
}

void
rseq_print(FILE* fstream, rseq_t* rseq)
{
  /* TODO does this need a lock? */
  fprintf(fstream,
          ">%s\n"
          "%s\n",
          rseq->head,
          rseq->seq);
}

struct derep_arg_t {
  int worker_num;
  int num_workers;
  int prefilter_len;
  tommy_array* seqs;
  unsigned long max_seq_len;
  /* tommy_array* unique_seq_indices; */
  unsigned long num_seqs;
};

struct derep_arg_t*
derep_arg_init(int worker_num,
               int num_workers,
               int prefilter_len,
               tommy_array* seqs,
               unsigned long max_seq_len,
               unsigned long num_seqs)
{
  struct derep_arg_t* derep_arg = malloc(sizeof(struct derep_arg_t));
  assert(derep_arg != NULL);

  derep_arg->worker_num = worker_num;
  derep_arg->num_workers = num_workers;
  derep_arg->prefilter_len = prefilter_len;
  derep_arg->seqs = seqs;
  derep_arg->max_seq_len = max_seq_len;
  derep_arg->num_seqs = num_seqs;

  return derep_arg;
}

void
derep_arg_destroy(struct derep_arg_t* derep_arg)
{
  free(derep_arg);
}

struct rseq_t*
derep_arg_get_seq(struct derep_arg_t* arg, unsigned long idx)
{
  return (struct rseq_t*)tommy_array_get(arg->seqs, idx);
}

struct fingerprint_t {
  uint64_t fingerprint;
  tommy_node node;
};

struct fingerprint_t*
fingerprint_init(uint64_t fingerprint)
{
  struct fingerprint_t* fp = malloc(sizeof(struct fingerprint_t));
  assert(fp != NULL);
  fp->fingerprint = fingerprint;

  return fp;
}

int
fingerprint_compare(const void* arg, const void* fp)
{
  return *(const uint64_t*)arg !=
    ((const struct fingerprint_t*)fp)->fingerprint;
}

void*
derep(void* derep_arg)
{

  /* TODO each thread will have this big array.... */
  int result = -2;
  unsigned long pattern_i = 0;
  unsigned long text_i = 0;

  struct derep_arg_t* arg = derep_arg;

  uint64_t* hash_vals = malloc(arg->max_seq_len * sizeof(uint64_t));
  assert(hash_vals != NULL);
  uint32_t num_hash_vals = 0;
  uint32_t hv_i = 0;
  uint64_t fingerprint = 0;

  struct fingerprint_t* fp;

  struct rseq_t* pattern;
  struct rseq_t* text;

  char* kmer = malloc((arg->prefilter_len + 1) * sizeof(char));
  assert(kmer != NULL);

  tommy_hashlin* text_hash = NULL;

  /* apple(WAIT); */
  /* tommy_array* unique_indices = malloc(sizeof(tommy_array)); */
  /* assert(unique_indices != NULL); */
  /* tommy_array_init(unique_indices); */


  for (text_i = arg->worker_num;
       text_i < arg->num_seqs - 1; /* go until next to last */
       text_i += arg->num_workers) {

    if (text_i % 10 == 0) {
      fprintf(stderr,
              "LOG -- Checking seq %lu\r",
              text_i);
    }

    text_hash = malloc(sizeof(tommy_hashlin));
    assert(text_hash != NULL);
    tommy_hashlin_init(text_hash);

    text = derep_arg_get_seq(arg, text_i);
    num_hash_vals = set_hash_vals(hash_vals, text->seq, text->len);

    for (hv_i = 0; hv_i < num_hash_vals; ++hv_i) {
      fp = tommy_hashlin_search(text_hash,
                                fingerprint_compare,
                                &fingerprint,
                                fingerprint);

      if (!fp) { /* new hash val */
        fp = fingerprint_init(hash_vals[hv_i]);
        tommy_hashlin_insert(text_hash,
                             &fp->node,
                             fp,
                             fp->fingerprint);
      }
    }

    for (pattern_i = text_i + 1; /* the next shortest pattern */
         pattern_i < arg->num_seqs;
         ++pattern_i) {

      fingerprint = *((uint64_t*)tommy_array_get(patterns, pattern_i));
      fp = tommy_hashlin_search(text_hash,
                                fingerprint_compare,
                                &fingerprint,
                                fingerprint);

      if (fp) { /* pattern_i is found in the text */
        /* passed the filter, check the whole seq against the text */
        pthread_mutex_lock(&mutex);
        ++passed_prefilter;
        pthread_mutex_unlock(&mutex);

        pattern = derep_arg_get_seq(arg, pattern_i);

        /* search the whole pattern */
        result = hash3_search((unsigned char*)pattern->seq,
                              pattern->len,
                              (unsigned char*)text->seq,
                              text->len);

        if (result > 0) { /* match! pattern is not unique */
          pthread_mutex_lock(&mutex);
          ++non_unique_seqs;
          non_unique_seq_index[pattern_i] = 1;
          pthread_mutex_unlock(&mutex);
        }
      }
    }

    tommy_hashlin_foreach(text_hash, free);
    tommy_hashlin_done(text_hash);
    free(text_hash);
  }


  free(arg);
  free(kmer);
  free(hash_vals);

  return NULL;
}

int main(int argc, char *argv[])
{
  if (argc < 3) {
    fprintf(stderr,
            "Version: %s --- Usage: %s <num threads> "
            "<contigs.fasta>\n",
            VERSION,
            argv[0]);

    return 1;
  }

  kseq_t* kseq;
  gzFile fp;

  fp = gzopen(argv[2], "r");
  if (!fp) {
    fprintf(stderr, "ERROR -- Could not open %s\n", argv[2]);
    return 2;
  }

  kseq = kseq_init(fp);
  assert(kseq);

  pthread_mutex_init(&mutex, NULL);

  int ret_code = 0;
  unsigned long i = 0;
  unsigned long idx = 0;
  long l = 0;

  uint64_t fingerprint;

  unsigned long max_seq_len;

  unsigned long num_seqs = 0;
  unsigned long num_workers = strtol(argv[1], NULL, 10);

  fprintf(stderr, "LOG -- num worker threads: %lu\n", num_workers);
  fprintf(stderr, "LOG -- infile: %s\n", argv[2]);

  tommy_array* unique_indices = malloc(num_workers * sizeof(tommy_array*));
  assert(unique_indices != NULL);

  struct rseq_t* rseq;
  struct derep_arg_t* derep_arg = NULL;

  pthread_t* threads = malloc(num_workers * sizeof(pthread_t));
  assert(threads != NULL);

  tommy_array* seqs = malloc(sizeof(tommy_array));
  assert(seqs != NULL);
  tommy_array_init(seqs);

  patterns = malloc(sizeof(tommy_array));
  assert(patterns != NULL);
  tommy_array_init(patterns);

  char* buf = malloc((KMER_LEN + 1) * sizeof(char));
  assert(buf != NULL);

  fprintf(stderr, "LOG -- reading seqs into memory\n");
  while ((l = kseq_read(kseq)) >= 0) {
    rseq = rseq_init(kseq);
    if (rseq->len > max_seq_len) {
      max_seq_len = rseq->len;
    }

    if (idx == 0) {
      rseq_print(stdout, rseq);
    }

    if ((++idx % 1000) == 0) {
      fprintf(stderr, "LOG -- Reading seq %lu\r", idx);
    }

    if (rseq->len < KMER_LEN) {
      rseq_print(stdout, rseq);
      rseq_destroy(rseq);
    } else {
      tommy_array_insert(seqs, rseq);

      for (i = 0; i < KMER_LEN; ++i) {
        buf[i] = rseq->seq[i];
      }
      buf[i] = '\0';

      fingerprint = rabin_fingerprint(buf, KMER_LEN);
      uint64_t* tmp = malloc(sizeof(uint64_t));
      assert(tmp != NULL);
      tmp[0] = fingerprint;
      tommy_array_insert(patterns, tmp);
    }
  }
  num_seqs = tommy_array_size(seqs);
  non_unique_seq_index = calloc(num_seqs, sizeof(int));
  assert(non_unique_seq_index != NULL);
  fprintf(stderr, "LOG -- Reading seq %lu\n", idx);


  for (i = 0; i < num_workers; ++i) {
    derep_arg = derep_arg_init(i, num_workers, KMER_LEN, seqs, max_seq_len, num_seqs);
    ret_code = pthread_create(&threads[i], NULL, derep, derep_arg);

    if (ret_code) {
      fprintf(stderr, "Error creating thread #%lu\n", i);

      return 3;
    }
  }

  for (i = 0; i < num_workers; ++i) {
    ret_code = pthread_join(threads[i], NULL);

    if (ret_code) {
      fprintf(stderr, "Error joining thread\n");

      return 4;
    }
  }
  /* TODO when short seqs are printed, this number wont match the
     number read. Could be confusing for user. */
  fprintf(stderr, "LOG -- Checking seq %lu\n", num_seqs);

  fprintf(stderr,
          "LOG -- %lu seqs passed the prefilter\n",
          passed_prefilter);

  fprintf(stderr,
          "LOG -- %lu seqs were not unique\n",
          non_unique_seqs);

  /* print out the non unique seqs */
  for (i = 0; i < num_seqs; ++i) {
    rseq = tommy_array_get(seqs, i);
    if (non_unique_seq_index[i] == 0) {
      rseq_print(stdout, rseq);
    }
    rseq_destroy(rseq);
  }

  tommy_array_done(seqs);
  free(seqs);
  free(threads);

  gzclose(fp);
  kseq_destroy(kseq);

  free(unique_indices);

  pthread_mutex_destroy(&mutex);

  free(buf);

  for (i = 0; i < tommy_array_size(patterns); ++i) {
    free(tommy_array_get(patterns, i));
  }
  tommy_array_done(patterns);
  free(patterns);

  free(non_unique_seq_index);

  return 0;

}
