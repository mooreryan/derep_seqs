#include <assert.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <zlib.h>

#include "vendor/hash3.h"
#include "vendor/kseq.h"
#include "vendor/tommyarray.h"

#define VERSION "0.6.0"
pthread_mutex_t mutex;
unsigned long passed_prefilter = 0;
unsigned long non_unique_seqs = 0;

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
  /* tommy_array* unique_seq_indices; */
  unsigned long num_seqs;
};

struct derep_arg_t*
derep_arg_init(int worker_num,
               int num_workers,
               int prefilter_len,
               tommy_array* seqs,
               unsigned long num_seqs)
{
  struct derep_arg_t* derep_arg = malloc(sizeof(struct derep_arg_t));
  assert(derep_arg != NULL);

  derep_arg->worker_num = worker_num;
  derep_arg->num_workers = num_workers;
  derep_arg->prefilter_len = prefilter_len;
  derep_arg->seqs = seqs;
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

void*
derep(void* derep_arg)
{

  int i = 0;
  int result = -2;
  unsigned long pattern_i = 0;
  unsigned long text_i = 0;

  struct derep_arg_t* arg = derep_arg;

  struct rseq_t* pattern;
  struct rseq_t* text;

  char* kmer = malloc((arg->prefilter_len + 1) * sizeof(char));
  assert(kmer != NULL);

  /* apple(WAIT); */
  /* tommy_array* unique_indices = malloc(sizeof(tommy_array)); */
  /* assert(unique_indices != NULL); */
  /* tommy_array_init(unique_indices); */


  for (pattern_i = 1 + arg->worker_num;
       pattern_i < arg->num_seqs;
       pattern_i += arg->num_workers) {

  /* for (pattern_i = 1; pattern_i < arg->num_seqs; ++pattern_i) { */
  /*   /\* should this thread do the work? *\/ */
  /*   if (pattern_i % arg->num_workers == arg->worker_num) { */
    if (pattern_i % 100 == 0) {
      fprintf(stderr,
              "LOG -- Checking seq %lu\r",
              pattern_i);
    }
    result = 0;
    pattern = derep_arg_get_seq(arg, pattern_i);

    for (text_i = 0; text_i < pattern_i; ++text_i) {
      text = derep_arg_get_seq(arg, text_i);

      if (text->len == pattern->len) {
        result = fast_compare(text->seq, pattern->seq, text->len);
        if (result == 0) { /* match */
          result = 1;
        } else {
          result = 0;
        }
      } else {

        for (i = 0; i < arg->prefilter_len; ++i) {
          kmer[i] = pattern->seq[i];
        }
        kmer[i] = '\0';

        /* search first derep_arg->prefilter_len chars of pattern */
        result = hash3_search((unsigned char*)kmer,
                              arg->prefilter_len,
                              (unsigned char*)text->seq,
                              text->len);

        if (result > 0) { /* prefilter match */
          pthread_mutex_lock(&mutex);
          ++passed_prefilter;
          pthread_mutex_unlock(&mutex);

          /* search the whole pattern */
          result = hash3_search((unsigned char*)pattern->seq,
                                pattern->len,
                                (unsigned char*)text->seq,
                                text->len);
          if (result > 0) { /* match! pattern is not unique */
            pthread_mutex_lock(&mutex);
            ++non_unique_seqs;
            pthread_mutex_unlock(&mutex);

            break;
          }
        }
      }
    }

    /* TODO is it faster to save indices and print at the end? */
    /* TODO does this need locking? */
    if (result == 0) { /* seq is unique */
      rseq_print(stdout, pattern);
    }
  }

  free(arg);
  free(kmer);

  return NULL;
}

int main(int argc, char *argv[])
{
  if (argc < 4) {
    fprintf(stderr,
            "Version: %s --- Usage: %s <num threads> <prefilter len> "
            "<contigs.fasta>\n",
            VERSION,
            argv[0]);

    return 1;
  }

  kseq_t* kseq;
  gzFile fp;

  fp = gzopen(argv[3], "r");
  if (!fp) {
    fprintf(stderr, "ERROR -- Could not open %s\n", argv[3]);
    return 2;
  }

  kseq = kseq_init(fp);
  assert(kseq);

  pthread_mutex_init(&mutex, NULL);

  int ret_code = 0;
  unsigned long i = 0;
  unsigned long idx = 0;
  long l = 0;

  unsigned long num_seqs = 0;
  unsigned long num_workers = strtol(argv[1], NULL, 10);
  int prefilter_len = strtol(argv[2], NULL, 10);

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

  fprintf(stderr, "LOG -- reading seqs into memory\n");
  while ((l = kseq_read(kseq)) >= 0) {
    rseq = rseq_init(kseq);

    if (idx == 0) {
      rseq_print(stdout, rseq);
    }

    if ((++idx % 1000) == 0) {
      fprintf(stderr, "LOG -- Reading seq %lu\r", idx);
    }

    if (rseq->len < prefilter_len) {
      rseq_print(stdout, rseq);
      rseq_destroy(rseq);
    } else {
      tommy_array_insert(seqs, rseq);
    }
  }
  num_seqs = tommy_array_size(seqs);
  fprintf(stderr, "LOG -- Reading seq %lu\n", idx);


  for (i = 0; i < num_workers; ++i) {
    derep_arg = derep_arg_init(i, num_workers, prefilter_len, seqs, num_seqs);
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

  for (i = 0; i < num_seqs; ++i) {
    rseq_destroy(tommy_array_get(seqs, i));
  }
  tommy_array_done(seqs);
  free(seqs);
  free(threads);

  gzclose(fp);
  kseq_destroy(kseq);

  free(unique_indices);

  pthread_mutex_destroy(&mutex);

  return 0;

}
