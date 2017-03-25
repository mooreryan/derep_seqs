#include "rabin_karp.h"

static uint64_t
power(int index)
{
  return llround(powl(PRIME, KMER_LEN - index - 1));
}

uint64_t
rabin_fingerprint(char* text, uint32_t text_len)
{
  uint64_t hash = 0;
  uint32_t i = 0;
  int chr = 0;

  int map_base['T' + 1];
  map_base['A'] = 0;
  map_base['C'] = 1;
  map_base['T'] = 2;
  map_base['G'] = 3;

  /* TODO check text_len > KMER_LEN */
  for (i = 0; i < KMER_LEN; ++i) {
    chr = map_base[text[i]];
    hash += chr * power(i);
  }

  return hash;
}

/* convert text A C T G to 0 1 2 3 first */
/* set the hash vals for each kmer in text */
uint32_t
set_hash_vals(uint64_t* hash_vals, char* text, uint32_t text_len)
{
  int map_base['T' + 1];
  map_base['A'] = 0;
  map_base['C'] = 1;
  map_base['T'] = 2;
  map_base['G'] = 3;

  uint32_t num_kmers = text_len - KMER_LEN + 1;
  uint32_t i = 0;

  uint64_t hash = 0;
  uint64_t tmp = 0;

  int chr = 0;
  int last_first_char = 0;
  int first_char = 0;
  int last_char = 0;

  /* calculate first hash val */
  for (i = 0; i < KMER_LEN; ++i) {
    chr = map_base[text[i]];
    tmp = chr * power(i);
    hash += tmp;

    /* if (i == 0) { /\* store the first char hash *\/ */
    /*   prev_first_char_hash = tmp; */
    /* } */
  }
  hash_vals[0] = hash;

  /* for (i = 1; i < num_kmers; ++i) { */
  /*   last_first_chr = map_base[text[i-1]]; */
  /*   first_chr = map_base[text[i]]; */
  /*   last_chr  = map_base[text[i + KMER_LEN - 1]]; */

  /*   /\* base * (prev full hash - prev first char hash) + new char hash *\/ */
  /*   hash = PRIME * (hash_vals[i - 1] - (last_first_chr * power(0))) + last_chr; */

  /*   hash_vals[i] = hash; */
  /*   prev_first_char_hash = first_chr * power(0); */
  /* } */

  /* for (i = 1; i < num_kmers; ++i) { */
  /*   hash = 0; */
  /*   for (j = 0; j < KMER_LEN; ++j) { */
  /*     chr = map_base[text[i + j]]; */
  /*     hash += chr * power(j); */
  /*   } */
  /*   hash_vals[i] = hash; */
  /* } */

  for (i = 1; i < num_kmers; ++i) {
    hash = 0;
    last_first_char = map_base[text[i-1]];
    first_char = map_base[text[i]];
    last_char = map_base[text[i + KMER_LEN - 1]];

    hash = (PRIME * (hash_vals[i-1] - (last_first_char * power(0)))) + last_char;
    hash_vals[i] = hash;
  }

  return num_kmers;
}

/* int main(int argc, char *argv[]) */
/* { */
/*   uint32_t i = 0; */
/*   uint32_t num_kmers = 0; */
/*   uint64_t* hash_vals = calloc(25, sizeof(uint64_t)); */
/*   assert(hash_vals != NULL); */

/*   int text_len = 4; */

/*   int* text = malloc(text_len * sizeof(int)); */
/*   assert(text != NULL); */

/*   text[0] = 0; */
/*   text[1] = 1; */
/*   text[2] = 2; */
/*   text[3] = 3; */

/*   num_kmers = set_hash_vals(hash_vals, text, text_len); */

/*   for (i = 0; i < num_kmers; ++i) { */
/*     printf("hash_vals[%u] == %llu\n", i, hash_vals[i]); */
/*   } */

/*   free(text); */
/*   free(hash_vals); */
/*   return 0; */
/* } */
