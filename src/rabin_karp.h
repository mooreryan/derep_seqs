#ifndef RABIN_KARP_H
#define RABIN_KARP_H

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* TODO the integer overflow calculation is wrong..17, 13 seems to work tho */
#define PRIME 7
#define KMER_LEN 18

uint64_t
rabin_fingerprint(char* text, uint32_t text_len);

uint32_t
set_hash_vals(uint64_t* hash_vals, char* text, uint32_t text_len);

#endif
