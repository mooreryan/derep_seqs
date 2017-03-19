#ifndef SSEF_H
#define SSEF_H

#include <emmintrin.h>
#include <stdlib.h>
#include <assert.h>
#include "tommyarray.h"
#include <string.h>
typedef union{
  __m128i* data16;
  unsigned char* data;
} TEXT;

typedef struct list {
  struct list *next;
  int pos;
} LIST;

int
ssef_search(unsigned char* x, int Plen, unsigned char *y, int Tlen);

#endif
