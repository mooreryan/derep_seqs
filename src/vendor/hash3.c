#include "hash3.h"
#include <stdio.h>

int hash3_search(unsigned char *x, int m, unsigned char *y, int n) {
  int count, j, i, sh, sh1, mMinus1, mMinus2, shift[WSIZE];
  unsigned char h;
  if (m<3) return -1;
  count = 0;
  mMinus1 = m-1;
  mMinus2 = m-2;

  /* Preprocessing */
  for (i = 0; i < WSIZE; ++i)
    shift[i] = mMinus2;

  h = x[0];
  h = ((h<<1) + x[1]);
  h = ((h<<1) + x[2]);
  shift[h] = m-RANK3;
  for (i=RANK3; i < mMinus1; ++i) {
    h = x[i-2];
    h = ((h<<1) + x[i-1]);
    h = ((h<<1) + x[i]);
    shift[h] = mMinus1-i;
  }
  h = x[i-2];
  h = ((h<<1) + x[i-1]);
  h = ((h<<1) + x[i]);
  sh1 = shift[h];
  shift[h] = 0;
  if (sh1==0) sh1=1;


  /* Searching */
  i = mMinus1;
  /* TODO this line blows up */
  unsigned char* x_y_cat = calloc((n + m + 2), sizeof(unsigned char));
  assert(x_y_cat != NULL);
  int ai = 0;
  int bi = 0;
  /* cat the strings */
  for (ai = 0; ai < n; ++ai) {
    x_y_cat[ai] = y[ai];
  }
  for (bi = 0; bi < m; ++bi) {
    x_y_cat[ai+bi] = x[bi];
  }
  x_y_cat[ai+bi] = '\0';
  /* strcat((char*)x_y_cat, (const char*)y); */
  /* strcat((char*)x_y_cat, (const char*)x); */
  /* sprintf((char*)x_y_cat, "%s%s", y, x); */
  /* memcpy(y+n, x, m); */
  while (1) {
    if (count > 0) {
      free(x_y_cat);
      return count;
    }
    sh = 1;
    while (sh != 0) {
      /* TODO almost 50% of L1 cache misses occur at this line */
      h = x_y_cat[i-2];

      h = ((h<<1) + x_y_cat[i-1]);
      h = ((h<<1) + x_y_cat[i]);
      sh = shift[h];
      i+=sh;
    }
    if (i < n) {
      j=0;
      while(j<m && x[j]==x_y_cat[i-mMinus1+j]) j++;
      if (j>=m) {
        OUTPUT(i-mMinus1);
      }
      i+=sh1;
    }
    else {
      free(x_y_cat);
      return count;
    }
  }
}
