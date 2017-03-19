/*
Copyright (c) 2012, Simone Faro and Thierry Lecroq.
http://www.dmi.unict.it/~faro/
http://www-igm.univ-mlv.fr/~lecroq/

This tool is free software; you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation; either version 3.0 of the License, or
(at your option) any later version.  (visit
http://www.gnu.org/licenses/gpl.html)

This tool is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA.
*/

#include "../pointers.h"
#include "ssef.h"

/* from https://www.dmi.unict.it/~faro/smart/algorithms.php?algorithm=SSEF&code=ssef */
/* TODO is that len with '\0' or without? */
int
ssef_search(unsigned char* x, int Plen, unsigned char *y, int Tlen)
{
  struct pointers_t* pointers = pointers_init();

  if (Plen < 32) {
    pointers_destroy(pointers);
    return -1;
  }

  LIST *flist[65536];
  LIST *t;
  memset(flist,0,sizeof(LIST*)*65536);
  /* __m128i tmp128; */
  TEXT T;
  T.data16 = (__m128i*) y;
  T.data = (unsigned char *) y;

  /* unsigned int K=7; */
  unsigned int count=0;
  unsigned int i,last,j;
  __m128i *ptr16;
  __m128i *lastchunk = &T.data16[Tlen/16];
  unsigned int filter;
  unsigned char* f = malloc(Plen);
  assert(f != NULL);

  last = (Plen/16) - 1;
  for (i=0;i<Plen;i++){
    f[i] = (x[i]&0x80)>>7;
  }
  count = 15;

  for (i=0;i<last*16;i++) {
    j = last*16-i;
    filter =
      f[j]         + f[j+1]*2     + f[j+2]*4      + f[j+3]*8      +
      f[j+4]*16    + f[j+5]*32    + f[j+6]*64     + f[j+7]*128    +
      f[j+8]*256   + f[j+9]*512   + f[j+10]*1024  + f[j+11]*2048  +
      f[j+12]*4096 + f[j+13]*8192 + f[j+14]*16384 + f[j+15]*32768 ;
    if (flist[filter]==0){
      /* TODO this might leak memory */
      flist[filter]       = (LIST*)malloc(sizeof(LIST));
      assert(flist[filter] != NULL);
      pointers_add_pointer(pointers, flist[filter]);
      flist[filter]->next = NULL;
      flist[filter]->pos  = i;
    } else {
      t = flist[filter];
      while(t->next!=NULL) t = t->next;
      /* TODO this leaks memory */
      t->next = (LIST*)malloc(sizeof(LIST));
      assert(t->next != NULL);
      pointers_add_pointer(pointers, t->next);
      t       = t->next;
      t->next = NULL;
      t->pos  = i;
    }
  }

  count = 0;
  ptr16 = &T.data16[last];

  while(ptr16 < lastchunk) {
    filter = _mm_movemask_epi8(*ptr16);

    if (flist[filter]) {
      i = ((ptr16 - &T.data16[0])-last)*16;
      t = flist[filter];
      while(t) {
        if (memcmp(x, &T.data[i+t->pos], Plen) == 0) {
          count++;

          /* TODO this should save some cycles since we only need one
             match for the pattern to be non-unique. Double check
             it. */

          pointers_destroy(pointers);
          free(f);

          return count;
        }
        t = t->next;
      }
    }
    ptr16 += last;
  }

  pointers_destroy(pointers);
  free(f);

  return count;
}
