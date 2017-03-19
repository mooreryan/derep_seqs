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

#include "pointers.h"

struct pointers_t*
pointers_init()
{
  struct pointers_t* pointers
    = malloc(sizeof(struct pointers_t));
  assert(pointers != NULL);

  pointers->pointers = malloc(sizeof(tommy_array));
  assert(pointers->pointers);
  tommy_array_init(pointers->pointers);
  pointers->num_pointers = 0;

  return pointers;
}

void
pointers_destroy(struct pointers_t* pointers)
{
  unsigned long i = 0;
  for (i = 0; i < pointers->num_pointers; ++i) {
    free(tommy_array_get(pointers->pointers, i));
  }
  tommy_array_done(pointers->pointers);
  free(pointers->pointers);

  free(pointers);
}

void
pointers_add_pointer(struct pointers_t* pointers, void* ptr)
{
  tommy_array_insert(pointers->pointers, ptr);
  ++pointers->num_pointers;
}
