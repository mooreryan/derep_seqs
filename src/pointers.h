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

#ifndef POINTERS_H
#define POINTERS_H

#include <stdlib.h>
#include <assert.h>
#include "vendor/tommyarray.h"

struct pointers_t {
  tommy_array* pointers;
  unsigned long num_pointers;
};

struct pointers_t*
pointers_init();

void
pointers_destroy(struct pointers_t* pointers);

void
pointers_add_pointer(struct pointers_t* pointers, void* ptr);

#endif
