CC = gcc
CFLAGS = -Wall -Wno-unused-function -g -O2
SRC = src
VENDOR = $(SRC)/vendor
BIN = bin
OBJS = $(VENDOR)/tommyarray.o $(VENDOR)/ssef.o $(VENDOR)/hash3.o \
       $(SRC)/pointers.o $(SRC)/rseq.o

.PHONY: all
all: derep_seqs

derep_seqs: $(OBJS)
	$(CC) $(CFLAGS) -lz -o $(BIN)/$@ $^ src/$@.c

.PHONY: clean
clean:
	-rm -r $(BIN)/derep_seqs $(BIN)/derep_seqs.dSYM $(OBJS)
