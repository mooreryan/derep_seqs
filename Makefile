CC = gcc
CFLAGS = -Wall -Wno-unused-function -g -O2
LDFLAGS = -lz
SRC = src
VENDOR = $(SRC)/vendor
BIN = bin
OBJS = $(VENDOR)/tommyarray.o $(VENDOR)/tommyhashlin.o \
       $(VENDOR)/tommyhash.o \
       $(VENDOR)/ssef.o $(VENDOR)/hash3.o \
       $(SRC)/pointers.o $(SRC)/rseq.o

.PHONY: all
all: derep_seqs

derep_seqs: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ src/$@.c $(LDFLAGS)

.PHONY: clean
clean:
	-rm -r $(BIN)/derep_seqs $(BIN)/derep_seqs.dSYM $(OBJS)
