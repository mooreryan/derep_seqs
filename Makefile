CC = gcc
CFLAGS = -Wall -Wno-char-subscripts -g -O2
LDFLAGS = -lm -lz -lpthread
SRC = src
VENDOR = $(SRC)/vendor
BIN = bin
OBJS = $(SRC)/rabin_karp.o \
       $(VENDOR)/tommyarray.o \
       $(VENDOR)/tommyhashlin.o \
       $(VENDOR)/tommyhash.o \
       $(VENDOR)/hash3.o


.PHONY: all
all: derep_seqs derep_seqs_non_rab

derep_seqs: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ src/$@.c $(LDFLAGS)

derep_seqs_non_rab: $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/$@ $^ src/$@.c $(LDFLAGS)

.PHONY: clean
clean:
	-rm -r $(BIN)/derep_seqs $(BIN)/derep_seqs.dSYM $(OBJS)
