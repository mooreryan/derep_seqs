# Derep Seqs

Dereplicate looooooong sequences fast!

If you want to get rid of duplicate long sequences (i.e. those sequences that are exact substrings of some other sequence), `derep_seqs` is the tool for you!

## Install

Download the source code (either with `git clone` or by downloading a release), `cd` into the source directory, and then use `make` to build it.

```
git clone https://github.com/mooreryan/derep_seqs.git
cd derep_seqs
make
```

This will install `derep_seqs` to the `bin` directory in the source directory. You can now move `derep_seqs` and `sort_fasta` to somewhere on your path if you'd like.

## Background

Substring matching...there are a million algorithms for this problem. This [paper](https://arxiv.org/pdf/1012.2547v1.pdf) does a cool job of testing 80 or so of them with different pattern sizes and alphabet lengths. For nucleotide and amino acid alphabets, the authors found the SSEF algorithm to be the best for patterns at least 256 characters long. So, `derep_seqs` uses this nifty algorithm!

## Example

The fasta file must be sorted by increasing sequence length. The program `sort_fasta` (included in the `bin` directory) will do this for you.

```
$ bin/derep_seqs <(bin/sort_fasta contigs.fasta) > contigs.derep.fa
```
That's it!
