# Derep Seqs

Dereplicate looooooong sequences!

If you want to get rid of duplicate long sequences (i.e. contigs that are exact substrings of some other contigs), `derep_seqs` is the tool for you!

## Install

Download the source code (either with `git clone` or by downloading a release), `cd` into the source directory, and then use `make` to build it.

```
git clone https://github.com/mooreryan/derep_seqs.git
cd derep_seqs
make
```

This will install `derep_seqs` to the `bin` directory in the source directory. You can now move `derep_seqs` and `sort_fasta` to somewhere on your path if you'd like.

## Usage

```
derep_seqs <num worker threads> <prefilter length> <seqs.fasta> > seqs.derep.fa
```

### Prefilter length

Prefilter length is the length of the first kmer of the pattern sequence to check. The longer it is, the fewer sequences pass the prefilter when they really are not substrings (making the program faster). On the other hand, the if prefilter length is *too* long, then this benefit is outweighed by the longer substring to check.

Around 60 seems to be a sweet spot.

## Example

The fasta file must be sorted by increasing sequence length. The program `sort_fasta` (included in the `bin` directory) will do this for you.

```
$ bin/derep_seqs 10 60 <(bin/sort_fasta contigs.fasta) > contigs.derep.fa
```
That's it!

## Error codes

- 0: Success
- 1: Argument error
- 2: Couldn't open a file
- 3: Error creating thread
- 4: Error joining thread

## Versions

- v0.1.0: First release
- v0.2.0: Sort on decreasing seq length. Use greedy algorithm. Prefilter. Use hash3 instead of SSEF.
- v0.3.0: Use hashing for prefiltering.
- v0.4.0: Don't store hash vals...uses way less memory :) but it's slow again :(
- v0.5.0: Use pthreads for multithreading!
- v0.6.0: Make prefilter length a tunable option
