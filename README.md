# kahn-ngs

A small package/script used to count the indexes/sequences observed in an
illumina DNA sequencing experiment.

The Kahn lab is using some very carefully crafted DNA fragments in order to
infer the relative energetics of inter/intra molecular ligation at single
nucleotide resolution.  (At least, this is my outsider's, ignorant
interpretation).

This script does the following:

1.  Reads a user-supplied table of sequencing indexes and parses their names in
    order to get some information about the various source DNA fragments.
2.  Iterates over every sequence in a set of fastq files and searches for likely
    hits of each index.
3.  Uses some logic to attempt to properly count the found indexes and infer
    from them the identity of the template(s) which generated the sequence.
4.  Prints the found indexes and their positions into the comment field of the
    fastq file(s).
5.  Creates and prints a set of summaries of the observed counts.

Along the way, it uses String::Approx to handle the matching, thus providing
some flexibility for mismatched indexes.

# Usage

The available options are:

* bicyc_csv:  Filename into which to write the counts of bimolecular results.
  (bimolecular_cyclized_lengths.csv by default)
* debug:  Print extra information at runtime (on by default)
* distance:  Levenstein distance for matching (0 by default)
* fourhitlin_csv:  Filename into which to write the counts of unimolecular
  results.  (unimolecular_linear_lengths.csv)
* indices:  Filename containing the set of indexes. (index.txt)
* input:  Input fastq file (test.fastq.gz)
* outdir:  Directory into which to write the outputs. (output)
* outfastq:  Prefix for the output fastq file containing the new comments. (out)
* spacer:  A constant added to each successful iteration to calculate the final
    sequence length.
* summary:  File into which to write the summary of the results. (summary.txt)
* unicyc_csv:  File to write the lengths of the unimolecular
  results. (unimolecular_cyclized_lengths.csv)
