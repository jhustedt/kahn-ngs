#!/usr/bin/env perl
use Modern::Perl;
use autodie qw":all";
use warnings qw"all";
use diagnostics;
use strict;

use Bio::Seq;
use Bio::SeqIO;
use Cwd qw(abs_path getcwd);
use File::Basename;
use File::Find;
use FileHandle;
use File::Path qw"make_path";
use File::Spec qw"rel2abs";
use Getopt::Long qw"GetOptions";
use PerlIO;
use PerlIO::gzip;
use String::Approx qw"amatch aindex";

$SIG{INT} = \&End_Handler;
$SIG{TERM} = \&End_Handler;

=head1 NAME

  sizecount.pl - A writeall script to sort sequences by index and infer the
    sizes of the template DNA.

=head1 SYNOPSIS

  This script only has a few useful options:
 --index  : The file defining the search strings and output files.
 --input  : Either a filename or directory containing the (relatively) raw
    fastq data.
 --outdir : Output directory, composed of two csv files, an output fastq file,
    and a summary file.
 --outfastq : Output fastq file, by default this is placed into the output
    directory. File extension not necessary.
 --substitution : String::Approx (Levenshtein) substitution distance, default is 0 (perfect match)
 --insertion: Levenshtein insertion distance
 --deletion: Levenshtein deletion distance
 --spacer : specify spacer length (default 72)
 --debug  : Print a bunch of debugging output

=head1 DESCRIPTION

  This was originally a demultiplexer for TNSeq data, which has a very specific and
  illumina-incompatible (bcl2fastq) format.  It was slightly reformatted to match
  the following query:

  1. The primary goal is to observe how many of every size DNA template was incorporated into
     a sequencing library.  The implicit inference: when (DNA-size % ~10), the energetic
     requirements for a DNA to loop back on itself decrease, increasing the likelihood that
     the molecule can cyclize and therefore end up in the sequencing library.  Ergo, more
     counts of a particular size is inversely related to the energy requirement.
  2. The input molecules were designed so that each size is distinguishable by its sequence.
  3. Four indices are used to distinguish molecule size (ranging 119-219 bp)
  4. The important caveat: inter-molecular ligation is in a fight with intra-molecular ligation.
     a. Ergo, we need to be aware of multiple indices for two reasons; first to tell the size
        range, second because we need to be able to observe when size-a linearly ligates
        to size-b and when size-a bimolecularly cyclizes with size-b.

  Thus the original implementation is useful but not ideal and has been changed to do the
  following:

  1. Rather than just dump the sequences into separate files, use the comment field
     to print the observed matches.
  2. When it does decide to write output files, it should have some logic to know the range(s)
     as well as the specific size(s) observed. The index file has been rewriten to include all
     relevant information within a single index file (for each diferent index type).
  3. Output a text file describing the number of observations of each-sized molecule that is
     cyclized or not cyclized, when observed with the appropriate number of hits for a single
     cycliztion event.
  4. Output individual CSV files for each cyclization state.

 Updates for 2018-11-13
 1) Verify that a sequence with three hits has the appropriate matching index before and after
    the cyclization site (this will tell us that it was likely unimolecular).
    a. added index sequences to allow for (1) via "stepsynth" and "stepcyc" implementation
    b. add logic to check by looking at index numbers for stepsynth and stepcyc & verifying match
 2) If a molecule is not unimolecular (as above) it is sorted into a bimolecular instead
    a. done for 4 hit, needs to be done for 6 hit still
    a. Does it matter that size-a ligated to size-b, or just that a molecule of each size
       participated in a bimolecular reaction?

 Updates for 2018-11-28
 1) Implemented amatch instead of direct string match for indices.
 2) the index number associated with stepsynth & stepcyc is now capable of being called upon to reference
    itself later for a logic check verifying they are the same index value.
    a. logic check separating unimolecular (matching) and bimolecular (not matching) stepsynth & stepcyc

 Updates for 2018-11-29
 1) fixed incorrect output sizes
 2) altered String::Approx to function as aindex() if params set and as index() if no params set
    a) this somewhat functions for aindex() but does NOT function for index()
        -aindex() no longer allows for unwanted insertions or deletions, but is not finding every match
         that it should

 Updates for 2018-12-02
 1) aindex still is not finding what I would expect, but index is better
 2) I have added counters for six hit bimolecular cyclization & two hit bimolecular cyclization.
    a) I don't do anything with these other than count their existence.
 3) I have added counters for excess hits and how many times they end up found.

 Updates for 2018-12-04
 1) added logic for six hit bimolecular and two hit bimolecular
    a) added logs for six hit and two hit bimolecular
 2) separated logs for bimolecular products - 2 hit, 4 hit, and 6 hit are reported in separate csv files

 Updates for 2019-03-12
 1) added matching for 3 hit unreacted (or dead end multimer).
 2) added matching for 1 hit unreacted (or dead end multimer).
 3) fixed 6-hit directionality
 4) added five hit bimolecular matching
 5) index file numbers changed for size class variant indices to include a leading zero in 47 and 77
    this is to make nomenclature uniform, all hits should report out as 3 digits size class, 2 digits
    variable index, 2 digits helical index. e.g. "0473010" - this allows for the index file input to take numbers directly and not need to reformat them here.
 6) full size description csv files will be generated with hits asigned to each index as approrpiate
 
 Needed to be fixed:
 1) fix aindex to find all direct matches and those with a single substitution

=cut

#options from Getopt::Long; defaults
my %options = (
    debug => 0,
    indices => 'index.txt',
    input => 'test.fastq.gz',
    outdir => 'output',
    outfastq => 'out',
    summary => 'summary.txt',
    unicyc_csv => 'unimolecular_ligated_lengths.csv',
    unicycfull_csv => 'unimolecular_ligated_lengths_full.csv',
    fourhitlin_csv => 'linear_fourhit_lengths.csv',
    fourhitlinfull_csv => 'linear_fourhit_lengths_full.csv',
    threehitlin_csv => 'linear_threehit_lengths.csv',
    threehitlinfull_csv => 'linear_threehit_lengths_full.csv',
    onehitlin_csv => 'linear_onehit_lengths.csv',
    bicyc4_csv => 'bimolecular4_ligated_lengths.csv',
    bicyc4full_csv => 'bimolecular4_ligated_lengths_full.csv',
    bicyc6_csv => 'bimolecular6_ligated_lengths.csv',
    bicyc6full_csv => 'bimolecular6_ligated_lengths_full.csv',
    bicyc5_csv => 'bimolecular5_ligated_lengths.csv',
    bicyc5full_csv => 'bimolecular5_ligated_lengths_full.csv',
    bicyc5frag_csv => 'bimolecular5_ligated_lengths_frag.csv',
    bicyc2_csv => 'bimolecular2_ligated_lengths.csv',
    bicyc2full_csv => 'bimolecular2_ligated_lengths_full.csv',
    spacer => 72,
    substitution => 0,
    insertion => 0,
    deletion => 0,
);
## This is a hash counting off how many times _every_ index is observed.
my $observed_reads = 0;
my %observations = (
    helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when 1-10 indices are observed.
my %singles = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %doubles = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %triples = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %quads = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %fives = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %sixes = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## in theory hits below here should not exist
my %sevens = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %eights = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %nines = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %tens = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
my %elevenup = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## Here I set up tags for counting all when all four are found unimolecular, when all four are found bimolecular, when all six are found biomolecular, and appropriate final lengths outputs for each
## currently this works by looking at the total summed length and not the individual makeup of the molecule.
my $datestring = localtime();
my $found_all_four = 0;
my $found_four_uni = 0;
my $found_four_lin = 0;
my $found_three_lin = 0;
my $found_one_lin = 0;
my $found_three_unknown = 0;
my $found_one_unknown = 0;
my $found_four_unknown = 0;
my %unicyclized4_final_lengths = ();
my %unicyclized4_full_final_lengths = ();
my %linear4_final_lengths = ();
my %linear4_full_final_lengths = ();
my %linear3_final_lengths = ();
my %linear3_full_final_lengths = ();
my %linear1_final_lengths = ();
my $found_four_bi = 0;
my %bicyclized4_final_lengths = ();
my %bicyclized4_full_final_lengths = ();
my $found_six_bi = 0;
my $found_six_unknown = 0;
my %bicyclized6_final_lengths = ();
my %bicyclized6_full_final_lengths = ();
my $found_five_bi = 0;
my $found_five_unknown = 0;
my %bicyclized5_final_lengths = ();
my %bicyclized5_full_final_lengths = ();
my %bicyclized5_frag_final_lengths = ();
my $found_two_bi = 0;
my $found_two_unknown = 0;
my %bicyclized2_final_lengths = ();
my %bicyclized2_full_final_lengths = ();
my $opt_result = GetOptions(
    "debug:i" => \$options{debug},
    "spacer:i" => \$options{spacer},
    "indices:s" => \$options{indices},
    "input:s" => \$options{input},
    "outdir:s" => \$options{outdir},
    "summary:s" => \$options{summary},
    "unicyc_csv:s" => \$options{unicyc_csv},
    "unicycfull_csv:s" => \$options{unicycfull_csv},
    "fourhitlin_csv:s" => \$options{fourhitlin_csv},
    "fourhitlinfull_csv:s" => \$options{fourhitlinfull_csv},
    "threehitlin_csv:s" => \$options{threehitlin_csv},
    "threehitlinfull_csv:s" => \$options{threehitlinfull_csv},
    "onehitlin_csv:s" => \$options{onehitlin_csv},
    "bicyc4_csv:s" => \$options{bicyc4_csv},
    "bicyc4full_csv:s" => \$options{bicyc4full_csv},
    "bicyc6_csv:s" => \$options{bicyc6_csv},
    "bicyc6full_csv:s" => \$options{bicyc6full_csv},
    "bicyc5_csv:s" => \$options{bicyc5_csv},
    "bicyc5full_csv:s" => \$options{bicyc5full_csv},
    "bicyc5frag_csv:s" => \$options{bicyc5frag_csv},
    "bicyc2_csv:s" => \$options{bicyc2_csv},
    "bicyc2full_csv:s" => \$options{bicyc2full_csv},
    "outfastq:s" => \$options{outfastq},
    "substitution:i" => \$options{substitution},
    "insertion:i" => \$options{insertion},
    "deletion:i" => \$options{deletion},
);
my $log = new FileHandle(">$options{outdir}/$options{summary}");
my $unicyc_csv = new FileHandle(">$options{outdir}/$options{unicyc_csv}");
my $unicycfull_csv = new FileHandle(">$options{outdir}/$options{unicycfull_csv}");
my $fourhitlin_csv = new FileHandle(">$options{outdir}/$options{fourhitlin_csv}");
my $fourhitlinfull_csv = new FileHandle(">$options{outdir}/$options{fourhitlinfull_csv}");
my $bicyc4_csv = new FileHandle(">$options{outdir}/$options{bicyc4_csv}");
my $bicyc4full_csv = new FileHandle(">$options{outdir}/$options{bicyc4full_csv}");
my $bicyc6_csv = new FileHandle(">$options{outdir}/$options{bicyc6_csv}");
my $bicyc6full_csv = new FileHandle(">$options{outdir}/$options{bicyc6full_csv}");
my $bicyc2_csv = new FileHandle(">$options{outdir}/$options{bicyc2_csv}");
my $bicyc2full_csv = new FileHandle(">$options{outdir}/$options{bicyc2full_csv}");
my $bicyc5_csv = new FileHandle(">$options{outdir}/$options{bicyc5_csv}");
my $bicyc5full_csv = new FileHandle(">$options{outdir}/$options{bicyc5full_csv}");
my $bicyc5frag_csv = new FileHandle(">$options{outdir}/$options{bicyc5frag_csv}");
my $threehitlin_csv = new FileHandle(">$options{outdir}/$options{threehitlin_csv}");
my $threehitlinfull_csv = new FileHandle(">$options{outdir}/$options{threehitlinfull_csv}");
my $onehitlin_csv = new FileHandle(">$options{outdir}/$options{onehitlin_csv}");
## This checks to see that the options for index, input, and output directory exist and provides an error if they do not.
if (!-r $options{input}) {
    die("The input file: $options{input} does not exist.");
}
if (!-r $options{indices}) {
    die("The index file: $options{indices} does not exist.");
}
if (!-d $options{outdir}) {
    die("The output directory $options{outdir} does not exist.");
}
my $index_hash = Read_indices();
## look at input file, if file, read, if dir, look in dir and read internal
## files
## Create an empty output fastq file into which we will copy the extant data and new fun comments.
my $abs_input = File::Spec->rel2abs($options{input});
my $fastq_output = qq"$options{outdir}/$options{outfastq}.fastq.gz";
my $abs_output = File::Spec->rel2abs($fastq_output);
my $out;
if ($abs_input eq $abs_output) {
    die("The input file and output file are the same.");
} else {
    ## generation of fastq file from another of Trey's scripts:
    ## my $with_ta = FileHandle->new("| gzip -f -9 > ${input_base}_ta.fastq.gz");
    $out = FileHandle->new("| gzip -f -9 > $fastq_output");
}
my $reads;
if (-f $options{input}) {
    $reads = Sort_File_Approx(
        input => $options{input},
        outdir => $options{outdir},
        index_hash => $index_hash,
    );
} elsif (-d $options{input}) {
    Sort_Dir(
        indir => $options{input},
        outdir => $options{outdir},
        index_hash => $index_hash,
    );
} else {
    die("I need a file/directory containing some sequence.") unless($options{input});
}
## general output, tells number of reads that went to each index
End_Handler();

=item Sort_File_Approx

  This function should look through every sequence for a reasonable match to
  the available indices.  If it gets some hits, record them.  If they are ambiguous,
  this should probably return where the possibilities lie.  If no match, return that.

=cut

sub Sort_File_Approx {
    my %args = @_;
    my $data = $args{index_hash};
    return(undef) unless ($args{input} =~ /\.fastq/);
    if (!-d $args{outdir}) {
        make_path($args{outdir});
    }
    my $inputted = FileHandle->new("zless $args{input} 2>/dev/null |");
    my $in = Bio::SeqIO->new(-fh => $inputted, -format => 'Fastq');
    my $count = 0;
    ## Here we define the order of lines within the fastq file and the contents of the comments.
  READS: while (my $in_seq = $in->next_dataset()) {
        $count++;
        $observed_reads++;
        my $id = $in_seq->{'-descriptor'};
        my $sequence = $in_seq->{'-seq'};
        my $qual = $in_seq->{'-raw_quality'};
        my $comment = $in_seq->{'-comment'};
        my $seqlen = length($sequence);
        if (!defined($comment)) {
            $comment = "";
        }
        my $found = 0;
        my @index_list = @{$data->{possibilities}};
        ## Reminder of what the data structure looks like: template => {
        ##         AAAATTTTCCC => { verb => 'add', direction => 'fwd',
        ##         number => 10, found => 0, total_observed => 0,
        ##         ambiguous_observed => 0, unique_observed => 0 }, },
        my $found_id = '';
        ## substring looking from 0-9, switched to regular expression match to
        ## scan full seq line my $match_substring = substr($sequence, 0, 9);
        ## I am creating a hash of observations for each sequence, and one for
        ## all sequences.
        my %observe = (
            helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
            helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,
            unknown => 0);
        my %positions = (
            helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
            helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,
            unknown => 0);
        my %numbers = (
            helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
            helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,
            unknown => 0);
        my $observed_indices = 0;
        foreach my $index (@index_list) {
            my $info = $data->{$index};
            ## direct sequence match, without error
            ## set up array for parameters for string::approx
            my $params = [ "I$options{insertion}","D$options{deletion}","S$options{substitution}" ];
            ##set up array for starts
            my @starts = ();
            ##looks at parameters & if not defined, continue with direct match, if params exist, do aindex
            if ( ! defined $options{insertion} && ! defined $options{deletion} &&
                     ! defined $options{substitution} ) {
                $sequence =~ m/$index/;
                @starts = @-;
            } else {
                @starts = aindexes($index, $params, ($sequence));
            }
            if (@starts) {
                $found++;
                $found_id = $index;
                for my $st (@starts) {
                    if ($options{debug} == 1) {
                        print "TESTME: @starts\n";
                    }
                    $comment .= "$st:$info->{name}:$info->{number}:$info->{direction} ";
                    ## where "st" is position within read, info name and number identify the index, and direction identifies fwd or rev
                    if ($info->{name} eq 'helical' && $info->{direction} eq 'fwd') {
                        $observations{helical_fwd}++;
                        $observe{helical_fwd}++;
                        $positions{helical_fwd} = $st;
                        $observed_indices++;
                        $numbers{helical_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'stepcyc' && $info->{direction} eq 'fwd') {
                        $observations{stepcyc_fwd}++;
                        $observe{stepcyc_fwd}++;
                        $positions{stepcyc_fwd} = $st;
                        $observed_indices++;
                        $numbers{stepcyc_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'stepsynth' && $info->{direction} eq 'fwd') {
                        $observations{stepsynth_fwd}++;
                        $observe{stepsynth_fwd}++;
                        $positions{stepsynth_fwd} = $st;
                        $observed_indices++;
                        $numbers{stepsynth_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'variable' && $info->{direction} eq 'fwd') {
                        $observations{variable_fwd}++;
                        $observe{variable_fwd}++;
                        $positions{variable_fwd} = $st;
                        $observed_indices++;
                        $numbers{variable_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'helical' && $info->{direction} eq 'rev') {
                        $observations{helical_rev}++;
                        $observe{helical_rev}++;
                        $positions{helical_rev} = $st;
                        $observed_indices++;
                        $numbers{helical_rev} = $info->{number};
                    } elsif ($info->{name} eq 'stepcyc' && $info->{direction} eq 'rev') {
                        $observations{stepcyc_rev}++;
                        $observe{stepcyc_rev}++;
                        $positions{stepcyc_rev} = $st;
                        $observed_indices++;
                        $numbers{stepcyc_rev} = $info->{number};
                    } elsif ($info->{name} eq 'stepsynth' && $info->{direction} eq 'rev') {
                        $observations{stepsynth_rev}++;
                        $observe{stepsynth_rev}++;
                        $positions{stepsynth_rev} = $st;
                        $observed_indices++;
                        $numbers{stepsynth_rev} = $info->{number};
                    } elsif ($info->{name} eq 'variable' && $info->{direction} eq 'rev') {
                        $observations{variable_rev}++;
                        $observe{variable_rev}++;
                        $positions{variable_rev} = $st;
                        $observed_indices++;
                        $numbers{variable_rev} = $info->{number};
                    }
                }
            }
        }         ## End each index
        ## to debug - print comment to standard out
        if ($options{debug} == 1) {
            print STDOUT "$comment ";
        }
        ## set counters to zero for each type we are looking for, only move forward if found
        my $fwd_valid = 0;
        my $rev_valid = 0;
        my $helical = 0;
        my $stepcyc = 0;
        my $stepsynth = 0;
        my $variable = 0;
        my $bimol_valid = 0;
        my $helicalfwd = 0;
        my $helicalrev = 0;
        my $variablefwd = 0;
        my $variablerev = 0;
        my $stepcycfwd = 0;
        my $stepcycrev = 0;
        my $stepsynthfwd = 0;
        my $stepsynthrev = 0;
        my $type = "yes";
        ## set hashes to zero for final size range if no molecule was found of that size
        my $short_size = [119..219];
        for my $shortsize (@{$short_size}) {
            if (!defined($unicyclized4_final_lengths{$shortsize})) {
                $unicyclized4_final_lengths{$shortsize} = 0;
            }
            if (!defined($linear4_final_lengths{$shortsize})) {
                $linear4_final_lengths{$shortsize} = 0;
            }
            if (!defined($linear3_final_lengths{$shortsize})) {
                $linear3_final_lengths{$shortsize} = 0;
            }
            if (!defined($bicyclized4_final_lengths{$shortsize})) {
                $bicyclized4_final_lengths{$shortsize} = 0;
            }
            if (!defined($bicyclized6_final_lengths{$shortsize})) {
                $bicyclized6_final_lengths{$shortsize} = 0;
            }
            if (!defined($bicyclized5_final_lengths{$shortsize})) {
                $bicyclized5_final_lengths{$shortsize} = 0;
            }
        }
        my $frag_size = [0..40];
        for my $fragsize (@{$frag_size}) {
            if (!defined($bicyclized5_final_lengths{$fragsize})) {
                $bicyclized5_final_lengths{$fragsize} = 0;
            }
        }
        ## set hashes to zero for full name if no molecule was found of that name
        ## this is not currently working, full_size and frag_full are not making hashes that I want them to
        my $step_sizes = [047,077,107];
        my $var_sizes = [00..30];
        my $hel_sizes = [00..10];
        my $full_size = [($step_sizes.$var_sizes.$hel_sizes)];
        for my $fullsize (@{$full_size}) {
            if (!defined($unicyclized4_full_final_lengths{$fullsize})) {
                $unicyclized4_full_final_lengths{$fullsize} = 0;
            }
            if (!defined($linear4_full_final_lengths{$fullsize})) {
                $linear4_full_final_lengths{$fullsize} = 0;
            }
            if (!defined($linear3_full_final_lengths{$fullsize})) {
                $linear3_full_final_lengths{$fullsize} = 0;
            }
            if (!defined($bicyclized4_full_final_lengths{$fullsize})) {
                $bicyclized4_full_final_lengths{$fullsize} = 0;
            }
            if (!defined($bicyclized6_full_final_lengths{$fullsize})) {
                $bicyclized6_full_final_lengths{$fullsize} = 0;
            }
            if (!defined($bicyclized5_full_final_lengths{$fullsize})) {
                $bicyclized5_full_final_lengths{$fullsize} = 0;
            }
        }
        my $frag_full = [($numbers{variable_fwd}.$numbers{helical_fwd})];
        for my $fragfull (@{$frag_full}) {
            if (!defined($bicyclized5_frag_final_lengths{$fragfull})) {
                $bicyclized5_frag_final_lengths{$fragfull} = 0;
            }
        }
        ## Here we look only at files that have cyclized, count them, and place the count into its own csv
        ## first check that four indices are observed & they are the four we expect to see for a unimolecular cyclization
        if ($observed_indices == 4 && $observe{stepsynth_fwd} > 0 && $observe{helical_fwd} > 0 &&
                $observe{stepcyc_fwd} > 0 && $observe{variable_fwd} > 0) {
            $fwd_valid = 1;
            $found_all_four++;
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir)  = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helical = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcyc = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variable = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynth = $name;
                }
            }
            ## I need a means of doing this that separates out into fwd/rev for bimolecular A-A
            my $final_size = $options{spacer} + $helical + $stepsynth + $variable;
            $comment .= "$count final size: $final_size ";
            ## Add a check for cyclized vs. not vs. unknown.
            ## I want to add a check for $info->{number} from stepsynth = $info->{number} from stepcyc. If this order & stepsynth number = stepcyc number $unicyclized = "yes"; else if this order & stepsynth number != stepcyc number $bimolcyclized = "yes"; then count out separately into csv files for unimolecular and bimolecular.
            if ($positions{stepcyc_fwd} < $positions{helical_fwd} &&
                    $positions{helical_fwd} < $positions{variable_fwd} &&
                    $positions{variable_fwd} < $positions{stepsynth_fwd} &&
                    $numbers{stepsynth_fwd} == $numbers{stepcyc_fwd}) {
                $type = "unimolecular";
                $found_four_uni++;
                if (!defined($unicyclized4_final_lengths{$final_size})) {
                    $unicyclized4_final_lengths{$final_size} = 1;
                } else {
                    $unicyclized4_final_lengths{$final_size}++;
                }
                if (!defined($unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
                ## if stepsynth & stepcyc do not match, but the order is still the same, this is a bimolecular A to B cyclization
            } elsif ($positions{stepcyc_fwd} < $positions{helical_fwd} &&
                         $positions{helical_fwd} < $positions{variable_fwd} &&
                         $positions{variable_fwd} < $positions{stepsynth_fwd} &&
                         $numbers{stepsynth_fwd} != $numbers{stepcyc_fwd}) {
                $type = "bimolecular-4";
                $found_four_bi++;
                if (!defined($bicyclized4_final_lengths{$final_size})) {
                    $bicyclized4_final_lengths{$final_size} = 1;
                } else {
                    $bicyclized4_final_lengths{$final_size}++;
                }
                if (!defined($bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
                ## if order is that of initial library we count up linear molecules
            } elsif ($positions{helical_fwd} < $positions{variable_fwd} &&
                         $positions{variable_fwd} < $positions{stepsynth_fwd} &&
                         $positions{stepsynth_fwd} < $positions{stepcyc_fwd}) {
                $type = "linear";
                $found_four_lin++;
                if (!defined($linear4_final_lengths{$final_size})) {
                    $linear4_final_lengths{$final_size} = 1;
                } else {
                    $linear4_final_lengths{$final_size}++;
                }
                if (!defined($linear4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $linear4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $linear4_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
            ## here we append the status of "cyclized" to the comment line, options for 4 hits are: unimolecular, bimolecular-4, linear, and unknown
            $comment .= "cyclized type: ${type} ";
        }
        ## now we look at the same thing as above, but for reverse molecules
        if ($observed_indices == 4 && $observe{stepsynth_rev} > 0 && $observe{helical_rev} > 0 &&
                $observe{stepcyc_rev} > 0 && $observe{variable_rev} > 0) {
            $rev_valid = 1;
            $found_all_four++;
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir)  = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'rev') {
                    $helical = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcyc = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variable = $name
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynth = $name
                }
            }
            my $final_size = $options{spacer} + $helical + $stepsynth + $variable;
            $comment .= "$count final size: $final_size ";
            ## Add a check for cyclized vs. not vs. unknown.
            if ($positions{stepcyc_rev} > $positions{helical_rev} &&
                    $positions{helical_rev} > $positions{variable_rev} &&
                    $positions{variable_rev} > $positions{stepsynth_rev} &&
                    $numbers{stepsynth_rev} == $numbers{stepcyc_rev}) {
                $type = "unimolecular";
                $found_four_uni++;
                if (!defined($unicyclized4_final_lengths{$final_size})) {
                    $unicyclized4_final_lengths{$final_size} = 1;
                } else {
                    $unicyclized4_final_lengths{$final_size}++;
                }
                if (!defined($unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
                ## if stepsynth & stepcyc do not match, but the order is still the same, this is a bimolecular A-B ligation
            } elsif ($positions{stepcyc_rev} > $positions{helical_rev} &&
                         $positions{helical_rev} > $positions{variable_rev} &&
                         $positions{variable_rev} > $positions{stepsynth_rev} &&
                         $numbers{stepsynth_rev} != $numbers{stepcyc_rev}) {
                $type = "bimolecular-4";
                $found_four_bi++;
                if (!defined($bicyclized4_final_lengths{$final_size})) {
                    $bicyclized4_final_lengths{$final_size} = 1;
                } else {
                    $bicyclized4_final_lengths{$final_size}++;
                }
                if (!defined($bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
            ## if order is that of initial library we count up linear molecules
            } elsif ($positions{helical_rev} > $positions{variable_rev} &&
                         $positions{variable_rev} > $positions{stepsynth_rev} &&
                         $positions{stepsynth_rev} > $positions{stepcyc_rev}) {
                $type = "linear4";
                $found_four_lin++;
                if (!defined($linear4_final_lengths{$final_size})) {
                    $linear4_final_lengths{$final_size} = 1;
                } else {
                    $linear4_final_lengths{$final_size}++;
                }
                if (!defined($linear4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $linear4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $linear4_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
                                ## here we append the status of "cyclized" to the comment line, options for 4 hits are: unimolecular, bimolecular-4, linear, and unknown
            $comment .= "cyclized type: ${type} ";
        }
        ## looking at triple index hits
        if ($observed_indices == 3 && $observe{helical_fwd} > 0 && $observe{variable_fwd} > 0 && $observe{stepsynth_fwd} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name;
                }
            }
            my $final_size = $options{spacer} + $helicalfwd + $stepsynthfwd + $variablefwd;
            $comment .= "$count final size: $final_size ";
            if ($positions{helical_fwd} < $positions{variable_fwd} &&
                $positions{variable_fwd} < $positions{stepsynth_fwd}) {
                $found_three_lin++;
                $type = "linear-3";
                if (!defined($linear3_final_lengths{$final_size})) {
                    $linear3_final_lengths{$final_size} = 1;
                } else {
                    $linear3_final_lengths{$final_size}++;
                }
                if (!defined($linear3_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $linear3_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $linear3_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown3hit";
                $found_three_unknown++;
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 3 && $observe{helical_rev} > 0 && $observe{variable_rev} > 0 && $observe{stepsynth_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name;
                }
            }
            my $final_size = $options{spacer} + $helicalrev + $stepsynthrev + $variablerev;
            $comment .= "$count final size: $final_size ";
            if ($positions{helical_rev} > $positions{variable_rev} &&
                $positions{variable_rev} > $positions{stepsynth_rev}) {
                $found_three_lin++;
                $type = "linear-3";
                if (!defined($linear3_final_lengths{$final_size})) {
                    $linear3_final_lengths{$final_size} = 1;
                } else {
                    $linear3_final_lengths{$final_size}++;
                }
                if (!defined($linear3_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $linear3_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $linear3_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown3hit";
                $found_three_unknown++;
            }
            $comment .= "type: ${type} ";
        }
        ## looking at single index hits
        if ($observed_indices == 1 && $observe{stepcyc_fwd} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcycfwd = $name;
                }
            }
            my $final_size = $stepcycfwd;
            $comment .= "$count final size: $final_size ";
            $found_one_lin++;
            $type = "linear-1";
            if (!defined($linear1_final_lengths{$final_size})) {
                $linear1_final_lengths{$final_size} = 1;
            } else {
                $linear1_final_lengths{$final_size}++;
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 1 && $observe{stepcyc_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcycrev = $name;
                }
            }
            my $final_size = $stepcycrev;
            $comment .= "$count final size: $final_size ";
            $found_one_lin++;
            $type = "linear-1";
            if (!defined($linear1_final_lengths{$final_size})) {
                $linear1_final_lengths{$final_size} = 1;
            } else {
                $linear1_final_lengths{$final_size}++;
            }
            $comment .= "type: ${type} ";
        }
        ## for 5 I am looking at a set with three and a fragment this can be:
        ## stepsynth_rev<variable_rev<helical_rev<helical_fwd<variable_fwd
        ## variable_rev<helical_rev<helical_fwd<variable_fwd<stepsynth_fwd
        ## the "size" for the helical+variable where there is no stepsynth will not include the spacer and will make these molecules range from "0-40" in size
        if ($observed_indices == 5 && $observe{stepsynth_fwd} > 0 && $observe{helical_fwd} > 0 &&
                $observe{variable_fwd} > 0 && $observe{helical_rev} > 0 && $observe{variable_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($positions, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name;
                }
            }
            my $final_bi5_fwd_size = $options{spacer} + $stepsynthfwd + $variablefwd + $helicalfwd;
            my $final_bi5_rev_size = $variablerev + $helicalrev;
            if ($positions{variable_rev} < $positions{helical_rev} &&
                    $positions{helical_rev} < $positions{helical_fwd} &&
                    $positions{helical_fwd} < $positions{variable_fwd} &&
                    $positions{variable_fwd} < $positions{stepsynth_fwd} ) {
                $found_five_bi++;
                $type = "bimolecular-5";
                $comment .= "$count final sizes: fwd: $final_bi5_fwd_size rev: $final_bi5_rev_size ";
                if (!defined($bicyclized5_final_lengths{final_bi5_fwd_size})) {
                    $bicyclized5_final_lengths{$final_bi5_fwd_size} = 1;
                } else {
                    $bicyclized5_final_lengths{final_bi5_fwd_size}++;
                }
                if (!defined($bicyclized5_final_lengths{final_bi5_rev_size})) {
                    $bicyclized5_final_lengths{$final_bi5_rev_size} = 1;
                } else {
                    $bicyclized5_final_lengths{final_bi5_rev_size}++;
                }
                if (!defined($bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
                if (!defined($bicyclized5_frag_final_lengths{$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $bicyclized5_frag_final_lengths{$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized5_frag_final_lengths{$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown5hit";
                $found_five_unknown++;
            }
            $comment .= "cyclized type: ${type} ";
        }
        if ($observed_indices == 5 && $observe{stepsynth_rev} > 0 && $observe{helical_fwd} > 0 &&
                $observe{variable_fwd} > 0 && $observe{helical_rev} > 0 && $observe{variable_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($positions, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name;
                }
            }
            my $final_bi5_fwd_size = $variablefwd + $helicalfwd;
            my $final_bi5_rev_size = $options{spacer} + $variablerev + $helicalrev + $stepsynthrev;
            if ($positions{stepsynth_rev} < $positions{variable_rev} &&
                    $positions{variable_rev} < $positions{helical_rev} &&
                    $positions{helical_rev} < $positions{helical_fwd} &&
                    $positions{helical_fwd} < $positions{variable_fwd}) {
                $found_five_bi++;
                $type = "bimolecular-5";
                $comment .= "$count final sizes: fwd: $final_bi5_fwd_size rev: $final_bi5_rev_size ";
                if (!defined($bicyclized5_final_lengths{final_bi5_fwd_size})) {
                    $bicyclized5_final_lengths{$final_bi5_fwd_size} = 1;
                } else {
                    $bicyclized5_final_lengths{final_bi5_fwd_size}++;
                }
                if (!defined($bicyclized5_final_lengths{final_bi5_rev_size})) {
                    $bicyclized5_final_lengths{$final_bi5_rev_size} = 1;
                } else {
                    $bicyclized5_final_lengths{final_bi5_rev_size}++;
                }
                if (!defined($bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
                if (!defined($bicyclized5_frag_final_lengths{$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $bicyclized5_frag_final_lengths{$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized5_frag_final_lengths{$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown5hit";
                $found_five_unknown++;
            }
            $comment .= "cyclized type: ${type} ";
        }
        ##
        if ($observed_indices == 6 && $observe{stepsynth_fwd} > 0 && $observe{helical_fwd} > 0 &&
                $observe{variable_fwd} > 0 && $observe{stepsynth_rev} > 0 && $observe{helical_rev} > 0 &&
                $observe{variable_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name
                }
            }
            my $final_bi6_fwd_size = $options{spacer} + $stepsynthfwd + $variablefwd + $helicalfwd;
            my $final_bi6_rev_size = $options{spacer} + $stepsynthrev + $variablerev + $helicalrev;
            if ($positions{stepsynth_fwd} > $positions{variable_fwd} &&
                    $positions{variable_fwd} > $positions{helical_fwd} &&
                    $positions{helical_fwd} > $positions{helical_rev} &&
                    $positions{helical_rev} > $positions{variable_rev} &&
                    $positions{variable_rev} > $positions{stepsynth_rev} ) {
                $found_six_bi++;
                $type = "bimolecular-6";
                $comment .= "$count final sizes: fwd: $final_bi6_fwd_size rev: $final_bi6_rev_size ";
                if (!defined($bicyclized6_final_lengths{$final_bi6_fwd_size})) {
                    $bicyclized6_final_lengths{$final_bi6_fwd_size} = 1;
                } else {
                    $bicyclized6_final_lengths{$final_bi6_fwd_size}++;
                }
                if (!defined(bicyclized6_final_lengths{$final_bi6_rev_size})) {
                    $bicyclized6_final_lengths{$final_bi6_rev_size} = 1;
                } else {
                    $bicyclized6_final_lengths{$final_bi6_rev_size}++;
                }
                if (!defined($bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}})) {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.$numbers{variable_rev}.$numbers{helical_rev}}++;
                }
                if (!defined($bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}})) {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.$numbers{variable_fwd}.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown6hit";
                $found_six_unknown++;
            }
            $comment .= "cyclized type: ${type} ";
        }
        ## Scan for bimolecular products from B-B ligations
        if ($observed_indices == 2 && $observe{stepcyc_fwd} > 0 && $observe{stepcyc_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcycfwd = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcycrev = $name;
                }
            }
            ## counting here just looks at the total size, as they can be separated using this in itself
            ## except we can't - 47+107 and 77+77 are equivilent. to separate this out, sizes need to be reported as what they are made of as well
            my $final_bi2_size = $stepcycfwd + $stepcycrev;
            if ($positions{stepcyc_fwd} < $positions{stepcyc_rev}) {
                $found_two_bi++;
                $type = "bimolecular-2";
                $comment .= "$count final size: $final_bi2_size ";
                if (!defined($bicyclized2_final_lengths{$final_bi2_size})) {
                    $bicyclized2_final_lengths{$final_bi2_size} = 1;
                } else {
                    $bicyclized2_final_lengths{$final_bi2_size}++;
                }
                if (!defined($bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.$numbers{stepcyc_rev}})) {
                    $bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.$numbers{stepcyc_rev}} = 1;
                } else {
                    $bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.$numbers{stepcyc_rev}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "cyclized type: ${type} ";
        } else {
            $comment .= "$count ";
        }
        ## append to the comment the number of observed indices.
        $comment .= " hits: ${observed_indices} ";
        ## Increment the singles->tens (& 11+) with the relevant observations.
        if ($observed_indices == 1) {
            $singles{sum}++;
            foreach my $k (keys %observe) {
                $singles{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 2) {
            $doubles{sum}++;
            foreach my $k (keys %observe) {
                $doubles{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 3) {
            $triples{sum}++;
            foreach my $k (keys %observe) {
                $triples{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 4) {
            $quads{sum}++;
            foreach my $k (keys %observe) {
                $quads{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 5) {
            $fives{sum}++;
            foreach my $k (keys %observe) {
                $fives{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 6) {
            $sixes{sum}++;
            foreach my $k (keys %observe) {
                $sixes{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 7) {
            $sevens{sum}++;
            foreach my $k (keys %observe) {
                $sevens{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 8) {
            $eights{sum}++;
            foreach my $k (keys %observe) {
                $eights{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 9) {
            $nines{sum}++;
            foreach my $k (keys %observe) {
                $nines{$k} += $observe{$k};
            }
        } elsif ($observed_indices == 10) {
            $tens{sum}++;
            foreach my $k (keys %observe) {
                $tens{$k} += $observe{$k};
            }
        } elsif ($observed_indices >= 11) {
            $elevenup{sum}++;
            foreach my $k (keys %observe) {
                $elevenup{$k} += $observe{$k};
            }
        }
        ## print STDOUT "$comment\n";
        my $fastq_string = qq"\@${id}
${sequence}
+${comment}
${qual}
";
        print $out $fastq_string;
    }
    return($data);
}
sub Sort_Dir {
    my %args = @_;
    my $cwd_dir = getcwd();
    my $searchdir = qq"$args{indir}";
    my $files = 0;
    unless ($searchdir =~ /^\//) {
        $searchdir = qq"${cwd_dir}/${searchdir}";
    }
    ## print "Searching: $searchdir  for files to read.\n";
    my @directory = ($searchdir);
    my @file_list = ();
    find(sub { push(@file_list, $File::Find::name)
                   if ($File::Find::name =~ /\.fastq\.gz/ and
                       $File::Find::name !~ /$args{outdir}/); }, @directory);
    my @approxes = ();
    foreach my $file (@file_list) {
        $files = $files++;
        next if ($file =~ /$options{outdir}/);
        ## This might be incorrect, CHECKME!
        $file =~ s/\/\.\//\//;
        my $approx = Sort_File_Approx(
            input => $file,
            outdir => $args{outdir},
            index_hash => $args{index_hash},
        );
        push(@approxes, $approx);
    }
    return(@approxes);
}

sub Read_indices {
    my %args = @_;

    my $indices = {
        total => {},
        unknown => {},
        possibilities => [],
    };

    $indices->{total} = {
        read => 0, written => 0, };
    $indices->{unknown} = {
        name => 'start_unknown', written => 0, };

    my $index_file = FileHandle->new("<$options{indices}");
    while (my $line = <$index_file>) {
        chomp $line;
        next if ($line =~ /^#/);
        next unless ($line =~ /^A|T|G|C|a|t|g|c/);

        my ($index_sequence, $phrase) = split(/\s+|\,|;/, $line);
        my ($name, $number, $direction) = split(/_/, $phrase);
        $indices->{$index_sequence} = {
            name => $name,
            direction => $direction,
            # I am formating each number to include three digits this is because my greatest size is three digits and I want everything to include all three (even if all zeros).
            # I think I can get away with not using this given how I have the numbering laid out. It is already formatted as I want.
            # number => sprintf("%03d",$number),
            number => $number,
            # original number input converted number to integer
            # number => int($number),
            total_observed => 0,
            ambiguous_observed => 0,
            unique_observed => 0
        };

        my @pos = @{$indices->{possibilities}};
        push(@pos, $index_sequence);
        $indices->{possibilities} = \@pos;
    }
    $index_file->close();
    return($indices);
}
## Here we create the overall summary file
sub End_Handler {
    print $log "Index used: $options{indices}\n";
    print $log "Input file: $options{input}\n";
    print $log "Date and time run $datestring\n";
    if ( ! defined $options{insertion} && ! defined $options{deletion} &&
            ! defined $options{substitution} ) {
        print $log "Direct match used";
    } else {
        print $log "String::Approx matching used: I$options{insertion},D$options{deletion},S$options{substitution}\n";
    }
    print $log "${observed_reads} reads were observed in total, of these:\n";
    foreach my $k (sort keys %observations) {
        if ($observations{$k} > 0 and $k ne 'sum') {
            ## Here we print how many times each individual index class was observed (regardless of if any other indices were found within that sequence). An index class is either the large variable region, small (helical) variable region, or synthesis variable region.
            print $log "     The read type: ${k} was observed: $observations{$k} times.\n";
        }
    }
    ## how many times single, double...eleven+ index sequences were found & how many times individual index classes were found within them.
    if ($singles{sum} > 0) {
        print $log "$singles{sum} single-index reads were observed, including:\n";
        foreach my $k (sort keys %singles) {
            if ($singles{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $singles{$k} times.\n";
            }
        }
    }
    if ($doubles{sum} > 0) {
        print $log "$doubles{sum} double-index reads were observed, including:\n";
        foreach my $k (sort keys %doubles) {
            if ($doubles{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $doubles{$k} times.\n";
            }
        }
    }
    if ($triples{sum} > 0) {
        print $log "$triples{sum} triple-index reads were observed, including:\n";
        foreach my $k (sort keys %triples) {
            if ($triples{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $triples{$k} times.\n";
            }
        }
    }
    if ($quads{sum} > 0) {
        print $log "$quads{sum} 4-index reads were observed, including:\n";
        foreach my $k (sort keys %quads) {
            if ($quads{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $quads{$k} times.\n";
            }
        }
    }
    if ($fives{sum} > 0) {
        print $log "$fives{sum} 5-index reads were observed, including:\n";
        foreach my $k (sort keys %fives) {
            if ($fives{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $fives{$k} times.\n";
            }
        }
    }
    if ($sixes{sum} > 0) {
        print $log "$sixes{sum} 6-index reads were observed, including:\n";
        foreach my $k (sort keys %sixes) {
            if ($sixes{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $sixes{$k} times.\n";
            }
        }
    }
    if ($sevens{sum} > 0) {
        print $log "$sevens{sum} 7-index reads were observed, including:\n";
        foreach my $k (sort keys %sevens) {
            if ($sevens{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $sevens{$k} times.\n";
            }
        }
    }
    if ($eights{sum} > 0) {
        print $log "$eights{sum} 8-index reads were observed, including:\n";
        foreach my $k (sort keys %eights) {
            if ($eights{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $eights{$k} times.\n";
            }
        }
    }
    if ($nines{sum} > 0) {
        print $log "$nines{sum} 9-index reads were observed, including:\n";
        foreach my $k (sort keys %nines) {
            if ($nines{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $nines{$k} times.\n";
            }
        }
    }
    if ($tens{sum} > 0) {
        print $log "$tens{sum} 10-index reads were observed, including:\n";
        foreach my $k (sort keys %tens) {
            if ($tens{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $tens{$k} times.\n";
            }
        }
    }
    if ($elevenup{sum} > 0) {
        print $log "$elevenup{sum} 11 or more-index reads were observed, including:\n";
        foreach my $k (sort keys %elevenup) {
            if ($elevenup{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $elevenup{$k} times.\n";
            }
        }
    }
    print $log "\n";
    print $log "${found_all_four} reads had a stepsynth+variable+helical+stepcyc:\n";
    print $log "${found_four_uni} reads had a stepsynth+variable+helical+stepcyc and cyclized:\n";
    foreach my $k (sort keys %unicyclized4_final_lengths) {
        print $log "    Size $k was found $unicyclized4_final_lengths{$k} times and cyclized.\n";
        if ($k ne "$options{spacer}") {
            print $unicyc_csv "$k,$unicyclized4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %unicyclized4_full_final_lengths) {
        print $unicycfull_csv "$k,$unicyclized4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_lin} reads had a stepcyc+stepsynth+variable+helical and not cyclized:\n";
    foreach my $k (sort keys %linear4_final_lengths) {
        print $log "    Size $k was found $linear4_final_lengths{$k} times and not cyclized.\n";
        if ($k ne "$options{spacer}") {
            print $fourhitlin_csv "$k,$linear4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %linear4_full_final_lengths) {
        print $fourhitlinfull_csv "$k,$linear4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_bi} reads had a biomolecular stepsynth+variable+helical+stepcyc:\n";
    foreach my $k (sort keys %bicyclized4_final_lengths) {
        print $log "    Size $k was found $bicyclized4_final_lengths{$k} times and cyclized (bimolecular).\n";
        if ($k ne "$options{spacer}") {
            print $bicyc4_csv "$k,$bicyclized4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %bicyclized4_full_final_lengths) {
        print $bicyc4full_csv "$k,$bicyclized4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_unknown} reads had four hits and were unknown\n";
    print $log "${found_three_lin} reads had three hits and were linear fragments stepsynth+variable+helical:\n";
    foreach my $k (sort keys %linear3_final_lengths) {
        print $log "    Size $k was found $linear3_final_lengths{$k} times.\n";
        if ($k ne "$options{spacer}") {
            print $threehitlin_csv "$k,$linear3_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %linear3_full_final_lengths) {
        print $threehitlinfull_csv "$k,$linear3_full_final_lengths{$k}\n";
    }
    print $log "${found_one_lin} reads had a single hit and were linear fragments:\n";
    foreach my $k (sort keys %linear1_final_lengths) {
        print $log "    Size $k was found $linear1_final_lengths{$k} times.\n";
        print $onehitlin_csv "$k,$linear1_final_lengths{$k}\n";
    }
    print $log "${found_five_bi} reads had five hits and were bimolecular\n";
    foreach my $k (sort keys %bicyclized5_final_lengths) {
        print $log "    Size $k was found $bicyclized5_final_lengths{$k} times and cyclized.\n";
        if ($k ne "$options{spacer}") {
            print $bicyc5_csv "$k,$bicyclized5_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %bicyclized5_full_final_lengths) {
        print $bicyc5full_csv "$k,$bicyclized5_full_final_lengths{$k}\n";
    }
    foreach my $k (sort keys %bicyclized5_frag_final_lengths) {
        print $bicyc5frag_csv "$k,$bicyclized5_frag_final_lengths{$k}\n";
    }
    print $log "${found_six_bi} reads had six hits and were bimolecular\n";
    foreach my $k (sort keys %bicyclized6_final_lengths) {
        print $log "    Size $k was found $bicyclized6_final_lengths{$k} times and cyclized.\n";
        if ($k ne "$options{spacer}") {
            print $bicyc6_csv "$k,$bicyclized6_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %bicyclized6_full_final_lengths) {
        print $bicyc6full_csv "$k,bicyclized6_full_final_lengths{$k}\n";
    }
    print $log "${found_six_unknown} reads had six hits and were unknown\n";
    print $log "${found_two_bi} reads had two hits and were bimolecular:\n";
    foreach my $k (sort keys %bicyclized2_final_lengths) {
        print $log "    Size $k was found $bicyclized2_final_lengths{$k} times and cyclized (bimolecular).\n";
        print $bicyc2_csv "$k,$bicyclized2_final_lengths{$k}\n";
    }
    foreach my $k (sort keys %bicyclized2_full_final_lengths) {
        print $bicyc2full_csv "$k,$bicyclized2_full_final_lengths{$k}\n";
    }
    $log->close();
    $unicyc_csv->close();
    $unicycfull_csv->close();
    $fourhitlin_csv->close();
    $fourhitlinfull_csv->close();
    $threehitlin_csv->close();
    $threehitlinfull_csv->close();
    $onehitlin_csv->close();
    $bicyc4_csv->close();
    $bicyc4full_csv->close();
    $bicyc6_csv->close();
    $bicyc6full_csv->close();
    $bicyc5_csv->close();
    $bicyc5full_csv->close();
    $bicyc5frag_csv->close();
    $bicyc2_csv->close();
    $bicyc2full_csv->close();
    $out->close();
    exit(0);
}


sub aindexes {
    my ($index, $params, $sequence) = @_;
    my @return = ();
    my $continue = 1;
    while ($continue) {
        my @ind = aindex($index, $params, ($sequence));
        my $res = $ind[0];
        if ($res == -1) {
            $continue = 0;
        } else {
            push(@return, $res);
            my $end_position = $res + length($index);
            $sequence = substr($sequence, $end_position);
        }
    }
    return(@return);
}
