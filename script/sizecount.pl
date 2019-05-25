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

  This script has a few useful options:
 --index  : The file defining the search strings and output files.
 --input  : Either a filename or directory containing the (relatively) raw
    fastq data.
 --outdir : Output directory, composed of two csv files, an output fastq file,
    and a summary file.
 --outfastq : Output fastq file, by default this is placed into the output
    directory. File extension not necessary.
 --substitution : String::Approx (Levenshtein) substitution distance, default
    is 0 (perfect match)
 --insertion: Levenshtein insertion distance
 --deletion: Levenshtein deletion distance
 --spacer : specify spacer length (default 72)
 --debug  : Print a bunch of debugging output
 All of the output scripts have default names, any individual name can be changed
 as an option, these are not listed here for brevity.

=head1 DESCRIPTION

  This was originally a demultiplexer for TNSeq data, which has a very specific and
  illumina-incompatible (bcl2fastq) format.  It was slightly reformatted to match
  the following query:
 
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
    onehitlincyc_csv => 'linear_onehit_cyc_lengths.csv',
    onehitlinsynth_csv => 'linear_onehit_synth_lengths.csv',
    bicyc4_csv => 'bimolecular4_ligated_lengths.csv',
    bicyc4full_csv => 'bimolecular4_ligated_lengths_full.csv',
    bicyc4varlib_csv => 'bimolecular4_ligated_variable_library_lengths.csv',
    bicyc4varlibfull_csv => 'bimolecular4_ligated_variable_library_lengths_full.csv',
    bicyc4steplib_csv => 'bimolecular4_ligated_step_library_lengths.csv',
    bicyc4steplibfull_csv => 'bimolecular4_ligated_step_library_lengths_full.csv',
    bicyc6_csv => 'bimolecular6_ligated_lengths.csv',
    bicyc6full_csv => 'bimolecular6_ligated_lengths_full.csv',
    bicyc5_csv => 'bimolecular5_ligated_lengths.csv',
    bicyc5full_csv => 'bimolecular5_ligated_lengths_full.csv',
    bicyc5frag_csv => 'bimolecular5_ligated_lengths_frag.csv',
    bicyc2_csv => 'bimolecular2_ligated_lengths.csv',
    bicyc2full_csv => 'bimolecular2_ligated_lengths_full.csv',
    varlibfull_csv => 'library_variable_full.csv',
    steplibfull_csv => 'library_step_full.csv',
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
my %sevenup = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## set up tags for counting and appropriate final lengths outputs for each match
my $datestring = localtime();
my $found_all_four = 0;
my $found_four_uni = 0;
my $found_four_lin = 0;
my $found_three_lin = 0;
my $found_one_lin_cyc = 0;
my $found_one_lin_synth = 0;
my $found_three_unknown = 0;
my $found_one_unknown = 0;
my $found_four_unknown = 0;
my %unicyclized4_final_lengths = ();
my %unicyclized4_full_final_lengths = ();
my %linear4_final_lengths = ();
my %linear4_full_final_lengths = ();
my %linear3_final_lengths = ();
my %linear3_full_final_lengths = ();
my %linear1_final_lengths_cyc = ();
my %linear1_final_lengths_synth = ();
my $found_four_bi = 0;
my %bicyclized4_final_lengths = ();
my %bicyclized4_full_final_lengths = ();
my $found_four_varlib_bi = 0;
my $found_four_steplib_bi = 0;
my %bicyclized4_var_lib_final_lengths = ();
my %bicyclized4_var_lib_full_final_lengths = ();
my %bicyclized4_step_lib_final_lengths = ();
my %bicyclized4_step_lib_full_final_lengths = ();
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
my $found_two_lib_var = 0;
my $found_two_lib_step = 0;
my %bicyclized2_final_lengths = ();
my %bicyclized2_full_final_lengths = ();
my %bi2_full_lib_lengths = ();
my %bi2_full_step_lengths = ();
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
    "onehitlincyc_csv:s" => \$options{onehitlincyc_csv},
    "onehitlinsynth_csv:s" => \$options{onehitlinsynth_csv},
    "bicyc4_csv:s" => \$options{bicyc4_csv},
    "bicyc4full_csv:s" => \$options{bicyc4full_csv},
    "bicyc4varlib_csv:s" => \$options{bicyc4varlib_csv},
    "bicyc4varlibfull_csv:s" => \$options{bicyc4varlibfull_csv},
    "bicyc4steplib_csv:s" => \$options{bicyc4steplib_csv},
    "bicyc4steplibfull_csv:s" => \$options{bicyc4steplibfull_csv},
    "bicyc6_csv:s" => \$options{bicyc6_csv},
    "bicyc6full_csv:s" => \$options{bicyc6full_csv},
    "bicyc5_csv:s" => \$options{bicyc5_csv},
    "bicyc5full_csv:s" => \$options{bicyc5full_csv},
    "bicyc5frag_csv:s" => \$options{bicyc5frag_csv},
    "bicyc2_csv:s" => \$options{bicyc2_csv},
    "bicyc2full_csv:s" => \$options{bicyc2full_csv},
    "varlibfull_csv:s" => \$options{varlibfull_csv},
    "steplibfull_csv:s" => \$options{steplibfull_csv},
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
my $bicyc4varlib_csv = new FileHandle(">$options{outdir}/$options{bicyc4varlib_csv}");
my $bicyc4varlibfull_csv = new FileHandle(">$options{outdir}/$options{bicyc4varlibfull_csv}");
my $bicyc4steplib_csv = new FileHandle(">$options{outdir}/$options{bicyc4steplib_csv}");
my $bicyc4steplibfull_csv = new FileHandle(">$options{outdir}/$options{bicyc4steplibfull_csv}");
my $bicyc6_csv = new FileHandle(">$options{outdir}/$options{bicyc6_csv}");
my $bicyc6full_csv = new FileHandle(">$options{outdir}/$options{bicyc6full_csv}");
my $bicyc2_csv = new FileHandle(">$options{outdir}/$options{bicyc2_csv}");
my $bicyc2full_csv = new FileHandle(">$options{outdir}/$options{bicyc2full_csv}");
my $bicyc5_csv = new FileHandle(">$options{outdir}/$options{bicyc5_csv}");
my $bicyc5full_csv = new FileHandle(">$options{outdir}/$options{bicyc5full_csv}");
my $bicyc5frag_csv = new FileHandle(">$options{outdir}/$options{bicyc5frag_csv}");
my $threehitlin_csv = new FileHandle(">$options{outdir}/$options{threehitlin_csv}");
my $threehitlinfull_csv = new FileHandle(">$options{outdir}/$options{threehitlinfull_csv}");
my $onehitlincyc_csv = new FileHandle(">$options{outdir}/$options{onehitlincyc_csv}");
my $onehitlinsynth_csv = new FileHandle(">$options{outdir}/$options{onehitlinsynth_csv}");
my $varlibfull_csv = new FileHandle(">$options{outdir}/$options{varlibfull_csv}");
my $steplibfull_csv = new FileHandle(">$options{outdir}/$options{steplibfull_csv}");

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

my $abs_input = File::Spec->rel2abs($options{input});
my $fastq_output = qq"$options{outdir}/$options{outfastq}.fastq.gz";
my $abs_output = File::Spec->rel2abs($fastq_output);
my $out;
if ($abs_input eq $abs_output) {
    die("The input file and output file are the same.");
} else {
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
End_Handler();

=item Sort_File_Approx

  This function should look through every sequence for a reasonable
  match to the available indices. If it gets some hits, record them.
  If they are ambiguous, this should return where the possibilities
  lie. If no match, return that.

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
            my @starts = ();
            if ( !defined $options{insertion} && !defined $options{deletion} &&
                    !defined $options{substitution} ) {
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
                    ## where "st" is position within read, info name
                    ## and number identify the index, and direction
                    ## identifies fwd or rev.
                    ## fwd & rev needed as adapter attachment to
                    ## flow cell is random & can occur in either
                    ## direction.
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
        my $frag_size = ['00','01','02','03','04','05','06','07','08','09',10..40];
        for my $fragsize (@{$frag_size}) {
            if (!defined($bicyclized5_final_lengths{"0".$fragsize})) {
                $bicyclized5_final_lengths{"0".$fragsize} = 0;
                # when pre-populating this hash I am doing so with
                # a leading zero such that the values are three digits,
                # this allows for them to end in the same hash as the
                # regular size library as the range of sizes here is
                # 000-040 and the regular size library is 119-219.
                # if a different size library is being used this
                # should be changed to a different hash entirely
            }
            if (!defined($bicyclized4_var_lib_final_lengths{"0".$fragsize})) {
                $bicyclized4_var_lib_final_lengths{"0".$fragsize} = 0;
            }
        }
        ## set hashes to zero for full name if no molecule was found of that name
        my $step_sizes = ['047','077','107'];
        my $var_sizes = ['00','01','02','03','04','05','06','07','08','09',10..30];
        my $hel_sizes = ['00','01','02','03','04','05','06','07','08','09','10'];
        my @full_sizes = ();
        foreach my $s (@{$step_sizes}) {
            foreach my $v (@{$var_sizes}) {
                foreach my $h (@{$hel_sizes}) {
                    my $entry = qq"${s},${v},${h}";
                    push (@full_sizes,$entry);
                }
            }
        }
        my @frag_full = ();
        foreach my $v (@{$var_sizes}) {
            foreach my $h (@{$hel_sizes}) {
                my $entry = qq"${v},${h}";
                push (@frag_full,$entry);
            }
        }
        foreach my $fullsize (@{full_sizes}) {
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
        for my $fragfull (@{frag_full}) {
            if (!defined($bicyclized5_frag_final_lengths{$fragfull})) {
                $bicyclized5_frag_final_lengths{$fragfull} = 0;
            }
            if (!defined($bicyclized4_var_lib_full_final_lengths{$fragfull})) {
                $bicyclized4_var_lib_full_final_lengths{$fragfull} = 0;
            }
            if (!defined($bi2_full_lib_lengths{$fragfull})) {
                $bi2_full_lib_lengths{$fragfull} = 0;
            }
        }
        ## start looking for matches & assigning them
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
            my $final_size = $options{spacer} + $helical + $stepsynth + $variable;
            $comment .= "final size: $final_size ";
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
                if (!defined($unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
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
                if (!defined($bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
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
                if (!defined($linear4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $linear4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $linear4_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
            $comment .= "type: ${type} ";
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
            $comment .= "final size: $final_size ";
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
                if (!defined($unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $unicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
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
                if (!defined($bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
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
                if (!defined($linear4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $linear4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $linear4_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
            $comment .= "type: ${type} ";
        }
        ## four hits that come purely from library fragments that did not
        ## have a live BstEII (synthesis junction) side
        if ($observed_indices == 4 && $observe{helical_fwd} > 0 && $observe{helical_rev} > 0 && $observe{variable_fwd} > 0 && $observe{variable_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name;
                }
            }
            my $var_fwd_size = $helicalfwd + $variablefwd;
            my $var_rev_size = $helicalrev + $variablerev;
            if ($positions{variable_rev} < $positions{helical_rev} && $positions{helical_rev} < $positions{helical_fwd} && $positions{helical_fwd} <$positions{variable_fwd}) {
                $type = "bicyc4varlib";
                $found_four_varlib_bi++;
                if ($var_fwd_size < 10) {
                    if (!defined($bicyclized4_var_lib_final_lengths{'00'.$var_fwd_size})) {
                        $bicyclized4_var_lib_final_lengths{'00'.$var_fwd_size} = 1;
                    } else {
                        $bicyclized4_var_lib_final_lengths{'00'.$var_fwd_size}++;
                    }
                } else {
                    if (!defined($bicyclized4_var_lib_final_lengths{'0'.$var_fwd_size})) {
                        $bicyclized4_var_lib_final_lengths{'0'.$var_fwd_size} = 1;
                    } else {
                        $bicyclized4_var_lib_final_lengths{'0'.$var_fwd_size}++;
                    }
                }
                if ($var_rev_size < 10) {
                    if (!defined($bicyclized4_var_lib_final_lengths{'00'.$var_rev_size})) {
                        $bicyclized4_var_lib_final_lengths{'00'.$var_rev_size} = 1;
                    } else {
                        $bicyclized4_var_lib_final_lengths{'00'.$var_rev_size}++;
                    }
                } else {
                    if (!defined($bicyclized4_var_lib_final_lengths{'0'.$var_rev_size})) {
                        $bicyclized4_var_lib_final_lengths{'0'.$var_rev_size} = 1;
                    } else {
                        $bicyclized4_var_lib_final_lengths{'0'.$var_rev_size}++;
                    }
                }
                if (!defined($bicyclized4_var_lib_full_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bicyclized4_var_lib_full_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized4_var_lib_full_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
                if (!defined($bicyclized4_var_lib_full_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $bicyclized4_var_lib_full_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized4_var_lib_full_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
                $comment .= "fwd size: $var_fwd_size rev size: $var_rev_size ";
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
            $comment .= "$type: ${type} ";
        }
        if ($observed_indices == 4 && $observe{stepsynth_fwd} > 0 && $observe{stepsynth_rev} > 0 && $observe{stepcyc_fwd} > 0 && $observe{stepcyc_fwd} > 0 && $numbers{stepsynth_fwd} == $numbers{stepcyc_fwd} && $numbers{stepsynth_rev} == $numbers{stepcyc_rev}) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcycfwd = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcycrev = $name;
                }
            }
            my $step_fwd_size = $stepsynthfwd;
            my $step_rev_size = $stepsynthrev;
            if ($positions{stepsynth_fwd} < $positions{stepcyc_fwd} && $positions{stepcyc_fwd} < $positions{stepcyc_rev} && $positions{stepcyc_rev} <$positions{stepsynth_rev}) {
                $type = "bicyc4steplib";
                $found_four_steplib_bi++;
                if (!defined($bicyclized4_step_lib_final_lengths{$step_fwd_size})) {
                    $bicyclized4_step_lib_final_lengths{$step_fwd_size} = 1;
                } else {
                    $bicyclized4_step_lib_final_lengths{$step_fwd_size}++;
                }
                if (!defined($bicyclized4_step_lib_final_lengths{$step_rev_size})) {
                    $bicyclized4_step_lib_final_lengths{$step_rev_size} = 1;
                } else {
                    $bicyclized4_step_lib_final_lengths{$step_rev_size}++;
                }
                if (!defined($bicyclized4_step_lib_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{stepsynth_rev}})) {
                    $bicyclized4_step_lib_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{stepsynth_rev}} = 1;
                } else {
                    $bicyclized4_step_lib_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{stepsynth_rev}}++;
                }
                $comment .= "fwd size: $step_fwd_size rev size: $step_rev_size ";
            } else {
                $type = "unknown4hit";
                $found_four_unknown++;
            }
            $comment .= "$type: ${type} ";
        }
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
            $comment .= "final size: $final_size ";
            if ($positions{helical_fwd} < $positions{variable_fwd} &&
                $positions{variable_fwd} < $positions{stepsynth_fwd}) {
                $found_three_lin++;
                $type = "linear-3";
                if (!defined($linear3_final_lengths{$final_size})) {
                    $linear3_final_lengths{$final_size} = 1;
                } else {
                    $linear3_final_lengths{$final_size}++;
                }
                if (!defined($linear3_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $linear3_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $linear3_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
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
            $comment .= "final size: $final_size ";
            if ($positions{helical_rev} > $positions{variable_rev} &&
                $positions{variable_rev} > $positions{stepsynth_rev}) {
                $found_three_lin++;
                $type = "linear-3";
                if (!defined($linear3_final_lengths{$final_size})) {
                    $linear3_final_lengths{$final_size} = 1;
                } else {
                    $linear3_final_lengths{$final_size}++;
                }
                if (!defined($linear3_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $linear3_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $linear3_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown3hit";
                $found_three_unknown++;
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 1 && $observe{stepcyc_fwd} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcycfwd = $name;
                }
            }
            my $final_size = $stepcycfwd;
            $comment .= "final size: $final_size ";
            $found_one_lin_cyc++;
            $type = "linear-1-cyc";
            if (!defined($linear1_final_lengths_cyc{$final_size})) {
                $linear1_final_lengths_cyc{$final_size} = 1;
            } else {
                $linear1_final_lengths_cyc{$final_size}++;
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
            $comment .= "final size: $final_size ";
            $found_one_lin_cyc++;
            $type = "linear-1-cyc";
            if (!defined($linear1_final_lengths_cyc{$final_size})) {
                $linear1_final_lengths_cyc{$final_size} = 1;
            } else {
                $linear1_final_lengths_cyc{$final_size}++;
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 1 && $observe{stepsynth_fwd} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name;
                }
            }
            my $final_size = $stepsynthfwd;
            $comment .= "final size: $final_size ";
            $found_one_lin_synth++;
            $type = "linear-1-synth";
            if (!defined($linear1_final_lengths_synth{$final_size})) {
                $linear1_final_lengths_synth{$final_size} = 1;
            } else {
                $linear1_final_lengths_synth{$final_size}++;
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 1 && $observe{stepsynth_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name;
                }
            }
            my $final_size = $stepsynthrev;
            $comment .= "final size: $final_size ";
            $found_one_lin_synth++;
            $type = "linear-1-synth";
            if (!defined($linear1_final_lengths_synth{$final_size})) {
                $linear1_final_lengths_synth{$final_size} = 1;
            } else {
                $linear1_final_lengths_synth{$final_size}++;
            }
            $comment .= "type: ${type} ";
        }
        ## for 5 I am looking at a set with three and a fragment this can be:
        ## stepsynth_rev < variable_rev < helical_rev < helical_fwd < variable_fwd
        ##                 variable_rev < helical_rev < helical_fwd < variable_fwd < stepsynth_fwd
        ## these should only exist in the event that a library variable fragment
        ## found a full molecule and multimerized
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
            ## the reverse size is variable plus helical, but this
            ## gives a minimum "size" of "0" up to "40" - this is
            ## not technically the size, the size is actually 072-112
            # I have chosen to leave this as the size shown because it
            # reflects a fragment and I want it to be clear that these
            # fragments should not be capable of cyclization., in theory
            # only appearing in a dead linear multimer (so the size is
            # less relevant so much as knowing that a molecule participated
            # in generation of a multimer) and which particular molecule it was.
            my $final_bi5_rev_size = $variablerev + $helicalrev;
            if ($positions{variable_rev} < $positions{helical_rev} &&
                    $positions{helical_rev} < $positions{helical_fwd} &&
                    $positions{helical_fwd} < $positions{variable_fwd} &&
                    $positions{variable_fwd} < $positions{stepsynth_fwd} ) {
                $found_five_bi++;
                $type = "bimolecular-5";
                $comment .= "fwd size: $final_bi5_fwd_size rev size: $final_bi5_rev_size ";
                if (!defined($bicyclized5_final_lengths{$final_bi5_fwd_size})) {
                    $bicyclized5_final_lengths{$final_bi5_fwd_size} = 1;
                } else {
                    $bicyclized5_final_lengths{$final_bi5_fwd_size}++;
                }
                # for the fragmented size I need it to have three
                # digits like the rest of the library, for values
                # 0..9 prepend two leading 0s, else 10..40 prepend
                # one leading zero.
                if ($final_bi5_rev_size < 10) {
                    if (!defined($bicyclized5_final_lengths{'00'.$final_bi5_rev_size})) {
                        $bicyclized5_final_lengths{'00'.$final_bi5_rev_size} = 1;
                    } else {
                        $bicyclized5_final_lengths{'00'.$final_bi5_rev_size}++;
                    }
                } else {
                    if (!defined($bicyclized5_final_lengths{'0'.$final_bi5_rev_size})) {
                        $bicyclized5_final_lengths{'0'.$final_bi5_rev_size} = 1;
                    } else {
                        $bicyclized5_final_lengths{'0'.$final_bi5_rev_size}++;
                    }
                }
                if (!defined($bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
                if (!defined($bicyclized5_frag_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $bicyclized5_frag_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized5_frag_final_lengths{$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
            } else {
                $type = "unknown5hit";
                $found_five_unknown++;
            }
            $comment .= "type: ${type} ";
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
                $comment .= "fwd size: $final_bi5_fwd_size rev size: $final_bi5_rev_size ";
                if ($final_bi5_fwd_size < 10) {
                    if (!defined($bicyclized5_final_lengths{'00'.$final_bi5_fwd_size})) {
                        $bicyclized5_final_lengths{'00'.$final_bi5_fwd_size} = 1;
                    } else {
                        $bicyclized5_final_lengths{'00'.$final_bi5_fwd_size}++;
                    }
                } else {
                    if (!defined($bicyclized5_final_lengths{'0'.$final_bi5_fwd_size})) {
                        $bicyclized5_final_lengths{'0'.$final_bi5_fwd_size} = 1;
                    } else {
                        $bicyclized5_final_lengths{'0'.$final_bi5_fwd_size}++;
                    }
                }
                if (!defined($bicyclized5_final_lengths{$final_bi5_rev_size})) {
                    $bicyclized5_final_lengths{$final_bi5_rev_size} = 1;
                } else {
                    $bicyclized5_final_lengths{$final_bi5_rev_size}++;
                }
                if (!defined($bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized5_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
                if (!defined($bicyclized5_frag_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bicyclized5_frag_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized5_frag_final_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown5hit";
                $found_five_unknown++;
            }
            $comment .= "type: ${type} ";
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
                $comment .= "fwd size: $final_bi6_fwd_size rev size: $final_bi6_rev_size ";
                if (!defined($bicyclized6_final_lengths{$final_bi6_fwd_size})) {
                    $bicyclized6_final_lengths{$final_bi6_fwd_size} = 1;
                } else {
                    $bicyclized6_final_lengths{$final_bi6_fwd_size}++;
                }
                if (!defined($bicyclized6_final_lengths{$final_bi6_rev_size})) {
                    $bicyclized6_final_lengths{$final_bi6_rev_size} = 1;
                } else {
                    $bicyclized6_final_lengths{$final_bi6_rev_size}++;
                }
                if (!defined($bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}})) {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}} = 1;
                } else {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_rev}.','.$numbers{variable_rev}.','.$numbers{helical_rev}}++;
                }
                if (!defined($bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bicyclized6_full_final_lengths{$numbers{stepsynth_fwd}.','.$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown6hit";
                $found_six_unknown++;
            }
            $comment .= "type: ${type} ";
        }
        ##
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
            my $final_bi2_size = $stepcycfwd + $stepcycrev;
            if ($positions{stepcyc_fwd} < $positions{stepcyc_rev}) {
                $found_two_bi++;
                $type = "bimolecular-2";
                $comment .= "final size: $final_bi2_size ";
                # truthfully the final size here is useless for
                # our current library size given the separated
                # component counting below. I am leaving it
                # in case it becomes useful for future projects
                if (!defined($bicyclized2_final_lengths{$final_bi2_size})) {
                    $bicyclized2_final_lengths{$final_bi2_size} = 1;
                } else {
                    $bicyclized2_final_lengths{$final_bi2_size}++;
                }
                if (!defined($bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.','.$numbers{stepcyc_rev}})) {
                    $bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.','.$numbers{stepcyc_rev}} = 1;
                } else {
                    $bicyclized2_full_final_lengths{$numbers{stepcyc_fwd}.','.$numbers{stepcyc_rev}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 2 && $observe{stepcyc_fwd} > 0 && $observe{stepsynth_fwd} > 0 && $numbers{stepcyc_fwd} == $numbers{stepsynth_fwd} ) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'fwd') {
                    $stepcycfwd = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $stepsynthfwd = $name;
                }
            }
            my $final_step_component = $stepcycfwd;
            if ($positions{stepcyc_fwd} > $positions{stepsynth_fwd}) {
                $found_two_lib_step++;
                $type = "Step-Component-2";
                $comment .= "component size: $final_step_component ";
                if (!defined($bi2_full_step_lengths{$numbers{stepcyc_fwd}})) {
                    $bi2_full_step_lengths{$numbers{stepcyc_fwd}} = 1;
                } else {
                    $bi2_full_step_lengths{$numbers{stepcyc_fwd}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 2 && $observe{stepcyc_rev} > 0 && $observe{stepsynth_rev} > 0 && $numbers{stepcyc_rev} == $numbers{stepsynth_rev} ) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcycrev = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynthrev = $name;
                }
            }
            my $final_step_component = $stepcycrev;
            if ($positions{stepcyc_rev} > $positions{stepsynth_rev}) {
                $found_two_lib_step++;
                $type = "Step-Component-2";
                $comment .= "component size: $final_step_component ";
                if (!defined($bi2_full_step_lengths{$numbers{stepcyc_rev}})) {
                    $bi2_full_step_lengths{$numbers{stepcyc_rev}} = 1;
                } else {
                    $bi2_full_step_lengths{$numbers{stepcyc_rev}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 2 && $observe{helical_fwd} > 0 && $observe{variable_fwd} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'fwd') {
                    $helicalfwd = $name;
                } elsif ($piece eq 'variable' && $dir eq 'fwd') {
                    $variablefwd = $name;
                }
            }
            my $final_bi2_lib_size = $helicalfwd + $variablefwd;
            if ($positions{helical_fwd} < $positions{variable_fwd}) {
                $found_two_lib_var++;
                $type = "Variable-Component-2";
                $comment .= "component size: $final_bi2_lib_size ";
                if (!defined($bi2_full_lib_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}})) {
                    $bi2_full_lib_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}} = 1;
                } else {
                    $bi2_full_lib_lengths{$numbers{variable_fwd}.','.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "type: ${type} ";
        }
        if ($observed_indices == 2 && $observe{helical_rev} > 0 && $observe{variable_rev} > 0) {
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir) = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'rev') {
                    $helicalrev = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variablerev = $name;
                }
            }
            my $final_bi2_lib_size = $helicalrev + $variablerev;
            if ($positions{helical_rev} > $positions{variable_rev}) {
                $found_two_lib_var++;
                $type = "Variable-Component-2";
                $comment .= "component size: $final_bi2_lib_size ";
                if (!defined($bi2_full_lib_lengths{$numbers{variable_rev}.','.$numbers{helical_fwd}})) {
                    $bi2_full_lib_lengths{$numbers{variable_rev}.','.$numbers{helical_fwd}} = 1;
                } else{
                    $bi2_full_lib_lengths{$numbers{variable_rev}.','.$numbers{helical_fwd}}++;
                }
            } else {
                $type = "unknown2hit";
            }
            $comment .= "type: ${type} ";
        }
        $comment .= "Count: $count  hits: ${observed_indices} ";
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
        } elsif ($observed_indices >= 7) {
            $sevenup{sum}++;
            foreach my $k (keys %observe) {
                $sevenup{$k} += $observe{$k};
            }
        }
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
    my @directory = ($searchdir);
    my @file_list = ();
    find(sub { push(@file_list, $File::Find::name)
                   if ($File::Find::name =~ /\.fastq\.gz/ and
                       $File::Find::name !~ /$args{outdir}/); }, @directory);
    my @approxes = ();
    foreach my $file (@file_list) {
        $files = $files++;
        next if ($file =~ /$options{outdir}/);
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
            # I am formating each number to include three digits
            # this is because my greatest size is three digits and
            # I want everything to include all three (even if all zeros).
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
    print $log "Date and time run: $datestring\n";
    if ( !defined $options{insertion} && !defined $options{deletion} &&
            !defined $options{substitution} ) {
        print $log "Direct match used";
    } else {
        print $log "String::Approx matching used: I$options{insertion},D$options{deletion},S$options{substitution}\n";
    }
    print $log "Summary of index hits and full matches below:\n";
    print $log "${observed_reads} reads were observed in total, of these:\n";
    foreach my $k (sort keys %observations) {
        if ($observations{$k} > 0 and $k ne 'sum') {
            print $log "     The read type: ${k} was observed: $observations{$k} times.\n";
        }
    }
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
    if ($sevenup{sum} > 0) {
        print $log "$sevenup{sum} 7 or more-index reads were observed, including:\n";
        foreach my $k (sort keys %sevenup) {
            if ($sevenup{$k} > 0 and $k ne 'sum') {
                print $log "     The read type: ${k} was observed: $sevenup{$k} times.\n";
            }
        }
    }
    print $log "\n";
    print $log "${found_one_lin_cyc} reads had a single hit and were stepcyc linear fragments\n";
    foreach my $k (sort keys %linear1_final_lengths_cyc) {
        print $onehitlincyc_csv "$k,$linear1_final_lengths_cyc{$k}\n";
    }
    print $log "${found_one_lin_synth} reads had a single hit and were stepsynth linear fragments\n";
    foreach my $k (sort keys %linear1_final_lengths_synth) {
        print $onehitlinsynth_csv "$k,$linear1_final_lengths_synth{$k}\n";
    }
    print $log "\n";
    print $log "${found_two_bi} reads had two hits and were bimolecular B-B\n";
    foreach my $k (sort keys %bicyclized2_final_lengths) {
        print $bicyc2_csv "$k,$bicyclized2_final_lengths{$k}\n";
    }
    foreach my $k (sort keys %bicyclized2_full_final_lengths) {
        print $bicyc2full_csv "$k,$bicyclized2_full_final_lengths{$k}\n";
    }
    print $log "${found_two_lib_var} reads had two hits and were library variable fragments\n";
    foreach my $k (sort keys %bi2_full_lib_lengths) {
        print $varlibfull_csv "$k,$bi2_full_lib_lengths{$k}\n";
    }
    print $log "${found_two_lib_step} reads had two hits and were library step fragments\n";
    foreach my $k (sort keys %bi2_full_step_lengths) {
        print $steplibfull_csv "$k,$bi2_full_step_lengths{$k}\n";
    }
    print $log "\n";
    print $log "${found_three_lin} reads had three hits and were linear fragments stepsynth+variable+helical\n";
    foreach my $k (sort keys %linear3_final_lengths) {
        if ($k ne "$options{spacer}") {
            print $threehitlin_csv "$k,$linear3_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %linear3_full_final_lengths) {
        print $threehitlinfull_csv "$k,$linear3_full_final_lengths{$k}\n";
    }
    print $log "\n";
    print $log "${found_four_uni} reads had four hits and were cyclized unimolecular\n";
    foreach my $k (sort keys %unicyclized4_final_lengths) {
        if ($k ne "$options{spacer}") {
            print $unicyc_csv "$k,$unicyclized4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %unicyclized4_full_final_lengths) {
        print $unicycfull_csv "$k,$unicyclized4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_lin} reads had four hits and were linear library molecules\n";
    foreach my $k (sort keys %linear4_final_lengths) {
        if ($k ne "$options{spacer}") {
            print $fourhitlin_csv "$k,$linear4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %linear4_full_final_lengths) {
        print $fourhitlinfull_csv "$k,$linear4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_bi} reads had four hits and were biomolecular A-B\n";
    foreach my $k (sort keys %bicyclized4_final_lengths) {
        if ($k ne "$options{spacer}") {
            print $bicyc4_csv "$k,$bicyclized4_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %bicyclized4_full_final_lengths) {
        print $bicyc4full_csv "$k,$bicyclized4_full_final_lengths{$k}\n";
    }
    print $log "${found_four_varlib_bi} reads had four hits and were bimolecular A-A library fragments\n";
    foreach my $k (sort keys %bicyclized4_var_lib_final_lengths) {
        print $bicyc4varlib_csv "$k,$bicyclized4_var_lib_final_lengths{$k}\n";
    }
    foreach my $k (sort keys %bicyclized4_var_lib_full_final_lengths) {
        print $bicyc4varlibfull_csv "$k,$bicyclized4_var_lib_full_final_lengths{$k}\n";
    }
    print $log "${found_four_steplib_bi} reads had four hits and were bimolecular B-B library fragments\n";
    foreach my $k (sort keys %bicyclized4_step_lib_final_lengths) {
        print $bicyc4steplib_csv "$k,$bicyclized4_step_lib_final_lengths{$k}\n";
    }
    foreach my $k (sort keys %bicyclized4_step_lib_full_final_lengths) {
        print $bicyc4steplibfull_csv "$k,$bicyclized4_step_lib_full_final_lengths{$k}\n";
    }
    print $log "\n";
    print $log "${found_five_bi} reads had five hits and were bimolecular A-A\n";
    foreach my $k (sort keys %bicyclized5_final_lengths) {
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
    print $log "\n";
    print $log "${found_six_bi} reads had six hits and were bimolecular A-A\n";
    foreach my $k (sort keys %bicyclized6_final_lengths) {
        if ($k ne "$options{spacer}") {
            print $bicyc6_csv "$k,$bicyclized6_final_lengths{$k}\n";
        }
    }
    foreach my $k (sort keys %bicyclized6_full_final_lengths) {
        print $bicyc6full_csv "$k,$bicyclized6_full_final_lengths{$k}\n";
    }
    $log->close();
    $unicyc_csv->close();
    $unicycfull_csv->close();
    $fourhitlin_csv->close();
    $fourhitlinfull_csv->close();
    $threehitlin_csv->close();
    $threehitlinfull_csv->close();
    $onehitlincyc_csv->close();
    $onehitlinsynth_csv->close();
    $bicyc4_csv->close();
    $bicyc4full_csv->close();
    $bicyc4varlib_csv->close();
    $bicyc4varlibfull_csv->close();
    $bicyc4steplib_csv->close();
    $bicyc4steplibfull_csv->close();
    $bicyc6_csv->close();
    $bicyc6full_csv->close();
    $bicyc5_csv->close();
    $bicyc5full_csv->close();
    $bicyc5frag_csv->close();
    $bicyc2_csv->close();
    $bicyc2full_csv->close();
    $varlibfull_csv->close();
    $steplibfull_csv->close();
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
