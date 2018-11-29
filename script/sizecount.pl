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

  sizes.pl - A writeall script to sort sequences by index and infer the sizes of the template DNA.

=head1 SYNOPSIS

  This script only has a few useful options:
 --index  : The file defining the search strings and output files.
 --input  : Either a filename or directory containing the (relatively) raw fastq data.
 --outdir : Output directory, composed of two csv files, an output fastq file, and a summary file.
 --outfastq : Output fastq file, by default this is placed into the output directory
 --debug  : Print a bunch of debugging output?

=head1 DESCRIPTION

  This was originally a demultiplexer for TNSeq data, which has a very specific and
  illumina-incompatible (bcl2fastq) format.  It was slightly reformatted to match Jason Hustedt and Jason Kahn's query which is as follows:

  1. The primary goal is to observe how many of every size DNA template was incorporated into
     a sequencing library.  The implicit inference: when (DNA-size % ~10), the energetic
     requirements for a DNA to loop back on itself decrease, increasing the likelihood that
     the molecule can circularize and therefore end up in the sequencing library.  Ergo, more
     counts of a particular size is inversely related to the energy requirement.
  2. The input molecules were designed so that each size is distinguishable by its sequence.
  3. One index is used to distinguish size-ranges (something like: 0-40, 41-70, 71-100, 100+)
  4. The important caveat: inter-molecular ligation is in a fight with intra-molecular ligation.
     a. Ergo, we need to be aware of multiple indices for two reasons; first to tell the size
        range, second because we need to be able to observe when size-a linearly ligates
        to size-b and when size-a bimolecularly cyclizes with size-b.

  Thus Trey's original implementation is useful but not ideal and has been changed to do the
  following:

  1. Rather than just dump the sequences into separate files, use the comment field
     to print the observed matches.
  2. When it does decide to write output files, it should have some logic to know the range(s)
     as well as the specific size(s) observed. The index file has been rewriten to include all relevant information within a single index file (for each diferent index type).
  3. Output a text file describing the number of observations of each-sized molecule that is cyclized or not cyclized, when observed with the appropriate number of hits for a single cycliztion event.
 
 The current implimentation (2018-11-13) needs to be updated to do the following:
 
 1. Verify that a sequence with three hits has the appropriate matching index before and after the cyclization site (this will tell us that it was likely unimolecular).
    a. add index sequences to allow for (1) --done?
    b. add logic to check
 2. Look at molecules with multiple hits (more than 3) and determine if they are bimolecular reactions.
    a. Does it matter that size-a ligated to size-b, or just that a molecule of each size participated in a bimolecular reaction?
 3. Separate out hits based not only on size (already done) and cyclization status (already done) but which molecule they came from. Ideally something like:
     47-30-10,165,2000
     Where the components say that the three index sizes obtained (47, 30, and 10) built a molecule of size 165, and it was found 2000 times. Currently this is ignored and the count is solely based on size, such that 47-30-10, 77-00-10, and 77-10-00 (and all iterations 77-01-09...77-09-01) are all included into the count of 165. While it is important that this total size of DNA cyclized (or did not cyclize), I cannot currently remove bias from library input if I don't know which of the individual molecules contirbuted to the overall count.
     Adding this logic will be complicated, but the initial searching of sequences and matching them to index files can already be done, at this point, we are looking only at the data output in the comment field. Should this be a separate program to run this type of logic?

=cut

#options from Getopt::Long; defaults
#we have added csv files for cyclized and not cyclized counts
my %options = (
    debug => 1,
    indices => 'proposed_index.txt',
    input => 'test.fastq.gz',
    outdir => 'output',
    outfastq => 'out.fastq.gz',
    summary => 'summary.txt',
    cyc_csv => 'cyclized_lengths.csv',
    not_csv => 'not_cyclized_lengths.csv',
    spacer => 72,
);
## This is a hash counting off how many times _every_ index is observed.
my $observed_reads = 0;
my %observations = (
    helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when one index is observed.
my %singles = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when two were found.
my %doubles = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when three were found.
my %triples = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when four were found.
my %quads = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when five were found.
my %fives = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This counts how often each index is observed when six were found.
my %sixes = (
    sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
    helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## Added to count for sevens (although in theory they should not exist)
my %sevens = (
sum => 0, helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0,);
## This lists how many of cyclized final lenghts, not cyclized final lengths, and all three were found.
my $found_all_four = 0;
my %cyclized_final_lengths = ();
my %notcyclized_final_lengths = ();

=item Getopt::Long invocation

  Invoke the command 'perldoc Getopt::Long' to learn about the many ways one may
  Pass command line options to your program.

=cut

my $opt_result = GetOptions(
    "debug:i" => \$options{debug},
    "spacer:i" => \$options{spacer},
    "indices:s" => \$options{indices},
    "input:s" => \$options{input},
    "outdir:s" => \$options{outdir},
    "summary:s" => \$options{summary},
    "cyc_csv:s" => \$options{cyc_csv},
    "not_csv:s" => \$options{not_csv},
    "outfastq:s" => \$options{outfastq},
);
my $log = new FileHandle(">$options{outdir}/$options{summary}");
my $cyc_csv = new FileHandle(">$options{outdir}/$options{cyc_csv}");
my $not_csv = new FileHandle(">$options{outdir}/$options{not_csv}");

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
my $fastq_output = qq"$options{outdir}/$options{outfastq}";
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

    ## Does FileHandle work here?
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
            helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0, unknown => 0);
        my %positions = (
            helical_fwd => 0, stepcyc_fwd => 0, variable_fwd => 0, stepsynth_fwd => 0,
            helical_rev => 0, stepcyc_rev => 0, variable_rev => 0, stepsynth_rev => 0, unknown => 0);
        my $observed_indices = 0;
        foreach my $index (@index_list) {
            my $info = $data->{$index};
            if ($sequence =~ m/$index/) {
                my @starts = @-;
                $found++;
                $found_id = $index;
                for my $st (@starts) {
                    ## old comment format
                    ## $comment .= "$st:$seqlen:$info->{name}:$info->{number}:$info->{direction} ";
                    $comment .= "$st:$info->{name}:$info->{number}:$info->{direction} ";
                    ## where "st" is position within read, info name and number identify the index, and direction identifies fwd or rev
                    if ($info->{name} eq 'helical' && $info->{direction} eq 'fwd') {
                        $observations{helical_fwd}++;
                        $observe{helical_fwd}++;
                        $positions{helical_fwd} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'stepcyc' && $info->{direction} eq 'fwd') {
                        $observations{stepcyc_fwd}++;
                        $observe{stepcyc_fwd}++;
                        $positions{stepcyc_fwd} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'stepsynth' && $info->{direction} eq 'fwd') {
                        $observations{stepsynth_fwd}++;
                        $observe{stepsynth_fwd}++;
                        $positions{stepsynth_fwd} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'variable' && $info->{direction} eq 'fwd') {
                        $observations{variable_fwd}++;
                        $observe{variable_fwd}++;
                        $positions{variable_fwd} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'helical' && $info->{direction} eq 'rev') {
                        $observations{helical_rev}++;
                        $observe{helical_rev}++;
                        $positions{helical_rev} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'stepcyc' && $info->{direction} eq 'rev') {
                        $observations{stepcyc_rev}++;
                        $observe{stepcyc_rev}++;
                        $positions{stepcyc_rev} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'stepsynth' && $info->{direction} eq 'rev') {
                        $observations{stepsynth_rev}++;
                        $observe{stepsynth_rev}++;
                        $positions{stepsynth_rev} = $st;
                        $observed_indices++;
                    } elsif ($info->{name} eq 'variable' && $info->{direction} eq 'rev') {
                        $observations{variable_rev}++;
                        $observe{variable_rev}++;
                        $positions{variable_rev} = $st;
                        $observed_indices++;
                    }
                }
            }
        } ## End each index
      ## for debugging print STDOUT - this is commented out for now.
        ## print STDOUT "$comment ";
        my $fwd_valid = 0;
        my $rev_valid = 0;
        my $helical = 0;
        my $stepcyc = 0;
      my $stepsynth = 0;
        my $variable = 0;
      ## Here we look only at files that have cyclized, count them, and place the count into its own csv
        my $cyclized = "yes";
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
                    $variable = $name
                } elsif ($piece eq 'stepsynth' && $dir eq 'fwd') {
                    $variable = $name
                }
            }
            my $final_size = $options{spacer} + $helical + $stepcyc + $variable;
            $comment .= "$count final size: $final_size ";
            ## Add a check for cyclized vs. not vs. unknown.
            if ($positions{stepcyc_fwd} < $positions{helical_fwd} &&
                    $positions{helical_fwd} < $positions{variable_fwd} &&
                    $positions{variable_fwd} < $positions{stepsynth_fwd}) {
                $cyclized = "yes";
                if (!defined($cyclized_final_lengths{$final_size})) {
                    $cyclized_final_lengths{$final_size} = 1;
                } else {
                    $cyclized_final_lengths{$final_size}++;
                }
            } elsif ($positions{helical_fwd} < $positions{variable_fwd} &&
                    $positions{variable_fwd} < $positions{stepsynth_fwd} &&
                    $positions{stepsynth_fwd} < $positions{stepcyc_fwd}) {
                $cyclized = "no";
                if (!defined($notcyclized_final_lengths{$final_size})) {
                    $notcyclized_final_lengths{$final_size} = 1;
                } else {
                    $notcyclized_final_lengths{$final_size}++;
                }
            } else {
                $cyclized = "unknown";
            }
            $comment .= "cyclized: ${cyclized} ";
        } elsif ($observed_indices == 4 && $observe{helical_rev} > 0 &&
      $observe{stepcyc_rev} > 0 && $observe{variable_rev} > 0 && $observe{stepsynth_rev} > 0) {
            $found_all_four++;
            $rev_valid = 1;
            my @pieces = split(/\s+/, $comment);
            for my $p (@pieces) {
                my ($position, $piece, $name, $dir)  = split(/:/, $p);
                if ($piece eq 'helical' && $dir eq 'rev') {
                    $helical = $name;
                } elsif ($piece eq 'stepcyc' && $dir eq 'rev') {
                    $stepcyc = $name;
                } elsif ($piece eq 'stepsynth' && $dir eq 'rev') {
                    $stepsynth = $name;
                } elsif ($piece eq 'variable' && $dir eq 'rev') {
                    $variable = $name;
                }
            }
            my $final_size = $options{spacer} + $helical + $stepcyc + $variable;
            $comment .= "$count final size: $final_size ";
            ## Add a check for cyclized vs. not vs. unknown.
            if ($positions{stepcyc_rev} > $positions{helical_rev} &&
          $positions{helical_rev} > $positions{variable_rev} && $positions{variable_rev} > $positions{stepsynth_rev}) {
                $cyclized = "yes";
                if (!defined($cyclized_final_lengths{$final_size})) {
                    $cyclized_final_lengths{$final_size} = 1;
                } else {
                    $cyclized_final_lengths{$final_size}++;
                }
            } elsif ($positions{helical_rev} > $positions{variable_rev} &&
          $positions{variable_rev} > $positions{stepsynth_rev} && $positions{stepsynth_rev} > $positions{stepcyc_rev}) {
                $cyclized = "no";
                if (!defined($notcyclized_final_lengths{$final_size})) {
                    $notcyclized_final_lengths{$final_size} = 1;
                } else {
                    $notcyclized_final_lengths{$final_size}++;
                }
            } else {
                $cyclized = "unknown";
            }
            $comment .= "cyclized: ${cyclized} ";

        } else {
            $comment .= "$count ";
        }

        ## append to the comment the number observed indices.
        $comment .= " hits: ${observed_indices} ";
        ## Increment the singles->sevens with the relevant observations.
        if ($observed_indices == 1) {
            $singles{sum}++;
            foreach my $k (keys %observe) { $singles{$k} += $observe{$k}; }
        } elsif ($observed_indices == 2) {
            $doubles{sum}++;
            foreach my $k (keys %observe) { $doubles{$k} += $observe{$k}; }
        } elsif ($observed_indices == 3) {
            $triples{sum}++;
            foreach my $k (keys %observe) { $triples{$k} += $observe{$k}; }
        } elsif ($observed_indices == 4) {
            $quads{sum}++;
            foreach my $k (keys %observe) { $quads{$k} += $observe{$k}; }
        } elsif ($observed_indices == 5) {
            $fives{sum}++;
            foreach my $k (keys %observe) { $fives{$k} += $observe{$k}; }
        } elsif ($observed_indices == 6) {
            $sixes{sum}++;
            foreach my $k (keys %observe) { $sixes{$k} += $observe{$k}; }
        } elsif ($observed_indices == 7) {
            $sevens{sum}++;
            foreach my $k (keys %observe) { $sevens{$k} += $observe{$k}; }
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
            number => int($number),
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
## Here we create the summary file
sub End_Handler {
    print $log "${observed_reads} reads were observed in total, of these:\n";
    foreach my $k (sort keys %observations) {
        if ($observations{$k} > 0 and $k ne 'sum') {
            ## Here we print how many times each individual index class was observed (regardless of if any other indices were found within that sequence). An index class is either the large variable region, small (helical) variable region, or synthesis variable region.
            print $log "The read type: ${k} was observed: $observations{$k} times.\n";
        }
    }
    ## how many times single, double...six index sequences were found & how many times individual index classes were found within them.
    if ($singles{sum} > 0) {
        print $log "$singles{sum} single-index reads were observed, including:\n";
        foreach my $k (sort keys %singles) {
            if ($singles{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $singles{$k} times.\n";
            }
        }
    }
    if ($doubles{sum} > 0) {
        print $log "$doubles{sum} double-index reads were observed, including:\n";
        foreach my $k (sort keys %doubles) {
            if ($doubles{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $doubles{$k} times.\n";
            }
        }
    }
    if ($triples{sum} > 0) {
        print $log "$triples{sum} triple-index reads were observed, including:\n";
        foreach my $k (sort keys %triples) {
            if ($triples{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $triples{$k} times.\n";
            }
        }
    }
    if ($quads{sum} > 0) {
        print $log "$quads{sum} 4-index reads were observed, including:\n";
        foreach my $k (sort keys %quads) {
            if ($quads{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $quads{$k} times.\n";
            }
        }
    }
    if ($fives{sum} > 0) {
        print $log "$fives{sum} 5-index reads were observed, including:\n";
        foreach my $k (sort keys %fives) {
            if ($fives{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $fives{$k} times.\n";
            }
        }
    }
    if ($sixes{sum} > 0) {
        print $log "$sixes{sum} 6-index reads were observed, including:\n";
        foreach my $k (sort keys %sixes) {
            if ($sixes{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $sixes{$k} times.\n";
            }
        }
    }
    if ($sevens{sum} > 0) {
        print $log "$sevens{sum} 7-index reads were observed, including:\n";
        foreach my $k (sort keys %sevens) {
            if ($sevens{$k} > 0 and $k ne 'sum') {
                print $log "The read type: ${k} was observed: $sevens{$k} times.\n";
            }
        }
    }
    print $log "\n";
    print $log "${found_all_four} reads had a helical+stepcyc+variable and cyclized:\n";
    foreach my $k (sort keys %cyclized_final_lengths) {
        print $log "Size $k was found $cyclized_final_lengths{$k} times and cyclized.\n";
        if ($k ne '72') {
            print $cyc_csv "$k,$cyclized_final_lengths{$k}\n";
        }
    }
    print $log "${found_all_four} reads had a helical+stepcyc+variable and not cyclized:\n";
    foreach my $k (sort keys %notcyclized_final_lengths) {
        print $log "Size $k was found $notcyclized_final_lengths{$k} times and not cyclized.\n";
        if ($k ne '72') {
            print $not_csv "$k,$notcyclized_final_lengths{$k}\n";
        }
    }
    $log->close();
    $cyc_csv->close();
    $not_csv->close();
    $out->close();
    exit(0);
}
