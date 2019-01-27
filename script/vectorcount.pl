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


=cut

#options from Getopt::Long; defaults
my %options = (
    debug => 0,
    indices => 'index.txt',
    input => 'test.fastq.gz',
    outdir => 'output',
    outfastq => 'out',
    summary => 'summary.txt',
    spacer => 72,
    substitution => 0,
    insertion => 0,
    deletion => 0,
);
## This is a hash counting off how many times _every_ index is observed.
my $observed_reads = 0;
my %observations = (
    ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %singles = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %doubles = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %triples = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %quads = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %fives = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %sixes = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %sevens = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %eights = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %nines = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %tens = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);
my %elevenup = (
    sum => 0, ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
    ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,);

my $opt_result = GetOptions(
    "debug:i" => \$options{debug},
    "spacer:i" => \$options{spacer},
    "indices:s" => \$options{indices},
    "input:s" => \$options{input},
    "outdir:s" => \$options{outdir},
    "summary:s" => \$options{summary},
    "outfastq:s" => \$options{outfastq},
    "substitution:i" => \$options{substitution},
    "insertion:i" => \$options{insertion},
    "deletion:i" => \$options{deletion},
);
my $log = new FileHandle(">$options{outdir}/$options{summary}");
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
            ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
            ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,
            unknown => 0);
        my %positions = (
            ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
            ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,
            unknown => 0);
        my %numbers = (
            ampr_fwd => 0, forigin_fwd => 0, variable_fwd => 0, pucorigin_fwd => 0,
            ampr_rev => 0, forigin_rev => 0, variable_rev => 0, pucorigin_rev => 0,
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
                    if ($info->{name} eq 'ampr' && $info->{direction} eq 'fwd') {
                        $observations{ampr_fwd}++;
                        $observe{ampr_fwd}++;
                        $positions{ampr_fwd} = $st;
                        $observed_indices++;
                        $numbers{ampr_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'forigin' && $info->{direction} eq 'fwd') {
                        $observations{forigin_fwd}++;
                        $observe{forigin_fwd}++;
                        $positions{forigin_fwd} = $st;
                        $observed_indices++;
                        $numbers{forigin_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'pucorigin' && $info->{direction} eq 'fwd') {
                        $observations{pucorigin_fwd}++;
                        $observe{pucorigin_fwd}++;
                        $positions{pucorigin_fwd} = $st;
                        $observed_indices++;
                        $numbers{pucorigin_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'variable' && $info->{direction} eq 'fwd') {
                        $observations{variable_fwd}++;
                        $observe{variable_fwd}++;
                        $positions{variable_fwd} = $st;
                        $observed_indices++;
                        $numbers{variable_fwd} = $info->{number};
                    } elsif ($info->{name} eq 'ampr' && $info->{direction} eq 'rev') {
                        $observations{ampr_rev}++;
                        $observe{ampr_rev}++;
                        $positions{ampr_rev} = $st;
                        $observed_indices++;
                        $numbers{ampr_rev} = $info->{number};
                    } elsif ($info->{name} eq 'forigin' && $info->{direction} eq 'rev') {
                        $observations{forigin_rev}++;
                        $observe{forigin_rev}++;
                        $positions{forigin_rev} = $st;
                        $observed_indices++;
                        $numbers{forigin_rev} = $info->{number};
                    } elsif ($info->{name} eq 'pucorigin' && $info->{direction} eq 'rev') {
                        $observations{pucorigin_rev}++;
                        $observe{pucorigin_rev}++;
                        $positions{pucorigin_rev} = $st;
                        $observed_indices++;
                        $numbers{pucorigin_rev} = $info->{number};
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
        my $ampr = 0;
        my $forigin = 0;
        my $pucorigin = 0;
        my $variable = 0;
        my $bimol_valid = 0;
        my $amprfwd = 0;
        my $amprrev = 0;
        my $variablefwd = 0;
        my $variablerev = 0;
        my $foriginfwd = 0;
        my $foriginrev = 0;
        my $pucoriginfwd = 0;
        my $pucoriginrev = 0;
        ## append to the comment the number of observed indices.
        $comment .= "$count ";
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
## Here we create the overall summary file
sub End_Handler {
    print $log "Index used: $options{indices}\n";
    print $log "Input file: $options{input}\n";
    if ( ! defined $options{insertion} && ! defined $options{deletion} &&
            ! defined $options{substitution} ) {
        print $log "Direct match used";
    } else {
        print $log "String::Approx matching used: I$options{insertion},D$options{deletion},S$options{substitution}\n";
    }
    print $log "${observed_reads} reads were observed in total, of these:\n";
    foreach my $k (sort keys %observations) {
        if ($observations{$k} > 0 and $k ne 'sum') {
            print $log "     The read type: ${k} was observed: $observations{$k} times.\n";
        }
    }
    ## how many times single, double...seven index sequences were found & how many times individual index classes were found within them.
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

    $log->close();
    $unicyc_csv->close();
    $fourhitlin_csv->close();
    $bicyc4_csv->close();
    $bicyc6_csv->close();
    $bicyc2_csv->close();
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
