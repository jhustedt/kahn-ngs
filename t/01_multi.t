# -*-Perl-*-
use Test::More qw"no_plan";
use File::Path qw"remove_tree";
use File::Copy qw"cp";
use String::Diff qw"diff";
use diagnostics;

system("mkdir -p t/multi");
my $run = system("./script/sizecount.pl --input t/multihit.fastq --outdir t/multi --ind share/index.txt");
my $res = open(STUFF, "zgrep helical t/multi/out.fastq.gz |");
my $result;
while (<STUFF>) {
    $result .= $_;
}
my $expected = qq"+135:helical:8:fwd 153:helical:8:fwd 73:stepcyc:77:fwd 162:stepcyc:77:fwd 162:stepcyc:77:fwd 401:stepsynth:77:fwd 368:variable:20:fwd 1  hits: 7 \n";

print "TESTME1:
$result
";

print "TESTME2:
$expected
";

ok($result eq $expected);
