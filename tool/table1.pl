# further process Table 1 from 19-12 ObamaCare Health Insurance Study
# input is tab-delimited with row labels in column 1, space-separated data in table 2
# emit table with same row labels, but data split into separate cells 
use strict;
use Carp;

my $IN='table1.txt';
my $OUT='table1.split.txt';

open(IN,$IN) || confess "Unable to open input file $IN: $!";
open(OUT,"> $OUT") || confess "Unable to open output file $OUT: $!";

while (<IN>) {
  chomp;
  my($label,$data)=split("\t");
  $data=~s/\"|,//g;
  my @data=split(' ',$data);
  print OUT join("\t",$label,@data),"\n";
}
close IN;
close OUT;
