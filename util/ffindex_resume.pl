#!/usr/bin/env perl

use strict;
use warnings;

my $input_idx = shift;
my $done_input_idx = shift||die ("input IDX, basename of done IDX\n");
my $done_input = $done_input_idx;
$done_input =~s/\.idx$/.db/;
my $counter=int(0);
my (%hash_done,%hash_all);

open (IN,"$input_idx") ||die;
while (my $ln=<IN>){
        $ln=~/^(\d+)\s/ ||next;
        $hash_all{$1} = $ln;
}
close IN;


my @files_done = glob("./$done_input_idx*");
foreach my $done (@files_done){
 open (IN,$done);
  while (my $ln=<IN>){
    $ln=~/^(\d+)\s/ ||next;
    next if $ln=~/\b1$/;
    $hash_done{$1}=1;
  }
 close IN;
}

print "All: ".scalar(keys(%hash_all))."\n";
print "Done: ".scalar(keys(%hash_done))."\n";

open (OUT,">".$input_idx.".notdone");

foreach my $i (sort {$a <=> $b} keys %hash_all){
 next if $hash_done{$i};
 $counter++;
 print OUT $hash_all{$i};
}

close OUT;

print "Found $counter unfinished\n";
if ($counter==0){
   unlink($input_idx.".notdone");
   system("cat $done_input | tr -d '\000' > $done_input.txt");
   exit;
}

