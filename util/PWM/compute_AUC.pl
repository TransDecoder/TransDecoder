#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw(min);

my $usage = "\n\n\tusage: $0 feature.scores.roc\n\n";

my $roc_file = $ARGV[0] or die $usage;



main: {

    my $auc_outfile = "$roc_file.auc";
    open (my $ofh, ">$auc_outfile") or die "Error, cannot write to $auc_outfile";
        
    my %data = &parse_roc($roc_file);
    
    foreach my $cat (sort keys %data) {
        
        my @vals = sort {$a->[0] <=> $b->[0]
                             ||
                             $a->[1] <=> $b->[1] } @{$data{$cat}};
        
        unshift (@vals, [0,0]);
        push (@vals, [1,1]);

        my $total_auc = 0;

        my $prev_vals_aref = shift @vals;
        while (@vals) {
            my $next_vals_aref = shift @vals;

            my ($x1, $y1) = @$prev_vals_aref;
            my ($x2, $y2) = @$next_vals_aref;

            ########
            #
            #     * (x2,y2)
            #    /|
            #   / |
     # (x1,y1) *  |  
            #  |  |
            #  |  |
            #  |  |        area = rectangle + triangle
            #----------         = (x2-x1) * y1
            #                               + 1/2 * (x2-x1) * abs(y2-y1)

            my $area = ( ($x2-$x1) * min($y1,$y2) ) # rect
                + ( 0.5 * ($x2-$x1) * abs($y2-$y1) ); # triangle

            $total_auc += $area;
            
            $prev_vals_aref = $next_vals_aref;
        }
        
        print $ofh join("\t", $cat, $total_auc) . "\n";
    }

    close $ofh;


    my $Rscript_file = "$auc_outfile.plot.Rscript";
    {
        open (my $ofh, ">$Rscript_file");
        print $ofh "data = read.table(\"$auc_outfile\", header=F)\n"
            . "colnames(data) = c('cat', 'auc')\n"
            . "pdf(\"$auc_outfile.plot.pdf\")\n"
            . "barplot(data[,2], las=2, names=data[,1], cex.names=0.4, ylim=c(0,1))\n"
            . "data = data[rev(order(data[,2])),]\n"
            . "barplot(data[,2], las=2, names=data[,1], cex.names=0.4, ylim=c(0,1))\n"
            . "library(ggplot2)\n"
            . "before_after_df = data.frame(t(simplify2array(strsplit(as.character(data\$cat), ','))))\n"
            . "before_after_df = apply(before_after_df, 1:2, as.numeric)\n"
            . "colnames(before_after_df) = c('before', 'after')\n"
            . "data = cbind(before_after_df, data)\n"
            . "ggplot(data, aes(x=before, y=after)) +  geom_point(aes(size=auc, color=auc))\n"
            . "dev.off()\n";
        close $ofh;
        
        system("Rscript $Rscript_file");

    }

    
    exit(0);
}

####
sub parse_roc {
    my ($roc_file) = @_;

    my %data;
    
    open (my $fh, $roc_file) or die "Error, cannot open file $roc_file";
    my $header = <$fh>;
    while (<$fh>) {
        chomp;

        my @x = split(/\t/);
        my $cat = $x[0];
        my $TPR = $x[6];
        my $FPR = $x[7];

        push (@{$data{$cat}}, [$FPR, $TPR]);

    }
    close $fh;

    return(%data);
}


