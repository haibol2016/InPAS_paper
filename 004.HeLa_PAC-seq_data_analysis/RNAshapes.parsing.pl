#!/bin/env perl

use warnings;
use strict;

my $max_struct = 3;
my $start = 0;
my $curr_seq;
my $num_struct = 1;
my @rna_structs;

while (<>)
{
    chomp;
    if (!$start && /^>/)
    {
        $start = 1;
        $curr_seq = $_;
        push @rna_structs, $curr_seq;
    }
    elsif (/[ATCG]+/ || (!/^>/ && $num_struct > $max_struct))
    {
        next;
    }
    elsif (!/^>/ && $num_struct <= $max_struct)
    {
        $num_struct++;
        my @items = split /\s+/;
        push @rna_structs, $items[1];
    }
    elsif (/^>/)
    {
        print join "\n", @rna_structs;
        print "\n";
        @rna_structs = ();
        $num_struct = 1;
        $curr_seq = $_;
        push @rna_structs, $curr_seq;
    }
}
print join "\n", @rna_structs;
print "\n";
