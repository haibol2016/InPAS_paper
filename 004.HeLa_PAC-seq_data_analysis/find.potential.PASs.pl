#!/bin/env perl

############# USAGE ######################
# samtools view  XXX.srt.bam |  perl find.potential.PASs.pl  - PAC-seq  > XXX.candidate.PASs.bedgraph
#
##########################################

use strict;
use warnings;
use aliased 'Hash::DefaultValue' => 'HDV';
use Data::Dumper;

die "Please specify the 3' sequencing protocol: PAC-seq or polyA-seq" if scalar @ARGV != 2;

my %strands;

if ($ARGV[1] eq "PAC-seq")
{
    %strands = (
                0 => "+",
                16 => "-",
                272 => "-",
                256 => '+');
}
elsif ($ARGV[1] eq "polyA-seq")  # reverse complementary strand sequenced
{
    %strands = (
                0 => "-",
                16 => "+",
                272 => "+",
                256 => '-');    
}


my @cigar_key = qw(M D = X N);

my $start = 0;
my $start_chr = '';
# my $curr_chr = '';
my @chr_records;
sub read_alignment;
sub process_chr;

while (<STDIN>)
{
    chomp;
    my @items = split /\t/;
    
    next if $ARGV[1] eq "PAC-seq" && $items[0] !~ /A$/;
    
    if (!$start)
    {
        $start = 1;
        $start_chr = $items[2];
        # $curr_chr = $items[2];
    }
    elsif ($items[2] eq $start_chr)
    {
        my $line = read_alignment(@items);  # anonymous array
        push @chr_records,  $line;
        # print Dumper(@chr_records);
        # $curr_chr = ${$line}[0];
    }
    else  ## new chromosome start, processing @chr_records
    {
        process_chr(\@chr_records);
    
        # re-initiate @chr_records;
        @chr_records = ();
        my $line = read_alignment(@items);
        push @chr_records,  $line;
        $start_chr = $items[2];
        # $curr_chr = ${$line}[0];
    }
}
process_chr(\@chr_records);
 
 
sub read_alignment
{
    my @items = @_;
    
    my $A_num = $ARGV[1] eq "PAC-seq" ? ($items[0] =~ /_(\d+)A/)[0] : 0;        
    my $strand = $strands{$items[1]};
      
    my $curr_chr = $items[2];
    
    my @nt_num = ($items[5] =~/(\d+)/g);
    my @nt_type = ($items[5] =~/([MIDNSHP=X])/g);

    ## make the clipping keys unique  
    if ($nt_type[0] eq "S" || $nt_type[0] eq "H")
    {
        $nt_type[0] = $nt_type[0]."_left";
    }
    if ($nt_type[$#nt_type] eq "S" || $nt_type[0] eq "H")
    {
        $nt_type[$#nt_type] = $nt_type[$#nt_type]."_right";
    }
    
    #print Dumper(@nt_num);
    #print Dumper(@nt_type);
    
    tie my %cigar, HDV, 0;
    my $three_end;
    
    for (my $i = 0; $i < scalar(@nt_type); $i++)
    {
        $cigar{$nt_type[$i]} +=  $nt_num[$i]; 
    }
    
    
    # print $strand." ".$items[3]." ".$items[5]."\n";
    #print Dumper(%cigar);
    
    if ($strand eq "-")
    {
        $three_end = $items[3];
    }
    else  ## + strand, add the distance the read spans: M, D, N, =, X
    {
        $three_end = $items[3];
        for my $k (@cigar_key)
        {
            $three_end += $cigar{$k}
        }  
    }
    [$curr_chr, $three_end, $strand, $A_num];
    
}

sub process_chr
{
    my @chr_record = @{$_[0]};
    
    tie my %pas, HDV, 0;  # number of reads
    tie my %pas_maxA, HDV, 0;  # maximal number of A's
    for my $rec (@chr_record)
    {
        my @elements = @{$rec};
        my $curr_pos = join ":", @elements[0..2];
        $pas{$curr_pos} += 1;
        
        my $curr_As = $elements[3];
        if ($curr_As > $pas_maxA{$curr_pos})
        {
            $pas_maxA{$curr_pos} = $curr_As;
        }
    }
    for my $pos (sort keys(%pas))
    {
        my @locus = split /:/, $pos;
        
        # CHR, POS, STRAND, # reads, # A's
        push @locus, $pas{$pos}, $pas_maxA{$pos};
        print join "\t", @locus;
        print "\n";
    }
}