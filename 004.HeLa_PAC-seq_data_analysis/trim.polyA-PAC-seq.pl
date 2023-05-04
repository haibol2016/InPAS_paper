#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $trim_polyG = 0;
my $trim_polyA = 1;
my $trim_polyT = 0;

my ($windowSize, $stepSize, $min_content, $min_length) = (15, 1, 90, 25);
sub get_base_content;
sub get_kmer_base_content;
sub get_trimming_site;
sub trim_read;


while (my $header =  <STDIN>)
{
   chomp $header;
   
   my @items = split /\s+/, $header;
   $header = $items[0];
   
   chomp (my $seq = <STDIN>);
   chomp (my $conj = <STDIN>);
   chomp (my $qual = <STDIN>);
   
   my $polyA_len = 0;
   
   # ### trimming non-polyA and  polyA at the 3' end only if the remaining length is >=25
   if (length($seq) >= $min_length)
   {
     if ($seq =~/A{10,}/ && $trim_polyA)
     {
           ## reverse sequence and quality to remove nucleotides after polyA (A% in 15 bp > 90*)
           $seq = reverse $seq;
           $qual = reverse $qual;
           
          ($seq, $qual) = @{trim_read($seq, $qual, "A", "5")};
           
           ## reverse sequence and quality back to the read directions
           $seq = reverse $seq;
           $qual = reverse $qual;
           
          ($seq, $qual, $polyA_len) = @{trim_read($seq, $qual, "A", "3")};
     }
   }
   
   ### output reads if the current length is >= $min_length
    if (length($seq) >= $min_length)
    {
        $conj = "+";
        ## add polyA length to header if $polyA_len !=0;
        $header = $header."_".$polyA_len."A" if $polyA_len;
        my @fastq = ($header, $seq, $conj, $qual);
        print join "\n", @fastq;
        print "\n";
    }
}

sub get_base_content
{
    my ($kmer, $bases) = @_;
    
    no strict;
    no warnings;
    my @nt_count = eval "\$kmer =~ tr/$bases//";
    my @nt_pct = map {$_ / length($kmer)} @nt_count;
    my %freq = ("$bases" => $nt_pct[0]);
    return \%freq;
}


sub get_kmer_base_contents
{
    my ($sequence, $windowSize, $stepSize, $base) = @_;
    my @kmers=();
    for( my $windowStart = 0 ; $windowStart <= (length($sequence) - $windowSize);
        $windowStart += $stepSize)
    {
      my $kmer = substr($sequence, $windowStart, $windowSize);
      
      my $kemr_nt_content = get_base_content($kmer, $base); ## reference to a hash
      push @kmers, {$kmer => $kemr_nt_content};  ## array of anonymous hash of anonymous hash
    }
    return \@kmers;
}


sub get_trimming_site
{
    ## $min_content: min G contents,considered as polyG stretch
    ## $base to trim
    my ($kmer_base_content_ref, $base, $min_content, $read_length, $stepSize) =  @_;
    my $highest_content = $min_content;
    my @base_pct = ();
    my @kmers = ();
    
    my @kmer_base_content_tmp = @$kmer_base_content_ref; ## dereference the first level to get
                                                         ## array of anonymous hash of anonymous hash
    ## add kmer sequence and content of base to be trimmed to two parallel arrays
    for my $kmer (@kmer_base_content_tmp)
    {
       push @kmers, keys(%{$kmer});
       
       # print Dumper(%{(values (%{$kmer}))[0]});  ### values(%{$kmer}) return an array of a
                                                   ### single element, which is an anonymous hash
       push @base_pct, ${(values (%{$kmer}))[0]}{$base};
    }

    
    my $trimming_start = $read_length;
        
    for (my $i=0; $i < scalar(@base_pct); $i++)
    {
        if ($kmers[$i] =~/^$base/ && $base_pct[$i]*100 >= $min_content)
        {        
            ## scan the following kemrs to make sure the minumal baseq% is satisfied
            $trimming_start = $i * $stepSize;   # first high A site
            return $trimming_start;
        }
    }
    return $trimming_start;
}

sub trim_read
{
    my ($seq, $qual, $trim_base, $end) = @_;
    my $read_length = length($seq);
    my $kmer_base_content_ref = get_kmer_base_contents($seq, $windowSize, $stepSize, $trim_base);
    my $triming_start = get_trimming_site($kmer_base_content_ref, $trim_base,
                                          $min_content, $read_length, $stepSize);
    if ($end eq "5")
    {
         $seq  = substr($seq, $triming_start);
         $qual = substr($qual, $triming_start);
          ## return anomymous array containing $seq and quality after trimming
         return [$seq, $qual];
    } elsif ($end eq "3")
    {
         $seq  = substr($seq, 0, $triming_start);
         $qual = substr($qual, 0, $triming_start);
    
         # number of polyA with mismatch (>90% A in every 15 bp)
          my $polyA_len = $read_length - length($seq);
    
         ## return anomymous array containing $seq and quality after trimming
         return [$seq, $qual, $polyA_len];
    }
}
