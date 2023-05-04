#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $trim_polyG = 0;
my $trim_polyA = 1;
my $trim_polyT = 1;

my ($windowSize, $stepSize, $min_content, $min_length) = (10, 1, 80, 25);
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
   
   ## trimming polyG from the 3' end, for Nextseq data
   if ($seq =~/G{10,}/ && $trim_polyG && length($seq) >= $min_length)
   {
      ($seq, $qual) = @{trim_read($seq, $qual, "G")};
   }
    
   # ### trimming polyT at the 5' end only if the remaining length is >=25
   if (length($seq) >= $min_length)
   {
     if ($seq =~/T{10,}/ && $trim_polyT)
     {
           ## reverse sequence and quality
           $seq = reverse $seq;
           $qual = reverse $qual;
           
           ($seq, $qual) = @{trim_read($seq, $qual, "A")};
           
           ## reverse sequence and quality back to the read directions
           $seq = reverse $seq;
           $qual = reverse $qual; 
     }
   }

    ### trimming polyA at the 3' end only if the remaining length is >=25 
    ### even after trimming polyG
    if (length($seq) >= $min_length)
    {
        if ($seq =~/A{10,}/ && $trim_polyA)
        {         
            ($seq, $qual) = @{trim_read($seq, $qual, "A")};
        }
    }
    
    ### trimming long reads to 50 bases if they are longer than 50
    # if (length($seq) >= 50)
    # {
      # $seq  = substr($seq, 0, 50);
      # $qual = substr($qual, 0, 50);
    # }
    
    ### output reads if the current length is >= $min_length
    if (length($seq) >= $min_length)
    {
        $conj = "+";
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
            ## scan the following kemrs to make sure the minimal baseq% is satisfied
            my $all_high = 1;
            for (my $j = $i +1; $j < scalar(@base_pct); $j++)
            {
                if ($base_pct[$j]*100 < $min_content)
                {
                    $all_high = 0;
                    last;
                }
            }
            if ($all_high)
            {
                $trimming_start = $i * $stepSize; ## trimming site
                ## print $trimming_start."\n";
                return $trimming_start;
            }  
        }
    }
    return $trimming_start;
}

sub trim_read
{
    my ($seq, $qual, $trim_base) = @_;
    my $read_length = length($seq);
    my $kmer_base_content_ref = get_kmer_base_contents($seq, $windowSize, $stepSize, $trim_base);
    my $triming_start = get_trimming_site($kmer_base_content_ref, $trim_base,
                                          $min_content, $read_length, $stepSize);
    $seq  = substr($seq, 0, $triming_start);
    $qual = substr($qual, 0, $triming_start);
    
    ## return anomymous array containing $seq and quality after trimming
    return [$seq, $qual];
}
