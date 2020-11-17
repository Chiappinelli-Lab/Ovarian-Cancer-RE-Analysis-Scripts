#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;
use Statistics::R;

# This script will take two, tab-delimited files as input.
# They should have this format: "<RE ID>\t<Count>\n".
# The script will calculate the enrichment of the REs in the group RE ID from 
# the sample relative to the reference file.
# In this case, the reference file will usually be the RepeatMasker counts.
# The statistical test will be the binomial test, performed with R.
# P-values will be corrected for multiple testing with Bonferroni correction.
# It will output the RE ID, the counts used, P-value, and Q-value (corrected P-value).

# Set variable defaults and usage declaration.
my $pvalfilter = 1;
my $fileextension = "";
my $usage = "USAGE: countRE.pl <Options: -f -n> <Ref. File> <Sample File>
	-f\t-filter\tMaximum q-value allowed for output. Optional. Default = 1 (no filter).
	-n\t-name\tString for suffix to output file name. Optional. Default = empty string.
	NOTE: The Bonferroni correction is applied to the p-value, creating a q-value. Filter is based on q-values.
	NOTE: The script uses the Statistics::R library. Please sync the Perl and library version.\n";

# Get user input
GetOptions("filter=f" => \$pvalfilter, "name=s" => \$fileextension) or die($usage);

# Check that the proper number of files were passed
if (@ARGV != 2)
{
        die("Incorrect number of files specified\n$usage");
}

# Get and/or make needed file names
my $reffilename = $ARGV[0];
my $samplefilename = $ARGV[1];
my $outfilename = join(".", substr($samplefilename, 0, -11), "enrich", "binom", $fileextension);

# Open necessary files
open (REF, $reffilename) or die("Couldn't open $reffilename, $!\n\n");
open (SAMPLE, $samplefilename) or die("Couldn't open $samplefilename, $!\n\n");
open (OUT, '>', $outfilename) or die("Couldn't open $outfilename, $!\n\n");

# Start a communication bridge between R and the script
my $R = Statistics::R->new();

# Global variables
my $refREtotal = 0;
my $sampleREtotal = 0;
my %refREcounts;
my %sampleREcounts;

# Make the reference hash and calculate total RE in reference
while(my $line = <REF>)
{
	$line =~ tr/\r/\n/; chomp($line); chomp($line); # Turn \r into newlines (\n) and remove newlines
	my @lineelements = split(/\t/, $line); 		# If correct format element 0 = RE ID and element 1 = RE count
	$refREcounts{$lineelements[0]} = $lineelements[1]; # This is m, see below
	$refREtotal += $lineelements[1]; 		   # This is N, see below
}
my @keys = keys(%refREcounts);
my $hsize = @keys;
print "Ref Hash Size: $hsize\n";

# Make the sample hash and calculate total RE in sample
while(my $line = <SAMPLE>)
{
        $line =~ tr/\r/\n/; chomp($line); chomp($line);	# Turn \r into newlines (\n) and remove newlines
	my @lineelements = split(/\t/, $line); 		# If correct format element 0 = RE ID and element 1 = RE count
        $sampleREcounts{$lineelements[0]} = $lineelements[1]; # This is s, see below
        $sampleREtotal += $lineelements[1]; 		      # This is n, see below
}
@keys = keys(%sampleREcounts);
$hsize = @keys;
print "Hash size: $hsize\n";

# Troubleshooting varaible for counting the loops below.
my $loopcount = 0;

# For each RE ID (hash key) in sample, calculate enrichament by hypergeometric test

# Notes on variables used:
# N = Number of RE in RepeatMasker
# m = Number of RE from a specific class/name (i.e. total number of successes)
# n = Number of RE in DMRs from a sample (i.e. the number of draws/trials)
# s = Number of RE in DMRs from a sample that are of the specific class/name (successes drawn)

# Notes on variables for R:
# q    = s     = number of successes drawn
# size = n     = number of draws, i.e. the sample size
# prob = m / N = number of non-successes in the total population
foreach my $reID (keys %sampleREcounts)
{
	# Set variables
	my $refcount = $refREcounts{$reID};		# This is m
	my $samplecount = $sampleREcounts{$reID};	# This is s
	my $probability = $refcount / $refREtotal;	# This is m / N, the probability of success
	$samplecount--;					# Decrement b/c lower.tail=F calcs P(X > q) not P(X >= q)
	
	# Put variables into R instance and run test	
	$R->set('success', $samplecount);
	$R->set('size', $sampleREtotal);
	$R->set('prob', $probability);
	$R->run(q'binompval <- pbinom(success, size, prob, lower.tail = FALSE)');
	my $pvalue = $R->get('binompval');
	my $qvalue = $pvalue * $hsize;			# This is the Bonferroni correction, 1 test for each member of the hash, so use the hash size calculated above
	if ($qvalue > 1) { $qvalue = 1; }		# Largest possible q-value is 1
	$samplecount++;					# Increment to get original value
	
	if($qvalue <= $pvalfilter)
	{
		print OUT "$reID\t$refREtotal\t$refcount\t$sampleREtotal\t$samplecount\t$pvalue\t$qvalue\n";
	}
	$loopcount++;
}
print "Foreach loop count: $loopcount\n";
