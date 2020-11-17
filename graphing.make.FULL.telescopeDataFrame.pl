#!/usr/bin/perl -w

use strict;
use warnings;
use Getopt::Long;

# This script will take a Telescope annotation file as input.
# There are two options: GTF and TSV
# Example GTF -- use if() to match only gene lines, then get the gene_id (array 8), chr (array 0), start (array 3), and stop (array 4)
# 	[jimcdonald@login3 results]$ head HERV_rmsk.hg38.v2.gtf
# 	### ERV316A3_1p36.33 ###
#	chr1    rmsk    gene    40736   42273   .       -       .       gene_id "ERV316A3_1p36.33"; transcript_id "ERV316A3_1p36.33"; locus "ERV316A3_1p36.33"; category "oneside"; intModel "ERV3-16A3_I-int"; locid "ERV3-16A3_I-int_0001"; model_cov "893"; model_pct "17.1";
#	chr1    rmsk    exon    40736   40878   260     -       .       gene_id "ERV316A3_1p36.33"; transcript_id "ERV316A3_1p36.33"; locus "ERV316A3_1p36.33"; repName "LTR16C"; geneRegion "ltr"; locid "ERV3-16A3_I-int_0001"; repClass "LTR"; repEnd "376"; repFamily "ERVL"; repLeft "232"; repStart "-113";
#	chr1    rmsk    exon    41394   42273   765     -       .       gene_id "ERV316A3_1p36.33"; transcript_id "ERV316A3_1p36.33"; locus "ERV316A3_1p36.33"; repName "ERV3-16A3_I-int"; geneRegion "internal"; locid "ERV3-16A3_I-int_0001"; repClass "LTR"; repEnd "4789"; repFamily "ERVL"; repLeft "3896"; repStart "-437";
# Example TSV -- grab the locus ID (array 0), chr (array 1), start (array 2), and end (array 3)
#	[jimcdonald@login3 results]$ head HERV_rmsk.hg38.v2.tsv
#	locus   chrom   start   end     strand  score   category        gene_id intModel        locid   model_cov       model_pct    transcript_id
#	ERV316A3_1p36.33        chr1    40736   42273   -       .       oneside ERV316A3_1p36.33        ERV3-16A3_I-int ERV3-16A3_I-int_0001 893     17.1    ERV316A3_1p36.33

# Currently all this information is in the ATAC-seq data file in the directory there.
# However, doing it this way lets the user specify a subset of the annotations.
# Also, the cumbersome ATAC-seq/LTR intersect output files could be trimmed down if needed.

# The script also requires several data directories to be specified:
# NOTE: ALL FILES MUST HAVE THE CELL LINE AND TREATMENT IN THEIR FILE NAME!!!
# 1. Directory of telescope expression data files in this format:
#	<locus name> <baseMean> <log2FoldChange> <lfcSE> <stat> <pvalue> <adj pvalue>
#	
#	Example:
#	[jimcdonald@login3 software]$ head /groups/chiappinellilab/processed_data/RNA-seq/all.lines.telescope.library001.hg38.results/results/Hey.treat4.AZA.ITF.telescope.DESeq2.tsv
#	"ERV316A3_1p36.33"      0       NA      NA      NA      NA      NA
#	"HML3_1p36.33"  0.69983861367698        -0.0456938707476243     4.4405930628    -0.0102900378623778     0.99178988254767     0.995828354145708 
#
# 2. Direcotry of methylation data from methylCRF in BED format:
#	<chr> <start> <end> <CpG ID> <5mC level>		
#	
#	Example:
#	[jimcdonald@login3 software]$ head /lustre/groups/chiappinellilab/processed_data/methylCRF/mCRF.A549.library001.treat2.ITF/mCRF.A549.library001.treat2.ITF.sorted_methylCRF.bed
#	chr1    10468   10470   chr1.1  0.95
#	chr1    10470   10472   chr1.2  0.98
#
# 3. A direcotry of age files from Matthew's ERV age analysis for telescope:
#	See the GitHub: https://github.com/mlbendall/HERV_age_estimation
#	Format:
#	<locus> <category> <subfamily> <p.dist> <jc.dist> <tn.dist> <k2p.dist> <t3p.dist> <cof.dist>
#
#	Example:
#	[jimcdonald@login3 software]$ head /lustre/groups/cbi/Projects/HERV_age_estimation.git/HERV4/distances.tsv
#	locus   category        subfamily       p.dist  jc.dist tn.dist k2p.dist        t3p.dist        cof.dist
#	HERV4_1p35.3    internal        NA      inf     inf     inf     inf     inf     inf
#	HERV4_1p31.1    internal        NA      inf     inf     inf     inf     inf     inf
#	HERV4_1q22a     oneside_3       MER51A  0.1717791       0.195088        0.2010842       0.2007222       0.2008134   0.189189189189
#
# 4. A directory of files that contain the overlap of the Telescope annotation data with the ATAC-seq peaks unique to each sample.
#	Format:
#	<chr> <LTR start> <LTR stop> <LTR strand> <score> <category> <gene ID (USE THIS LOCUS NAME)> <intModel (internal ERV model)?> <unique intModel locus ID?> <model_cov (intModel coverage?)> <model_pct> <transcript ID> <chr> <Peak start> <Peak end> <Peak ID> <unused: "."> <unused: "."> <open state start> <open state end> <color code> <num. sub-regions> <list sub-regions> <starts ea. sub-region> <Peak Score> <unused: -1> <unused: -1>
#	I want to use this to join together the gene ID (array element 6) with the ATAC-seq peak score (array element 24).
#	The overlap between the peaks is in array element 27.
#
#	Example:
#	[jimcdonald@login3 software]$ head /lustre/groups/chiappinellilab/processed_data/paper.graphs/telescope.EXPR.5mC.ATACSeq.correlation/intersect.teleLTRs.uniquePeaks.Kuramochi.ATACseq.library001.treat3.AZA.bed
#	chr1    40736   42273   -       .       oneside ERV316A3_1p36.33        ERV3-16A3_I-int ERV3-16A3_I-int_0001    893 17.1     ERV316A3_1p36.33        .       -1      -1      .       .       .       .       .       .       .       .   ..       .       .       0
#	chr1    895544  898699  -       .       oneside HML3_1p36.33    HERVK9-int      HERVK9-int_0002 2681    44.5    HML3_1p36.33 chr1    897725  899055  Peak_285196     .       .       898105  898615  255,0,0 3       1,510,1 0,380,1329  20.0     -1      -1      974
#	chr1    1275713 1280788 -       .       prototype       MER4B_1p36.33   MER4B-int       MER4B-int_0001  1740    25.7MER4B_1p36.33    chr1    1280760 1282950 Peak_285391     .       .       1280950 1282520 255,0,0 3       1,1570,1    0,190,2189       141.0   -1      -1      28
#	chr1    1412252 1418852 -       .       prototype       HARLEQUIN_1p36.33       Harlequin-int   Harlequin-int_0001  5443     78.9    HARLEQUIN_1p36.33       chr1    1417800 1420340 Peak_285449     .       .       1417900 1420030 255,0,0      3       1,2130,1        0,100,2539      37.0    -1      -1      1052	

# The script will gather the appropriate data from each of the files in each of the directories.
# It will then go through the list of LTR loci from Telescope and print out the joined information
# from all files for that:

# Output format:
# <LTR Locus Name> <Cell Line> <Treatment> <LTR Age> <Expr log2FC relative to Mock> <ATAC-seq Peak Score> <Avg Mock 5mC> <Avg Treat 5mC> <Avg Treat - Avg Mock 5mC> <Min 5mC> <Max 5mC>
# This should make a nice tibble in R for analysis.

# Set variable defaults and usage declaration.
my ($exprDir, $methylDir, $distDir, $atacDir) = ("", "", "", "");
my $annotationType = "TSV"; 
my $outputFile = "R.ready.telescope.ALL.data.frame.txt";
my ($locusElement, $chrElement, $startElement, $endElement) = -1;
my $usage = "USAGE: graphing.make.EXPRdataFrame.pl <Options: -e, -m, -d, -a, -t, -n> <Annotation File>

	Options:
	-e\t-exprDir\tDirectory of expression files. Required.
	-m\t-methylDir\tDirectory of methylation data. Required.
	-d\t-distDir\tDirectory for distance files. Required.
	-a\t-atacDir\tDirectory for ATACseq files. Required.
	-t\t-typeAnnotation\tType of annotation file. Options: GTF or TSV. TSV is default. Required.
	-n\t-name\tFile name for output. Optional.
		Default is R.ready.telescope.data.txt in working directory.

	Output will contain data from accross all 4 input data types:
	<LTR Locus Name> <Cell Line> <Treatment> <LTR Age> <Expr log2FC relative to Mock> 
	<ATAC-seq Peak Score> <Avg Mock 5mC> <Avg Treat 5mC> <Avg Treat - Avg Mock 5mC> <Min 5mC> <Max 5mC>
	
	Purpose of output: input into R as a data frame/tibble for ggplot2 graphs so that
	the data can be analyzed as a tibble.";

# Get user input
GetOptions("exprDir=s" => \$exprDir,
	   "methylDir=s" => \$methylDir,
	   "distDir=s" => \$distDir,
	   "atacDir=s" => \$atacDir,
	   "typeAnnotation=s" => \$annotationType,
	   "name=s" => \$outputFile) or die("GetOptions failed.\n$usage\n");

# Check that all needed directories were specified
if ($exprDir eq "" | $methylDir eq "" | $distDir eq "" | $atacDir eq "") { 
	die("One or more data directories not specified.\n$usage"); }

# Set the correct array elements for the annotation file given by the user
if ($annotationType eq "GTF") 
{
	$locusElement = 8; $chrElement = 0; $startElement = 3; $endElement = 4; 
	print "Annotation type = GTF\n";
	die("WARNING: GTF code does not work. Please use a TSV file in the meantime. Sorry for the inconvenience.\n");
}
elsif ($annotationType eq "TSV")
{ $locusElement = 0; $chrElement = 1; $startElement = 2; $endElement = 3; print "Annotation type = TSV\n"; }
elsif ($annotationType eq "")
{ die ("Annotation type not specified.\n$usage"); }
else
{ die ("Annotation type incorrectly specified.\n$usage"); }

# Check that there was an annotation file
if (@ARGV != 1) { die ("Annotation file not specified, $usage\n\n"); }

# Declare the main data hashes and two hashes to hold the cell lines and treatments used
my (%expr_hash, %mC_hash, %dist_hash, %atac_hash, %atac_overlap, %annotations, %seen_lines, %seen_treats);

# Open input/output files
opendir (EXPRDIR, $exprDir) or die ("Couldn't open $exprDir, $!\n\n");
opendir (MCDIR, $methylDir) or die ("Couldn't open $methylDir, $!\n\n");
opendir (DISTDIR, $distDir) or die ("Couldn't open $distDir, $!\n\n");
opendir (ATACDIR, $atacDir) or die ("Couldn't open $atacDir, $!\n\n");
open (ANNOTATION, $ARGV[0]) or die ("Coultn't open $ARGV[0], $!\n\n");
open (OUT, '>', $outputFile) or die ("Couldn't open $outputFile, $!\n\n");

########
# MAIN #
########

# Make the expression hashes
while(my $inFile = readdir EXPRDIR)
{
	if ($inFile eq "." | $inFile eq "..") { next; } # Colonial One directories often contain these and I want to skip them
        $inFile = join("", $exprDir, $inFile);
	my $isdone = "";
	$isdone = fill_EXPR_hash($inFile);
	print STDERR "$isdone\n";
}

# Make the methylation hashes
while(my $inFile = readdir MCDIR)
{
	if ($inFile eq "." | $inFile eq "..") { next; } # Colonial One directories often contain these and I want to skip them
        $inFile = join("", $methylDir, $inFile);
	my $isdone = "";
        $isdone = fill_MC_hash($inFile);
        print STDERR "$isdone\n";
}

# Make the Distance hashes
while(my $inFile = readdir DISTDIR)
{
	if ($inFile eq "." | $inFile eq "..") { next; } # Colonial One directories often contain these and I want to skip them
	my $isdone = "";
        if ($inFile =~ /distances/) # A few other files will unavoidable be in the directory so skip those
	{
		$inFile = join("", $distDir, $inFile);
	        $isdone = fill_DIST_hash($inFile);
        }
	else { next; }
	print STDERR "$isdone\n";
}

# Make the ATACseq hashes
while(my $inFile = readdir ATACDIR)
{
	if ($inFile eq "." | $inFile eq "..") { next; } # Colonial One directories often contain these and I want to skip them
	$inFile = join("", $atacDir, $inFile);
	my $isdone = "";
        $isdone = fill_ATAC_hash($inFile);
        print STDERR "$isdone\n";
}

# Hash print for troubleshooting
#print "Seen lines and treatments for troubleshooting:\n";
#foreach my $troublekey (keys %seen_lines)
#{
#       print "$troublekey = $seen_lines{$troublekey}\n";
#}
#foreach my $troublekey (keys %seen_treats)
#{
#       print "$troublekey = $seen_treats{$troublekey}\n";
#}

# Process the data for each region in the annotation file:
while (my $line = <ANNOTATION>)
{
	if ($annotationType eq "TSV")
	{
		# Skip a header by searching for two column names in Matthew's headers
		if ($line =~ /strand/ && $line =~ /intModel/) { next; }

		# Remove rogue new lines and split so manipulation of the elements is easy.
		$line =~ tr/\r/\n/; chomp($line); chomp($line);
		my @lineelements = "";
		@lineelements = split(/\t/, $line);
	
		my $locus = $lineelements[$locusElement];
		my $chr = $lineelements[$chrElement];
		my $start = $lineelements[$startElement];
		my $end = $lineelements[$endElement];
		my $coordinates = join("-", $start, $end);
		$coordinates = join(":", $chr, $coordinates);
		$annotations{$locus} = $coordinates;
	}
	elsif ($annotationType eq "GTF") 
	{
		# Remove rogue new lines and split so manipulation of the elements is easy.
                $line =~ tr/\r/\n/; chomp($line); chomp($line);
                my @lineelements = "";
                @lineelements = split(/\t/, $line);

		# The GTF has headers and lines for exons (really sub-sections) of the LTRs I want to skip those
		if ($line =~ /###/) { next; }
		if ($lineelements[2] eq "gene")
		{
			my $locus = "";
	                if ($lineelements[$locusElement] =~ /^gene_id "([^"]+)";/) { $locus = $1; } # THIS IS THE TROUBLE LINE. THE GREEDY REGEX GETS THE ENTIRE ANNOTATION PART
        	        my $chr = $lineelements[$chrElement];
                	my $start = $lineelements[$startElement];
	                my $end = $lineelements[$endElement];
        	        my $coordinates = join("-", $start, $end);
	                $coordinates = join(":", $chr, $coordinates);
        	        $annotations{$locus} = $coordinates;
		}
	}	
}

print "\nStarting final output!\n\n";

foreach my $cell_line (keys %seen_lines )
{
	foreach my $treatment (keys %seen_treats )
	{
		if ($treatment eq "Mock") { print "Skipped a Mock treatment\n"; next; } # Skip mock cell lines. They only have 5mC data here and that is hard coded in the 5mC subroutine below
		foreach my $annotated_locus (keys %annotations)
		{
			my $output_expr = $expr_hash{$cell_line}{$treatment}{$annotated_locus};
			my $output_dist = $dist_hash{$annotated_locus};
			my ($treat_mC_avg, $mock_mC_avg, $mC_difference, $treat_mC_min, $treat_mC_max, $mock_mC_min, $mock_mC_max, $mC_count) = methylation_math($cell_line, $treatment, $annotations{$annotated_locus});
			my $output_atac_score = $atac_hash{$cell_line}{$treatment}{$annotated_locus};
			if ($output_expr eq "")
                        { print STDERR "Warning! No expression data for $cell_line $treatment $annotated_locus. Data set may not contain all combinations of cell line and treatment.\n"; }
			if ($output_atac_score eq "")
                        { print STDERR "Warning! No ATAC-seq peak score data for $cell_line $treatment $annotated_locus. Data set may not contain all combinations of cell line and treatment.\n"; }
			print OUT "$annotated_locus\t$cell_line\t$treatment\t$annotations{$annotated_locus}\t$output_dist\t$output_expr\t$output_atac_score\t$treat_mC_avg\t$mock_mC_avg\t$mC_difference\t$treat_mC_min\t$treat_mC_max\t$mock_mC_min\t$mock_mC_max\t$mC_count\n";
		}
	}
}



###############
# SUBROUTINES #
###############

sub fill_EXPR_hash
{
	# Open the file
	my $inFile = $_[0];
	open (INFILE, $inFile) or die("Couldn't open $inFile, $!\n\n");
	#print "File = $inFile\n";	

	# Reset the variables and get the new names for the cell line and treatment
        my ($cell_line, $treatment, $locus_name, $hashkey, $lastHashKeys, $lastEXPRval) = ("", "", "", "", "", "");
        ($cell_line, $treatment) = get_cellLine_treatment($inFile);
	#print "line = $cell_line\ntreat = $treatment\n";
	
	# Go through the file and put the expression information into the hash
	while(my $line = <INFILE>)
	{
		# Remove rogue new lines and split so manipulation of the elements is easy.
		$line =~ tr/\r/\n/; chomp($line); chomp($line);
		my @lineelements = "";		
		@lineelements = split(/\t/, $line);

		# Split the first field to get the locus name. This gets rid of the quotes and checks it's formatted correctly.
		# I can also get the family name from this with $2 if I want later, but not right now.
		if ($lineelements[0] =~ /^"(.*)_(.*)"/) { $locus_name = join("_", $1, $2); }
                else { print "$lineelements[0] is unusual. Skipping.\n"; next; }

		# Add it to the hash. Expression data in the 3rd column.
		# I need one global hash so I add cell line/treatment to keep keys unique
		#my $hashkey = join(".", $locus_name, $cell_line, $treatment); # Old line from before I knew about multi-dimensional hashes.
		$expr_hash{$cell_line}{$treatment}{$locus_name} = $lineelements[2];

		# Not sure why, but despige larger scope of $hashkey, it doesn't work in the return string. 
		# I have to do this to return the key value as if it had limited scope like @lineelements.
		# Whatever. It works, I'm not pausing to discover why. Gotta keep moving.
		$lastHashKeys = "$cell_line, $treatment, $locus_name";
		$lastEXPRval = $lineelements[2];
        }	

	#print "Hash values\n";
	#foreach my $derefkey (keys %{ $expr_hash{$cell_line}{$treatment} }) { print "Key = $derefkey. Value = $expr_hash{$cell_line}{$treatment}{$derefkey}\n"; }

	# Return an update for the user
	return("Finished processing $cell_line $treatment expression data.\nLast entry: Keys = $lastHashKeys. Value = $lastEXPRval.\n");
}

sub fill_MC_hash
{
        # Open the file
        my $inFile = $_[0];
        open (INFILE, $inFile) or die("Couldn't open $inFile, $!\n\n");
        #print "File = $inFile\n";

        # Reset the variables and get the new names for the cell line and treatment
        my ($cell_line, $treatment, $chr_name, $start, $last_chr, $last_start, $last_mCper) = ("", "", "", "", "", "", "");
        ($cell_line, $treatment) = get_cellLine_treatment($inFile);
        #print "line = $cell_line\ntreat = $treatment\n";

        # Go through the file and put the expression information into the hash
        while(my $line = <INFILE>)
        {
                # Remove rogue new lines and split so manipulation of the elements is easy.
		# Then get and store the data I need in a multidimensional hash
                $line =~ tr/\r/\n/; chomp($line); chomp($line);
                my @lineelements = "";
                @lineelements = split(/\t/, $line);
		$chr_name = $lineelements[0];
		$start = $lineelements[1];		
                $mC_hash{$cell_line}{$treatment}{$chr_name}{$start} = $lineelements[4];

		# Store the variables for the return report to the user
                $last_chr = $chr_name;
                $last_start = $start;
		$last_mCper = $lineelements[4];
        }

	#print "Hash values\n";
	#foreach my $derefkey (keys %{ $mC_hash{$cell_line}{$treatment}{$chr_name} }) { print "$mC_hash{$cell_line}{$treatment}{$chr_name}{$derefkey}\n"; }

	# Return an update for the user
        return("Finished processing $cell_line $treatment 5mC data.\nLast entry: chr = $last_chr, start = $last_mCper, %5mC = $last_mCper.\n");
}

sub methylation_math # Pass to this: 1) cell line, 2) treatment, 3) chromosome, 4) start, 5) end
{
	#print "Starting methylation_math subroutine\n";

	# Reset variables
	my ($cell_line, $treatment, $coordinates, $chr_name, $start, $end) = ("", "", "", "", "", "");
	my ($treat_mC_min, $mock_mC_min) = (2, 2);
	my ($treat_mC_max, $mock_mC_max) = (-1, -1); 
	my ($treat_mC_sum, $mock_mC_sum, $mC_count, $treat_mC_avg, $mock_mC_avg) = (0, 0, 0, 0, 0);
	my $mC_difference = -5;
	
	#print "Initialized mock min = $mock_mC_min. Initialized mock max = $mock_mC_max.\n";
	
	# Get passed variables
	$cell_line = $_[0]; $treatment = $_[1]; $coordinates = $_[2];
	if ($coordinates =~ /(.*):(.*)-(.*)/)
	{ $chr_name = $1; $start = $2; $end = $3; }

	# Use a for loop to go over all the coordinates in the region
	for (my $i = $start; $i < $end + 1; $i++)
	{
		if (exists $mC_hash{$cell_line}{$treatment}{$chr_name}{$i})
		{
			my $treat_is_Mock = "Mock";
			my $treat_mC_value = $mC_hash{$cell_line}{$treatment}{$chr_name}{$i};
			my $mock_mC_value = $mC_hash{$cell_line}{$treat_is_Mock}{$chr_name}{$i};
			if ($treat_mC_value > 1 || $treat_mC_value < 0 || $mock_mC_value > 1 || $mock_mC_value < 0)
			{ print STDERR "ERROR: $cell_line $treatment $chr_name $start mC value > 1 or < 0. Check this methylation value!\n\n"; }
			$treat_mC_sum += $treat_mC_value;
			$mock_mC_sum += $mock_mC_value;
			$mC_count++;
			if ($treat_mC_value > $treat_mC_max) { $treat_mC_max = $treat_mC_value; }
			if ($treat_mC_value < $treat_mC_min) { $treat_mC_min = $treat_mC_value; }
			if ($mock_mC_value > $mock_mC_max) { $mock_mC_max = $mock_mC_value; }
			if ($mock_mC_value < $mock_mC_min) { $mock_mC_min = $mock_mC_value; }
		}
	}

	# Calculate averages and differences
	if ($mC_count == 0) 
	{
		print STDERR "Warning! no CpGs located between $cell_line $treatment $coordinates. Possible explanation: data for this line expected, none given.\n";
		return ($treat_mC_avg, $mock_mC_avg, $mC_difference, $treat_mC_min, $treat_mC_max, $mock_mC_min, $mock_mC_max, $mC_count);
	}
	$treat_mC_avg = $treat_mC_sum / $mC_count;
	$mock_mC_avg = $mock_mC_sum / $mC_count;
	$mC_difference = $treat_mC_avg - $mock_mC_avg;

	return ($treat_mC_avg, $mock_mC_avg, $mC_difference, $treat_mC_min, $treat_mC_max, $mock_mC_min, $mock_mC_max, $mC_count);
}

sub fill_DIST_hash
{
        # Open the file
        my $inFile = $_[0];
        open (INFILE, $inFile) or die("Couldn't open $inFile, $!\n\n");
        #print "File = $inFile\n";

        # Reset the variables and get the new names for the cell line and treatment
        my ($locus_name, $k2p_dist, $last_locus, $last_dist) = ("", "", "", "");

        # Go through the file and put the expression information into the hash
        while(my $line = <INFILE>)
        {
                # Remove rogue new lines and split so manipulation of the elements is easy.
                # Then get and store the data I need in a multidimensional hash
		# The Kimura 2-parameter distance (k2p) is in the 6th column (0-based). Can change column for different distance parameters
                $line =~ tr/\r/\n/; chomp($line); chomp($line);
		if ($line =~ /locus/) { next; } # Skip header lines
                my @lineelements = "";
                @lineelements = split(/\t/, $line);
                $locus_name = $lineelements[0];
		$dist_hash{$locus_name} = $lineelements[6];

                # Store the variables for the return report to the user
                $last_locus = $locus_name;
                $last_dist = $lineelements[6];
        }

        # Return an update for the user
        return("Finished processing $inFile data.\nLast entry: Locus = $last_locus, distance = $last_dist.\n");
}

sub fill_ATAC_hash
{
        # Open the file
        my $inFile = $_[0];
        open (INFILE, $inFile) or die("Couldn't open $inFile, $!\n\n");
        #print "File = $inFile\n";

        # Reset the variables and get the new names for the cell line and treatment
        my ($cell_line, $treatment, $locus_name, $peak_score, $lastHashKeys, $last_score) = ("", "", "", "", "", "");
        ($cell_line, $treatment) = get_cellLine_treatment($inFile);
        #print "line = $cell_line\ntreat = $treatment\n";

        # Go through the file and put the expression information into the hash
        while(my $line = <INFILE>)
        {
                # Remove rogue new lines and split so manipulation of the elements is easy.
                # Then get and store the data I need in a multidimensional hash
                # The Kimura 2-parameter distance (k2p) is in the 6th column (0-based). Can change column for different distance parameters
                $line =~ tr/\r/\n/; chomp($line); chomp($line);
                my @lineelements = "";
                @lineelements = split(/\t/, $line);
                $locus_name = $lineelements[6];
		if ($lineelements[24] eq ".") { $lineelements[24] = "NA"; }
		#if (exists $atac_hash{$cell_line}{$treatment}{$locus_name}) # Code for taking the highest peak score
		#{
		#	if ($lineelements[24] < $atac_hash{$cell_line}{$treatment}{$locus_name}) { $lineelements[24] = $atac_hash{$cell_line}{$treatment}{$locus_name}; }
		#}
		if (exists $atac_hash{$cell_line}{$treatment}{$locus_name}) # Code for taking the peak score with the largest overlap
		{
			if ($lineelements[27] < $atac_overlap{$cell_line}{$treatment}{$locus_name}) 
			{
				$lineelements[24] = $atac_hash{$cell_line}{$treatment}{$locus_name};
				$lineelements[27] = $atac_overlap{$cell_line}{$treatment}{$locus_name};
			}
		}
                $atac_hash{$cell_line}{$treatment}{$locus_name} = $lineelements[24];
		$atac_overlap{$cell_line}{$treatment}{$locus_name} = $lineelements[27];

                # Store the variables for the return report to the user
                $lastHashKeys = "$cell_line $treatment $locus_name";
                $last_score = $lineelements[24];
        }

	#print "Hash values\n";
        #foreach my $derefkey (keys %{ $atac_hash{$cell_line}{$treatment} }) { print "Key = $derefkey. Value = $atac_hash{$cell_line}{$treatment}{$derefkey}\n"; }

        # Return an update for the user
        return("Finished processing $inFile ATAC-seq data.\nLast entry: Keys = $lastHashKeys, distance = $last_score.\n");	
}

sub get_cellLine_treatment # Call this to get the cell line and treatment names from a file name
{
	my ($filename, $linename, $treatname) = ("", "", ""); # Clear variables
	
	$filename = $_[0];	

	# Set the cell line name	
	if ($filename =~ /A2780/)	{ $linename = "A2780"; }
	if ($filename =~ /A549/)	{ $linename = "A549";  }
	if ($filename =~ /H1299/)	{ $linename = "H1299"; }
	if ($filename =~ /H23/)		{ $linename = "H23";   }
	if ($filename =~ /Hey/)		{ $linename = "Hey"; }
	if ($filename =~ /Kuramochi/)	{ $linename = "Kuramochi"; }
	if ($filename =~ /TykNu/) 	{ $linename = "TykNu"; }
	unless (exists $seen_lines{$linename}) { $seen_lines{$linename} = "Seen this line"; } 

	# Set the treatment name
	if ($filename =~ /treat1/) { $treatname = "Mock"; }
	if ($filename =~ /treat2/) { $treatname = "ITF"; }
	if ($filename =~ /treat3/) { $treatname = "AZA"; }
	if ($filename =~ /treat4/) { $treatname = "AZA/ITF"; }
	unless (exists $seen_treats{$treatname}) { $seen_treats{$treatname} = "Seen this treatment"; } 

	return ($linename, $treatname);
}
