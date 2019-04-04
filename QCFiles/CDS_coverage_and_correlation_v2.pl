#!/usr/bin/perl
#Created by Prakitchai Chotewutmontri, 15 June 2017
#Postdoc, Barkan Lab, Institute of Molecular Biology, U of Oregon
#
#original script: percent_CDS_coverage_v3.pl
#
#CDS_coverage_and_correlation_v1.pl
#
#Version 2 -- 6/19/18
# - change for Nick maize v4 pipeline
#
################################################################################
#OBJECTIVE
#  Quality control for read coverage in the CDS
# - check number of nucleotide in the CDS (or defined region) that are covered by read (with cut-off = 0 as build-in parameter)_
# - check the correlation coefficient of each CDS (ore defined region)
#
#Strategy
# - for the correlation coefficient, predetermined the strand specific read coverage from the maize tissue in the light (midday rep2 and 3,
#and reilluminated rep 2 and 3)--computed into averaged read coverage per million (total nuclear CDS reads). This files are defined as standard read
#coverage. These samples were prepared using BIOO sRNA kit version 3.
# - generate the strand specific read coverage of the input, compare to the standard coverage to compute the correlation coefficient. At the
#same time also keep track of the read covered region. -- THIS STEP NEED "samtools"
#
#INPUT1
#  A .bam file of the read alignment
#
#INPUT2
#  A file, my defined format containing CDS info
#
#INPUT3
#  Standard read coverage for FORWARD strand
#<nt position><read per million>
#
#INPUT4
#  Standard read coverage for REVERSE strand
#<nt position><read per million>
#
#OUTPUT1 -- number of covered nucleotide positions
#<gene><CDS_length><covered_CDS_length>
#
#OUTPUT2 -- correlation coefficient 
#<gene><Pearson correlation coefficient>

################################################################################
#MAIN PROGRAMS
################################################################################
use strict;
use warnings;

#defined coverage cut-off
my $COVERAGECUTOFF = 0; #this is the value cut-off for coverage

#keep the standard coverage values
my %forval = (); #key = nucleotide position, value = read coverage per million
my %revval = ();

#input coverage
my %infor = ();
my %inrev = ();

#global var
my $line = "";
my @temp = ();


my %nt_coverage = ();

my $minlen = 0;
my @readcount = ();

my $i = 0;
my $samline = "";
my $samcount = 0;
my @chr = ();
my %chrlen = ();
my $CDScount = 0;

my %gene_readcount = ();
my %gene_inicount = ();
my %gene_stallcount = ();
my %gene_bp = ();
my %gene_bp_offset = ();
my $totalcount = 0;
my @gene_list = ();

my $cigar = ""; #keep CIGAR string operation
my %lengthcount = (); #This count all encounter reads in the defined region
                      #may be duplicate if a read shows up in more than one regions

my %stallpos = ();


my %gene_stall_spot= ();
my %gene_stall_reads=();

#usage
my $USAGE = "usage: $0 <.BAM file> <CDS list> <Standard FORWARD strand cover> <Standard REVERSE strand cover>\n\n";

#check argument
unless (@ARGV == 4) {
	print $USAGE;
	exit -1;
}

#store argv
my $bamfname = $ARGV[0];
my $cdsfname = $ARGV[1];
my $forfname = $ARGV[2];
my $revfname = $ARGV[3];
my $outfname1 = $bamfname.".QC_nt-coverage.txt"; ##version 2
my $outfname2 = $bamfname.".QC_Pearson_correlation.txt";  ##version 2

#open infiles
unless ( open(BAMFILE,"<$bamfname") ) {
	print "Can't open $bamfname\n\n";
	exit -1;
}
close (BAMFILE); #just checking if the file is here

unless ( open(CDSFILE,"<$cdsfname") ) {
	print "Can't open $cdsfname\n\n";
	exit -1;
}

unless ( open(FORFILE,"<$forfname") ) {
	print "Can't open $forfname\n\n";
	exit -1;
}

unless ( open(REVFILE,"<$revfname") ) {
	print "Can't open $revfname\n\n";
	exit -1;
}

#open outfile, first time, overwrite it
unless ( open(COVERFILE,">$outfname1") ) {
	print "Can't create $outfname1\n\n";
}

#open outfile, first time, overwrite it
unless ( open(PEARSONFILE,">$outfname2") ) {
	print "Can't create $outfname2\n\n";
}

#------------------------------------------------------------------------------
#	(1) Read in the standard coverage files
#------------------------------------------------------------------------------
while ($line = <FORFILE>) { #read line, each line = 1 mRNA
	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line);
	$forval{$temp[0]} = $temp[1];
}
close (FORFILE);

while ($line = <REVFILE>) { #read line, each line = 1 mRNA
	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line);
	$revval{$temp[0]} = $temp[1];
}
close (REVFILE);

#------------------------------------------------------------------------------
#	(2) Generate the strand specific coverage of the input file and store them in the memory
#------------------------------------------------------------------------------
# split strand
my $forbam = $bamfname.".forward.bam";
my $revbam = $bamfname.".reverse.bam";
my $forextrun = "samtools view -b -F 16 $bamfname > $forbam";
my $revextrun = "samtools view -b -f 16 $bamfname > $revbam";
#execute samtools
system($forextrun);
system($revextrun);
#make sure the bam files are made
my $fbamsize_old = (stat $forbam)[7];
my $fbamcheck  = 0;
while ( $fbamcheck == 0 ) { ##running this loop until the file finish writing
	sleep 5;
	my $fbamsize_new = (stat $forbam)[7];
	if ($fbamsize_new == $fbamsize_old) {
		$fbamcheck++;
	} else {
		$fbamsize_old = $fbamsize_new;
	}
}
#make sure the bam files are made
my $rbamsize_old = (stat $revbam)[7];
my $rbamcheck  = 0;
while ( $rbamcheck == 0) { ##running this loop until the file finish writing
	sleep 5;
	my $rbamsize_new = (stat $revbam)[7];
	if ($rbamsize_new == $rbamsize_old) {
		$rbamcheck++;
	} else {
		$rbamsize_old = $rbamsize_new;
	}
}

#get the depth
my $fordepth = $bamfname.".forward.depth";
my $revdepth = $bamfname.".reverse.depth";
my $fordepthrun = "samtools depth -a -d 1000000 $forbam > $fordepth";
my $revdepthrun = "samtools depth -a -d 1000000 $revbam > $revdepth";
#execute samtools
system($fordepthrun);
system($revdepthrun);

#read in the depth into memory
my $fsize_old = (stat $fordepth)[7];
my $fcheck  = 0;
while ( $fcheck == 0) { ##running this loop until the file finish writing
	sleep 5;
	my $fsize_new = (stat $fordepth)[7];
	if ($fsize_new == $fsize_old) {
		$fcheck++;
	} else {
		$fsize_old = $fsize_new;
	}
}
unless ( open(DEPTHFOR,"<$fordepth") ) {
	print "Can't open $fordepth\n\n";
	exit -1;
}
while ($line = <DEPTHFOR>) { #read line, each line = 1 mRNA
	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line); #<chromosome><position><depth>
	$infor{$temp[1]} = $temp[2];
}
close (DEPTHFOR);

my $rsize_old = (stat $fordepth)[7];
my $rcheck  = 0;
while ( $rcheck == 0) { ##running this loop until the file finish writing
	sleep 5;
	my $rsize_new = (stat $revdepth)[7];
	if ($rsize_new == $rsize_old) {
		$rcheck++;
	} else {
		$rsize_old = $rsize_new;
	}
}
unless ( open(DEPTHREV,"<$revdepth") ) {
	print "Can't open $revdepth\n\n";
	exit -1;
}
while ($line = <DEPTHREV>) { #read line, each line = 1 mRNA
	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line); #<chromosome><position><depth>
	$inrev{$temp[1]} = $temp[2];
}
close (DEPTHREV);

##remove the intermediate files
unlink ($forbam);
unlink ($revbam);
unlink ($fordepth);
unlink ($revdepth);


#------------------------------------------------------------------------------
#	(3) Go to every gene, count the covered nucleotide & Pearson correlation
#   (4) write out the output
#------------------------------------------------------------------------------
#
#Pearson correlation coefficient = r = (n*(sum(x*y)) - sum(x)*sum(y))/square root((n*sum(x^2)-(sum(x))^2)*(n*sum(y^2)-(sum(y))^2))
#
#
#print "\n\nFINAL COUNT: $CDScount CDS genes, found $totalcount reads\n";
print COVERFILE "Gene\tCDS length (nt)\tCovered CDS length (nt)\n";
print PEARSONFILE "Gene\tPearson correlation coefficient\n";

#CDS list file is already opened
while ($line = <CDSFILE>) { #read line, each line = 1 mRNA
	#count CDS
	$CDScount++;

#{chromosome}	{gene}	{attributes}	{#CDS}	{CDS1Left,Right;CDS2Left,Right}	{strand info}
#9	mRNA	+	ID=GRMZM2G581216;Name=GRMZM2G581216;biotype=transposable_element	1	19970,20092

	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line);
	
	#get the CDS infomation
	my $CDS_nt = 0;
	my $compare_nt = 0;
	my $covered_nt = 0;
	my $g_name = $temp[1]; #####version 3#####
	my $numCDS = $temp[3]; #####version 3#####
	my @CDSinfo = split (/\;/, $temp[4]); #coordinate of the CDS #####version 3#####	
	my @strandinfo = split (/\;/, $temp[5]); #strand info of the CDS #####version 3#####
	
	#my $num = 0; #same value as $CDS_nt
	my $sumx = 0; #input
	my $sumy = 0; #standard
	my $sumxy = 0;
	my $sumxsquare = 0;
	my $sumysquare = 0;
	
	# go through each CDS in the gene and count the read
	#####################################
	for (my $ci = 0; $ci < $numCDS; $ci++) {
		my @aCDS = split(/\,/, $CDSinfo[$ci]);
		$CDS_nt = $CDS_nt + $aCDS[1] - $aCDS[0] + 1; #count bp
		
		my $CDSstrand = $strandinfo[$ci];
		
		#go to every nt in CDS
		for (my $nt=$aCDS[0]; $nt<=$aCDS[1]; $nt++) {
			my $readcount = 0;
			my $xval = 0;
			my $yval = 0;
			
			if ($CDSstrand eq '+') {
				#check the std val and keep in the list
				if (exists($forval{$nt})) {
					$yval = $forval{$nt};
					$compare_nt++;
				} #not exists, then value =0 -> do nothing form Pearson calculation
				#check the inout val and keep in the list
				if (exists($infor{$nt})) {
					$readcount = $infor{$nt};
					$xval = $infor{$nt};
				} #not exists, then value =0 -> do nothing form Pearson calculation
			} else { #negative strand
				#check the std val and keep in the list
				if (exists($revval{$nt})) {
					$yval = $revval{$nt};
					$compare_nt++;
				} #not exists, then value =0 -> do nothing form Pearson calculation
				#check the inout val and keep in the list
				if (exists($inrev{$nt})) {
					$readcount = $inrev{$nt};
					$xval = $inrev{$nt};
				} #not exists, then value =0 -> do nothing form Pearson calculation
			}
			# cal the Pearson value
			$sumx = $sumx + $xval; #input
			$sumy = $sumy + $yval; #standard
			$sumxy = $sumxy + ($xval*$yval);
			my $xx = $xval*$xval;
			my $yy = $yval*$yval;
			$sumxsquare = $sumxsquare + $xx;
			$sumysquare = $sumysquare + $yy;
			#print "nt\t$nt\t$xval\t$yval\t$sumx\t$sumy\t$sumxy\t$sumxsquare\t$sumysquare\n";
			
			#check if this nt pass cover cutoff
			if ($readcount > $COVERAGECUTOFF) {
        		$covered_nt++; #count covered nt position
        	}
        } #for every nucleotide

	} #for every CDS
    
    #calculate Pearson correlation
    #r = (n*(sum(x*y)) - sum(x)*sum(y))/square root((n*sum(x^2)-(sum(x))^2)*(n*sum(y^2)-(sum(y))^2))
    my $denominator = sqrt((($compare_nt*$sumxsquare)-($sumx*$sumx))*(($compare_nt*$sumysquare)-($sumy*$sumy)));
    my $rcoef;
    if ( $denominator == 0) {
    	$rcoef = "N/A";
    } else {
    	$rcoef = ($compare_nt*$sumxy - $sumx*$sumy)/$denominator;
    }    
    
    # write out outfile
	#####################################
	print COVERFILE "$g_name\t$CDS_nt\t$covered_nt\n";
	print PEARSONFILE "$g_name\t$rcoef\n";
	
} #end while <CDSFILE>

			
close (CDSFILE);
close (COVERFILE);
close (PEARSONFILE);


exit (0);