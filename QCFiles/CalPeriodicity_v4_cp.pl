#!/usr/bin/perl
#Created by Prakitchai Chotewutmontri, 11 Aug 2015
#Postdoc, Barkan Lab, Institute of Molecular Biology, U of Oregon
#
#VERSION 3 - 11 Aug 2015
#  There is no prior version -- called ver 3 to keep it consistent with other scripts in ver 3
#
#NOTE 8/13/15
#  Find out the the UnionCDS NUC list file show the exon coordinate from
#      LastMin,LastMax;...;FirstMin,FirstMax
#  Though that it was
#      FirstMin,FirstMax;...;LastMin,LastMax
#  Change this program to reflex this and also made change to CP CDS list file
#      Alice_082714_CDS_Pt_withNum_v4_for_periodicity.txt ==same way ass the UnionCDS list instead of FirstMin,FirstMax;...;LastMin,LastMax
#
#VERSION cp
#  Try to count periodicity from footprint of 20 to 40
#  -In cp the offset of the length and the P-site is very linear
#RF len	offset
#20	12
#21	13
#22	14
#23	15
#24	16
#25	17
#26	18
#27	19
#28	20
#29	21
#30	22
#31	23
#32	24
#33	25
#34	26
#35	27
#36	28
#37	29
#38	30
#39	31
#40	32
#
#Version 4 - 6/20/18
#  - cp coordinate is always list first CDS to last CDS (exon 1 to n) in both positive and negative strand.
#  - re-define the left and right position for the samtools extraction -- instead of using START and STOP codon and assume it as min/max position (this give rps12 wrong because tran-splice), use the min and max values in the CDS coordinate instead.
#
################################################################################
#OBJECTIVE
#  Calculate the periodicity of the ribosome footprint
#  Considering the codon in the P-site of the ribosome, how many time the ribosome stop at
#position 0, 1, 2 nt of the translating codon.
#  Based on prior plotting of the 5' and 3' end of the ribosome footprints at START and STOP codon,
#we could use the footprint size information to link each ribosome footprint to a specific codon that it is translating.
#  eg. In cytosolic footprint, the 31-nt footprints have 5'end located at -13 nt position compare to 0 position of the codon.
#      In CP footprint, the 31-nt footprints have 5'end located at -23 nt position to 0 position of the codon.
#      In MT footprint, the 28-nt footprints have 5'end located at -19 nt position to 0 position of the codon (the evidence here is not that solid...seem to be different in START & STOP plots)
#
# Fix on 3/12/15
#   only count the read in the same orientation as the CDS
# Modified on 3/17/15
#   Add offset position to calculate RPKM
#   Sum of the offset cannot be larger than 72,117,66 nt for CP, Mito, Nuc, respectively.
#     These numbers are the smallest gene lengths. If offset is bigger then there will be nothing to count.
#   Fix the offset using strand info
#
# Fix on 4/28/15
#   use strand info of individual CDS to compare with read orientation
#
#INPUT1
#  A .bam file of the read alignment
#
#INPUT2
#  A file, my defined format containing CDS info
#
#OBSOLETE (old version)
#chromosome	gene	strand	attributes #CDS	CDS1start,stop;CDS2start,stop
#9	mRNA	+	ID=GRMZM2G581216;Name=GRMZM2G581216;biotype=transposable_element	1	19970,20092
#
#Version 3
#{chromosome}	{gene}	{attributes}	{#CDS}	{CDS1start,stop;CDS2start,stop}	{strand info}
#9	GRMZM2G581216	ID=GRMZM2G581216;Name=GRMZM2G581216;biotype=transposable_element	1	19970,20092	+
#
#OUTPUT
#  A new defined file
#-----------------------------
#gene	length(kbp)	#reads	RPKM
#psbA	1062	xxxx	yyyyyy
#

################################################################################
#MAIN PROGRAMS
################################################################################
use strict;
use warnings;

#Subroutine prototypes
sub CIGARMOVE;


#global var
my $minlen = 0;
my @readcount = ();
#set up table for periodicity count per RF length
#set up the offset table (key=RF length, value=offset)
my %periodicitycount = ();
my %offset = ();
my $minoffset = 12;
for (my $len=20; $len<=40; $len++) {
	$offset{$len}= $minoffset;
	$minoffset++;
	for(my $pos=0; $pos<=2; $pos++) {
		$periodicitycount{$len}{$pos} = 0;
	}
}

my $readcount = 0;

my $i = 0;
my $line = "";
my $samline = "";
my $samcount = 0;
my @temp = ();
my @chr = ();
my %chrlen = ();
my $CDScount = 0;

my $offset_start = 0;
my $offset_stop = 0;


my %gene_readcount = ();
my %gene_bp = ();
my %gene_bp_offset = ();
my $totalcount = 0;
my @gene_list = ();

my $cigar = ""; #keep CIGAR string operation
my %lengthcount = (); #This count all encounter reads in the defined region
                      #may be duplicate if a read shows up in more than one regions


#usage - 8/11/15
my $USAGE = "usage: $0 <input .BAM file> <input CDS list> ".
            "<output prefix>\n\n".
			"INFO: CYT footprint: peak at 31 nt, 13 nt offset\n".
		    "      CP footprint:  peak at 31 nt, 23 nt offset\n".
		    "      MT footprint:  peak at 28 nt, 19 nt offset\n\n";

#check argument - 8/11/15
unless (@ARGV == 3) {
	print $USAGE;
	exit -1;
}

#store argv - 8/11/15
my $bamfname = $ARGV[0];
my $cdsfname = $ARGV[1];
my $outfname = $ARGV[2];
my $outfnamedetail = $outfname."_detail.txt";
my $outfnamesummary = $outfname."_summary.txt";

#open infiles - 8/11/15
unless ( open(BAMFILE,"<$bamfname") ) {
	print "Can't open $bamfname\n\n";
	exit -1;
}
close (BAMFILE); #just checking if the file is here

unless ( open(CDSFILE,"<$cdsfname") ) {
	print "Can't open $$cdsfname\n\n";
	exit -1;
}

#open outfile, first time, overwrite it
unless ( open(OUTDETAIL,">$outfnamedetail") ) {
	print "Can't create $outfnamedetail\n\n";
}

#open outfile, first time, overwrite it
unless ( open(OUTSUM,">$outfnamesummary") ) {
	print "Can't create $outfnamesummary\n\n";
}

#(1) Go through every gene is the list - 8/11/15
#    -to avoid repeat counting, should use the union CDS list (all transcript variants of a gene are combined). But
#     the union CDS may not keep all the codon in the correct frame.
#        --going around this, I will only used the gene in the union CDS list where
#          the overall length of all the exons is a factor of 3. This is not a perfect fix but at least this kind of gene seems to have
#          the frame in the correct places.
#
#------------------------------------------------------------------------------

print OUTDETAIL "Determine periodicity from RF length of 22 to 40 with the offset of 12 to 32\n";
print OUTDETAIL "RF length";
for (my $a=22; $a<=40; $a++) {
	print OUTDETAIL "\t$a\t\t";
}
print OUTDETAIL "\n";
print OUTDETAIL "GENE";
for (my $a=22; $a<=40; $a++) {
	print OUTDETAIL "\tPos0\tPos1\tPos2";
}
print OUTDETAIL "\n";

#CDS list file is already opened

while ($line = <CDSFILE>) { #read line, each line = 1 mRNA
	#count CDS
	$CDScount++;

#chromosome	gene	strand	attributes #CDS	CDS1start,stop;CDS2start,stop
#9	mRNA	+	ID=GRMZM2G581216;Name=GRMZM2G581216;biotype=transposable_element	1	19970,20092

#Version 3
#{chromosome}	{gene}	{attributes}	{#CDS}	{CDS1start,stop;CDS2start,stop}	{strand info}

	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line);
	
	#get the CDS infomation
	my $g_count = 0;
	my $g_bp = 0;
	my $g_bp_offset = 0;
	
	#my $g_name = $temp[3];
	my $g_name = $temp[1]; #####version 3#####
	
	#print "$g_name\n";
	
	#my $numCDS = $temp[4];
	my $numCDS = $temp[3]; #####version 3#####
	
	#my @CDSinfo = split (/\;/, $temp[5]); #coordinate of the CDS
	my @CDSinfo = split (/\;/, $temp[4]); #coordinate of the CDS #####version 3#####
	
	my @strandinfo = split (/\;/, $temp[5]); #strand info of the CDS #####version 3#####
	
	#8/11/15 check the total length of CDS if it is a factor of three
	my $factorofthree = 0;
	my $totalCDSlength = 0;
	for (my $cdsi=0; $cdsi<$numCDS; $cdsi++) {
		my @exon = split(/\,/, $CDSinfo[$cdsi]);
		$totalCDSlength = $totalCDSlength + $exon[1] - $exon[0] + 1;
	}
	if ($totalCDSlength%3==0) {
		$factorofthree = 1;
		#print "\tpass factor of three\n";
	} else {
		print "$g_name\t***FAILED factor of three in 3nt periodicity analysis\n";
	}
	
	#8/11/15
	#keep the bin info for each CDS nt position
	############################################
	my %CDSntPOSbin = (); #keep the bin info for each CDS nt
	my $CDSntpos = 0;
	if ($factorofthree == 1) { #pass criteria for factor of three length
		if ($strandinfo[0] eq '+') { #strand is forward
			for (my $cdsi=0; $cdsi<$numCDS; $cdsi++) {
				my @exon = split(/\,/, $CDSinfo[$cdsi]);
				#print "\texon = $exon[0], $exon[1]\n";
				for (my $genomepos=$exon[0]; $genomepos<=$exon[1]; $genomepos++) {
					$CDSntpos++; #count the CDS nt position
					$CDSntPOSbin{$genomepos} = ($CDSntpos-1)%3;
					#print "genomepos=$genomepos\t frame=".(($CDSntpos-1)%3)."\n";
				}
			}
		} else { #negative strand
		#      LastMin,LastMax;...;FirstMin,FirstMax ===8/13/15 
			#for (my $cdsi=$numCDS-1; $cdsi>=0; $cdsi--) {
			for (my $cdsi=0; $cdsi<$numCDS; $cdsi++) { #version 4 6/20/18 
				my @exon = split(/\,/, $CDSinfo[$cdsi]);
				#print "\texon = $exon[1], $exon[0]\n";
				for (my $genomepos=$exon[1]; $genomepos>=$exon[0]; $genomepos--) {
					$CDSntpos++; #count the CDS nt position
					$CDSntPOSbin{$genomepos} = ($CDSntpos-1)%3;
					#print "genomepos=$genomepos\t frame=".(($CDSntpos-1)%3)."\n";
				}			
			}				
		}
	}
	
	#get START and STOP position
	#my $startgenomepos = 0;
	#my $stopgenomepos = 0;
	#if ($strandinfo[0] eq '+') {
	#	my @exon = split(/\,/, $CDSinfo[0]);
	#	$startgenomepos = $exon[0];
	#	@exon = split(/\,/, $CDSinfo[$numCDS-1]);
	#	$stopgenomepos = $exon[1];
	#} else { #      LastMin,LastMax;...;FirstMin,FirstMax ===8/13/15
	#	my @exon = split(/\,/, $CDSinfo[0]);
	#	#$stopgenomepos = $exon[0];
	#	$startgenomepos = $exon[1]; #VERSION 4
	#	#print "STOP = $stopgenomepos\n";
	#	@exon = split(/\,/, $CDSinfo[$numCDS-1]);
	#	#$startgenomepos = $exon[1];	
	#	$stopgenomepos = $exon[0];	#VERSION 4
	#	#print "START = $startgenomepos\n";	
	#}
	#calculate the offset position
	my $left=0; #min CDS position for SAMTOOLS
	my $right=0; #max CDS position for SAMTOOLS
	#if ($strandinfo[0] eq '+') {
	#	$left = $startgenomepos - 32;
	#	$right = $stopgenomepos - 12;
	#} else {
	#	$right = $startgenomepos + 32;
	#	$left = $stopgenomepos + 12;	
	#}
	#
	#version4 - below
	my @firstCDS = split(/\,/, $CDSinfo[0]);
	$left = $firstCDS[0];
	$right = $firstCDS[0];
	for (my $cdsi=0; $cdsi<$numCDS; $cdsi++) { #VERSION 4
		my @exon = split(/\,/, $CDSinfo[$cdsi]);
		if ($exon[0] < $left) { $left = $exon[0]; }
		if ($exon[0] > $right) { $right = $exon[0]; }
		if ($exon[1] < $left) { $left = $exon[1]; }
		if ($exon[1] > $right) { $right = $exon[1]; }
		
	}	
	#print "$g_name: $left to $right\n";
	
	#processing CDS info from offset parameters
	#print "\tGET region $left to $right\n";
	#get all the reads cover the whole gene CDS region
    open (SAMOUT, "samtools view ".
    			  "$bamfname ".
    			  "$temp[0]\:".
            	  $left."-".$right." |")
    or die "Can't run samtools: $!\n";
       		 
    #for every read in SAMOUT (each line)
        	
    #SAMTOOLS OUTPUT
	#SRR966474.103965279	0	Chr1	99979	8	30M	*	0	0	GGGCTTCATTTGTTTCTAGAGTGATTCATG	CCCFFFFFHHHHFIJJJJJJJFHJIIIIJJ	AS:i:-10	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0T0C28	YT:Z:UU
    #SRR966474.33107309	0	Chr1	99985	42	29M	*	0	0	CATTTGTTTCTAGAGTGATTCATGGATTC	@@@DDDDDFDDDDIB2<:IIG@HCDCFHB	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:29	YT:Z:UU
    #for every read in SAMOUT (each line)
    
    #set up table for periodicity count per RF length
	#set up the offset table (key=RF length, value=offset)
	my %genecount = ();
	for (my $len=20; $len<=40; $len++) {
		for(my $pos=0; $pos<=2; $pos++) {
			$genecount{$len}{$pos} = 0;
		}
	}
    my $totalsamoutline =0;
    my $countfailedlength =0;
    my $countfailedsign =0;
    my $countfailedposition =0;
    my $countgoodread =0;
	while ( $samline = <SAMOUT> ) {
		$totalsamoutline++;
		chomp ($samline); #remove endline char
		my @col = (); #clear array
		@col = split (/\t/, $samline);
		my $cnum = scalar (@col); #number of column, can be differ, good to track
			
		#re-calculate POS and READ length based on the CIGAR operation
		#The file may contain untrimmed sequences
		
		#checking the strand of the read
		my $sign = "+"; 
		my $bitcheck = $col[1] & 16; #check strand, - only if contain flag = 16
		if ($bitcheck == 16) { $sign = "-"; }

		#re-determine the 5' position of each read and also the length without soft-clipping region
		my ($newpos, $matchlen, $newlen);
		($newpos, $matchlen, $newlen) = CIGARMOVE($sign, $col[3], $col[5], $col[9] );
		
		#print "\t\tOLE POS=$col[3], NEW POS=$newpos\n";
		#newpos = 5' pos of both - or + strand read
		
		#Counting the read into each periodicity bin 0,1,2
		#1. only use the read with the specify footprint length
		if ($newlen >= 20 && $newlen <= 40) {
			#2. only use the read with the same strand as the gene
			if ($sign eq $strandinfo[0]) {
				#3. only use if the codon position is in the CDS
				my $ntposition;
				if ($strandinfo[0] eq '+') {
					$ntposition = $newpos + $offset{$newlen};
				} else {
					$ntposition = $newpos - $offset{$newlen};
				}
				if (exists($CDSntPOSbin{$ntposition})) { #found this position
					#count based on the codon postion
					$periodicitycount{$newlen}{$CDSntPOSbin{$ntposition}} = $periodicitycount{$newlen}{$CDSntPOSbin{$ntposition}} + 1;
					$genecount{$newlen}{$CDSntPOSbin{$ntposition}} = $genecount{$newlen}{$CDSntPOSbin{$ntposition}} + 1;
					$readcount++;
					$countgoodread++;
				} else {
					$countfailedposition++;
				}
			} else {
				$countfailedsign++;
			}
		} else {
			$countfailedlength++;
		}
		
	} #end SAMOUT while loop
	#print "\tJust analyzed total of $totalsamoutline reads\n";
	#print "\tRead with wrong length $countfailedlength\n";
	#print "\tRead with wrong strand $countfailedsign\n";
	#print "\tRead with wrong position $countfailedposition\n";
	#print "\tGood reads = $countgoodread\n";
    close(SAMOUT);
    print OUTDETAIL "$g_name";
    for (my $l=20; $l<=40; $l++) {
		for(my $pos=0; $pos<=2; $pos++) {
			print OUTDETAIL "\t$genecount{$l}{$pos}";
		}
	}
	print OUTDETAIL "\n";
    			  
} #end while <CDSFILE>
			


#(2) Write output
#------------------------------------------------------------------------------
#general info
print OUTSUM "FINAL COUNT:\t$CDScount\tCDS genes\n";
print OUTSUM "Footprint length from 20 to 40 nt\n";
print OUTSUM "Offset from 12 to 32 nt\n";
print OUTSUM "Total found reads:\t$readcount\treads\n";

#periodicity distribution
print OUTSUM "Read length\tPos0\tPOS1\tPOS2\n";
for (my $l=20; $l<=40; $l++) {
	print OUTSUM "$l";
	for(my $pos=0; $pos<=2; $pos++) {
		print OUTSUM "\t$periodicitycount{$l}{$pos}";
	}
	print OUTSUM "\n";
}

close (CDSFILE);
close (OUTDETAIL);
close (OUTSUM);

exit (0);



############################################################################
#Subroutine programs
############################################################################
#
sub CIGARMOVE () {
	#input value
	my $strand = $_[0]; #direction of the read (+ or -)
	my $pos = $_[1]; #left position
	my $cig = $_[2]; #cigar string
	my $seq = $_[3]; #seq
	
	my $len = length ($seq);
	my $newpos = 0;
	my $mlen = 0; #length of matched sequence based on the reference genome
	my $newlen = 0; #read length without soft clipping
	
	#Need to find the first N-terminal position of alignment match
	
	#CIGAR operation
	#M  alignment match (can be a sequence match or mismatch)
	#I  insertion to the reference
	#D  deletion from the reference
	#N  skipped region from the reference
	#S  soft clipping (clipped sequences present in SEQ)
	#H  hard clipping (clipped sequences NOT present in SEQ)
	#P  padding (silent deletion from padded reference)
	#=  sequence match
	#X  sequence mismatch
	
	#The first position of the alignment
		#This is what I think:
	    #- cannot be an Insertion, Deletion (Gap) or Skipped region because it will be impossible to determine POS on the ref genome?
	    #- Hard clipping is not show anyways then it cannot be used to determine POS
	    #- Padding??? we don't used padded reference
	    #- The First Position CAN ONLY BE: M,S,=,X
	    
	    #- The LAST Position CAN ONLY BE: M,I,S,=,X
	
	#READ and keep CIGAR string
	my @val = ();
	my @letter = ();
	while ($cig =~ /(\d+)([MIDNSHP=X])/g) {
		push (@val, $1); #keep value in the array
		push (@letter, $2); #keep operation letter in the array
	}
	
	#Calculate the length of matched in term of the reference genome
	for (my $p = 0; $p < scalar(@val); $p++) {
		if (($letter[$p] eq 'M') || ($letter[$p] eq 'D') || ($letter[$p] eq '=') ||
		    ($letter[$p] eq 'X') || ($letter[$p] eq 'N')) {
			$mlen = $mlen + $val[$p];
		}
	}
	
	if ($strand eq '+') {
		#the position in .BAM is the first matched position
		#for the forward read, this is the 5' end position of the matched region
		$newpos = $pos;
	} elsif ($strand eq '-') {
		#the position in .BAM is the first matched position
		#for the reverse read, this is the 3' end position of the matched region
		#the 5'end = pos + matched length -1
		$newpos = $pos+$mlen-1;
	}
	
	#Calculate the length from M,I,=,X =>read length (S is not part of read)
	for (my $r = 0; $r < scalar(@val); $r++) {
		if (($letter[$r] eq 'M') || ($letter[$r] eq 'I') || ($letter[$r] eq '=') ||
		    ($letter[$r] eq 'X')) {
			$newlen = $newlen + $val[$r];
		}
	}
    return ($newpos,$mlen,$newlen);
}

