#!/usr/bin/perl
#Created by Prakitchai Chotewutmontri, 2 Mar 2015
#Postdoc, Barkan Lab, Institute of Molecular Biology, U of Oregon
#
#VERSION 3 - 28 April 2015
#  Update to recognize the CDS list input file version 3.
#  This format contain strand info for individual CDS
#
#  Update chromosome name in CDS list file to matched that of .bam file
#
################################################################################
#OBJECTIVE
#  Determine read length distribution within CDS
#     - using .bam (read aligned to the genome)
#     - CDS coordinate file
#
#  Modified from CountReadLengthInCDS_chrMt_v2.pl
#	  - call CICARMOVE() with strand determined from FLAG value of SAM file
#	  - re-write CICARMOVE()
#  Fixed 3/12/15
#	  - calculate position first match of the negative strand read
#     - fix the check position of read to overlap with CDS..differ for positive and negative strand
#     - only count the reads in the same strand as CDS
#
#VERSION 3 - 28 April 2015
#  Update to recognize the CDS list input file version 3.
#  This format contain strand info for individual CDS
#
#INPUT1
#  A .bam file of the read alignment
#
#INPUT2
#  A file, my defined format containing CDS info
#
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
#position	count
#-24		0
#-23		0
#...
#0			5000
#1			300
#2			350
#...
#62			57
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

my $i = 0;
my $line = "";
my $samline = "";
my @temp = ();
my @chr = ();
my %chrlen = ();
my $CDScount = 0;
my $readcount = 0;


my $cigar = ""; #keep CIGAR string operation
my %lengthcount = (); #This count all encounter reads in the defined region
                      #may be duplicate if a read shows up in more than one regions


#usage
my $USAGE = "usage: $0 <input .BAM file> <input CDS list> ".
            "<output filename>\n\n";

#check argument
unless (@ARGV == 3) {
	print $USAGE;
	exit -1;
}

#store argv
my $bamfname = $ARGV[0];
my $cdsfname = $ARGV[1];
my $outfname = $ARGV[2];

#open infiles
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
unless ( open(OUTFILE,">$outfname") ) {
	print "Can't create $outfname\n\n";
}


#(1) Go through every mRNA in the list, count read length
#       this should include all the transcript variants!!?!?!
#------------------------------------------------------------------------------

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
	my $gene_readcount = 0;
	#my $numCDS = $temp[4];
	my $numCDS = $temp[3]; #####version 3#####
	
	#my @CDSinfo = split (/\;/, $temp[5]); #coordinate of the CDS
	my @CDSinfo = split (/\;/, $temp[4]); #coordinate of the CDS #####version 3#####
	
	my @strandinfo = split (/\;/, $temp[5]); #strand info of the CDS #####version 3#####
	
	#print "Gene=$temp[3]\tNumberCDS=$numCDS\n";

		
	#Go through each CDS===================================================
	for (my $ci = 0; $ci < $numCDS; $ci++) {
		my %readid = (); #keep the id of reads, count once per gene.
		my @aCDS = split(/\,/, $CDSinfo[$ci]);
		#print "\t\tCDS = $aCDS[0]-$aCDS[1]\n";
		
		#SAMTOOLS OUTPUT
		#SRR966474.103965279	0	Chr1	99979	8	30M	*	0	0	GGGCTTCATTTGTTTCTAGAGTGATTCATG	CCCFFFFFHHHHFIJJJJJJJFHJIIIIJJ	AS:i:-10	XN:i:0	XM:i:2	XO:i:0	XG:i:0	NM:i:2	MD:Z:0T0C28	YT:Z:UU
        #SRR966474.33107309	0	Chr1	99985	42	29M	*	0	0	CATTTGTTTCTAGAGTGATTCATGGATTC	@@@DDDDDFDDDDIB2<:IIG@HCDCFHB	AS:i:0	XN:i:0	XM:i:0	XO:i:0	XG:i:0	NM:i:0	MD:Z:29	YT:Z:UU
        
        #get all read cover this CDS region
		open (SAMOUT, "samtools view ".
              		  "$bamfname ".
              		  "$temp[0]\:".
              		  $aCDS[0]."-".$aCDS[1]." |")
        or die "Can't run samtools: $!\n";
        
        #for every read in SAMOUT (each line)
		while ( $samline = <SAMOUT> ) {
			chomp ($samline); #remove endline char
			my @col = (); #clear array
			@col = split (/\t/, $samline);
			my $cnum = scalar (@col); #number of column, can be differ, good to track
			
			my $flag = "+";
			my $bitcheck = $col[1] & 16; #check strand, - only if contain flag = 16
			if ($bitcheck == 16) { $flag = "-"; }
			
			#only use the same strand as the gene CDS
			#if ($flag eq $temp[2]) {
			if ($flag eq $strandinfo[$ci]) { #####version 3#####
			
				#re-calculate POS and READ length based on the CIGAR operation
				#The file may contain untrimmed sequences
				my ($newpos, $newlen) = CIGARMOVE( $flag, $col[3], $col[5], $col[9] );
				
				#check if the aligned region is still overlap with the CDS region
				my $overlap = 0; #not detect overlap yet
				my $leftpos = 0;
				my $rightpos = 0;
				if ( $flag eq "+" ) {
					$leftpos = $newpos;
					$rightpos = $newpos+$newlen;
				} else { # - strand
					$leftpos = $newpos-$newlen;
					$rightpos = $newpos;			
				}
				for (my $n=$leftpos; $n<=$rightpos; $n++) {
					if (($n >=$aCDS[0]) && ($n <= $aCDS[1])) { #pos is in CDS
						$overlap = 1; #YES it is overlapped
						last; #no longer need to check the rest
					}
				}
			
				#if it is overlapped
				if ($overlap) {
					#check if this read is used?
					if ( !(exists($readid{$col[0]}))) {
						#not exists = new read
						$readcount++;
						$gene_readcount++;
						#print "$readcount FOUND\n";
						#count the length
						if ( exists ($lengthcount{$newlen})){
							$lengthcount{$newlen} = $lengthcount{$newlen} +1;
						} else {
							$lengthcount{$newlen} = 1;
						}
						#mark it in the hash
						$readid{$col[0]} = 1; #this id was used
					} else { #else if exists....don't count
						#print "this read was used\n";
					}
				} # else if not overlap....not use
			} #end check same strand as gene
		} #end while <SAMOUT>
		
		close (SAMOUT);
	
	} #end go through each CDS loop
	
	#print "Gene=$temp[3]\thas\t$gene_readcount\treads\n";
	
	if ($CDScount % 200 == 0) {
		#print "At $CDScount CDSs there are $readcount counted reads\n";
	}
} #end while <CDSFILE>
			


#(2) Write output
#------------------------------------------------------------------------------
#general info
#print "\n\nFINAL COUNT: $CDScount CDSs, found $readcount reads\n";
print OUTFILE "FINAL COUNT: $CDScount CDSs, found $readcount reads\n";

#print length distribution
#print "Length\tCount\n";
print OUTFILE "Length\tCount\n";

for my $ll ( sort keys %lengthcount ) {
	#print "$ll\t$lengthcount{$ll}\n";
	print OUTFILE "$ll\t$lengthcount{$ll}\n";
}
#print "\n\n\n";

close (CDSFILE);
close (OUTFILE);

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
	my $newlen = 0;
	
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
	
	#Calculate the length from M,I,=,X =>read length (S is not part of read)
	for (my $r = 0; $r < scalar(@val); $r++) {
		if (($letter[$r] eq 'M') || ($letter[$r] eq 'I') || ($letter[$r] eq '=') ||
		    ($letter[$r] eq 'X')) {
			$newlen = $newlen + $val[$r];
		}
	}

	#determine the 5' end position (of the read) of the matched region
	my $posadj = 0;
	
	if ($strand eq '+') {
		$newpos = $pos; #in positive strand: the pos = the first matched position
	} elsif ($strand eq '-') {
		#in negative strand: the pos = the last matched position
		$newpos = $pos+$newlen-1; #fixed on 3/12/15
	}

    return ($newpos,$newlen);
}

