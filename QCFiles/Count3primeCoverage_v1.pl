#!/usr/bin/perl
#Created by Prakitchai Chotewutmontri, 11 May 2014
#Postdoc, Barkan Lab, Institute of Molecular Biology, U of Oregon
#
#Original file = Count3primeCoverage_ChrPt_v3.pl
#
################################################################################
#OBJECTIVE
#  Determine the distribution of 3'end of the reads along the genome
#  This version will work both for Mt and Pt
#
#INPUT1
#  A .bam file of the read alignment
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
my $ntcount = 0;
my $readcount = 0;


my $cigar = ""; #keep CIGAR string operation
my %lengthcount = (); #This count all encounter reads in the defined region
                      #may be duplicate if a read shows up in more than one regions

my %rlength = ();
my %cover = ();


#usage
my $USAGE = "usage: $0 <input .BAM file> <output filename> <wig filename>\n\n";

#check argument
unless (@ARGV == 3) {
	print $USAGE;
	exit -1;
}

#store argv
my $bamfname = $ARGV[0];
my $outfname = $ARGV[1];
my $wigfname = $ARGV[2];
my $wigextname = "";

#open infiles
unless ( open(BAMFILE,"<$bamfname") ) {
	print "Can't open $bamfname\n\n";
	exit -1;
}
close (BAMFILE); #just checking if the file is here


#open outfile, first time, overwrite it
unless ( open(OUTFILE,">$outfname") ) {
	print "Can't create $outfname\n\n";
}

$wigextname = $wigfname."_20-24.wig"; #for reads with length 20-24
#open outfile, first time, overwrite it
unless ( open(WIG1,">$wigextname ") ) {
	print "Can't create $outfname\n\n";
}
$wigextname = $wigfname."_26-34.wig"; #for reads with length 26-34
#open outfile, first time, overwrite it
unless ( open(WIG2,">$wigextname ") ) {
	print "Can't create $outfname\n\n";
}
$wigextname = $wigfname."_sum.wig"; #for reads with length 16-40
#open outfile, first time, overwrite it
unless ( open(WIG3,">$wigextname ") ) {
	print "Can't create $outfname\n\n";
}
#(0) Check .bam file what genome is this
#-- need to have .bai file 
#-- need to be single chromosome
#------------------------------------------------------------------------------
#retrive reads that cover each position
open (BAMIDXSTATS, "samtools idxstats $bamfname |")
or die "Can't run samtools: $!\n";
#Mt	569630	1085155	0
#*	0	0	13555311
#
#X86563	140384	3477745	0
#*	0	0	16121217
my $idxstats = <BAMIDXSTATS>;
my @tmp = split (/\t/, $idxstats);
my $chr = $tmp[0];
my $clength = $tmp[1];
close (BAMIDXSTATS);

#(1) Go through every nt in the chromosome
#------------------------------------------------------------------------------

#print outfile header
print OUTFILE "position";
for ($i = 16; $i<=40; $i++) {
	print OUTFILE "\t$i";
}
print OUTFILE "\n";

#print WIG files header
print WIG1 "track type=wiggle_0\t\n";
print WIG1 "variableStep chrom=$chr span=1\t\n";
print WIG2 "track type=wiggle_0\t\n";
print WIG2 "variableStep chrom=$chr span=1\t\n";
print WIG3 "track type=wiggle_0\t\n";
print WIG3 "variableStep chrom=$chr span=1\t\n";

my $minnt=0;
my $maxnt=0;

for (my $n = 1; $n <= int($clength/5000)+1; $n++) { #read line, each line = 1 mRNA
	$minnt = $maxnt+1;
	$maxnt = $maxnt+5000;
	if ($maxnt > $clength){
		$maxnt = $clength;
	}

	#retrive reads that cover each position
	open (SAMOUT, "samtools view ".
                  "$bamfname ".
              	  "$chr\:".$minnt."-".$maxnt." |")
    or die "Can't run samtools: $!\n";
        
    #for every read in SAMOUT (each line)
	while ( $samline = <SAMOUT> ) {
		chomp ($samline); #remove endline char
		my @col = (); #clear array
		@col = split (/\t/, $samline);
		my $cnum = scalar (@col); #number of column, can be differ, good to track
			
		#re-calculate POS and READ length based on the CIGAR operation
		#The file may contain untrimmed sequences
		my ($newpos, $matchlen, $newlen);
		
		my $sign = "+";
		my $bitcheck = $col[1] & 16; #check strand, - only if contain flag = 16
		if ($bitcheck == 16) { $sign = "-"; }
		
		($newpos, $matchlen, $newlen) = CIGARMOVE($sign, $col[3], $col[5], $col[9] );
		
		#newpos = 5' pos of either - or + strand read
		
		my $pos3end = 0;
		if ($sign eq "+") {
			$pos3end = $newpos+$matchlen-1;
		} else { #$sign = "-"
			#$pos3end = $newpos-$matchlen+1;
			$pos3end = $col[3];
		}
		
		#store data -- read count
		if (exists $cover{$pos3end}{$newlen}) {
			$cover{$pos3end}{$newlen} = $cover{$pos3end}{$newlen}+1;
		}else{
			$cover{$pos3end}{$newlen} = 1;
		}
		#store readlenth data
		if (exists $rlength{$newlen}) {
			$rlength{$newlen} = $rlength{$newlen}+1;
		}else{
			$rlength{$newlen} = 1;
		}	
		$readcount++;
	} #end SAMOUT while loop
		
	close (SAMOUT);
				
	#print outfile
	for (my $pi=$minnt; $pi<=$maxnt; $pi++){
		print OUTFILE "$pi";
		print WIG1 "$pi";
		print WIG2 "$pi";
		print WIG3 "$pi";
		my $wig1sum = 0; #read length 20-24
		my $wig2sum = 0; #read length 26-34
		my $wig3sum = 0; #all read length
		for (my $nr = 16; $nr<=40; $nr++) {
			my $count = 0;
			if (exists $cover{$pi}{$nr} ) {
				$count = $cover{$pi}{$nr};
				$wig3sum = $wig3sum + $count; #add all the count to WIG3
				if (($nr ge 20) && ($nr le 24)){
					$wig1sum = $wig1sum + $count;
				}
				if (($nr ge 26) && ($nr le 34)){
					$wig2sum = $wig2sum + $count;
				}
			}
			print OUTFILE "\t$count";
		}
		print WIG1 "\t$wig1sum\n";
		print WIG2 "\t$wig2sum\n";
		print WIG3 "\t$wig3sum\n";
		print OUTFILE "\n";
	}
	
	print "At position $maxnt, there are $readcount counted reads\n";

} #end for loop --nt position

			


#(2) Write output
#------------------------------------------------------------------------------
#general info
print "\n\nFINAL COUNT: $maxnt positions, found $readcount reads\n";
print "\n\n\n";

close (OUTFILE);
close (WIG1);
close (WIG2);
close (WIG3);

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


