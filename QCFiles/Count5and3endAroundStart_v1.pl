#!/usr/bin/perl
#Created by Prakitchai Chotewutmontri, 11 May 2015
#Postdoc, Barkan Lab, Institute of Molecular Biology, U of Oregon
#
#Original program = Count3endAroundStart_v4.pl
#
#THIS PROGRAM ASSUME THAT THE CDS COORDINATES ARE ARRANGED FROM CDS1-to-CDSn
#  With the unionCDS -> this arrangement may be off...they get order from min to max position.
#  Most positive strand will be correct but ALL negative strand will be WRONG
#
####CAUTION!!!!!############
#ONLY use NON-UNION CDS with this program!!!!!!!!!!!!!!!!!!!!!!!!
#
################################################################################
#OBJECTIVE
#  Count number of 5 and 3 prime end around AUG start site
#		- user define the region
#			eg. from -24 to 62
#		- check if the CDS cover the region, if not filter it out
#			eg. CDS should be at 0 to 62 = 63 nt
#				The ribosome have offset of 12-nt.
#				For the 62th pos to get cover, the ribosome offset will be 62 to 73.
#				And that will be the 74, 75, 76 pos for the codon.
#				The last codon in CDS = STOP codon
#				Thus CDS need to be at least 76+1+3 = 80 nt or (the last pos + 18 nt)
#				!!if not filter out, it may distort the results.
#		- Hopefully, all the -24 to -1 position are not part of other genes!?! 
#		- ONLY USE READS with a defined length!!
#		- only count each START/STOP region once. Gene with multiple transcripts get count once at each START/STOP.
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
my %readcount5p = (); #keep read count based on the 5end locations of the reads, key = read length
my %readcount3p = (); #keep read count based on the 3end locations of the reads, key = read length

my $i = 0;
my $j = 0;
my $line = "";
my $samline = "";
my @temp = ();
my @chr = ();
my %chrlen = ();
my $gcount = 0;
my $nummatch5p = 0;
my $nummatch3p = 0;

my %USECOOR  = (); #keep track of counted coordinates #####version 4#####


my $cigar = ""; #keep CIGAR string operation

#usage
my $USAGE = "usage: $0 <input .BAM file> <input CDS list> ".
            "<5' outfile> <5' out sum file> ".
            "<3' outfile> <3' out sum file> ".
            "<negative pos, no sign> <positive pos, from zero> ".
            "<min read length> <max read length>\n\n";

#check argument
unless (@ARGV == 10) {
	print $USAGE;
	exit -1;
}

#store argv
my $bamfname = $ARGV[0];
my $cdsfname = $ARGV[1];
my $outfname5p = $ARGV[2];
my $outsumfname5p = $ARGV[3];
my $outfname3p = $ARGV[4];
my $outsumfname3p = $ARGV[5];
my $lpos = $ARGV[6];
my $rpos = $ARGV[7];
my $minreadlength = $ARGV[8];
my $maxreadlength = $ARGV[9];
my $totalpos = $lpos+$rpos+1;

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

#open outfile, first time, append it
unless ( open(OUTFILE5P,">>$outfname5p") ) {
	print "Can't create $outfname5p\n\n";
}

#open outfile, first time, append it
unless ( open(OUTSUMFILE5P,">>$outsumfname5p") ) {
	print "Can't create $outsumfname5p\n\n";
}

#open outfile, first time, append it
unless ( open(OUTFILE3P,">>$outfname3p") ) {
	print "Can't create $outfname3p\n\n";
}

#open outfile, first time, append it
unless ( open(OUTSUMFILE3P,">>$outsumfname3p") ) {
	print "Can't create $outsumfname3p\n\n";
}

#(0) Initialize the hash for storing the read count of each position
#------------------------------------------------------------------------------
for ($j = $minreadlength; $j <= $maxreadlength; $j++) {
	for ($i = 0; $i < $totalpos ; $i++) {
		$readcount5p{$j}{$i} = 0;
		$readcount3p{$j}{$i} = 0;
	}
}


#print header for the output
my $timestamp = localtime();
print OUTFILE5P "==================================================================================\n";
print OUTFILE5P "The Count5and3endAroundStart_v1.pl was run at $timestamp \n";
print OUTFILE5P "Parameters: left,right pos = -$lpos,$rpos; min, max read length = $minreadlength,$maxreadlength\n";
print OUTFILE5P "GENE\tREAD LENGTH(nt)\\POSITIONS\t";
for ($i = 0; $i < $totalpos; $i++) {
	my $pval = $i-$lpos;
	print OUTFILE5P "$pval\t";
}
print OUTFILE5P "\n";

print OUTFILE3P "==================================================================================\n";
print OUTFILE3P "The Count5and3endAroundStart_v1.pl was run at $timestamp \n";
print OUTFILE3P "Parameters: left,right pos = -$lpos,$rpos; min, max read length = $minreadlength,$maxreadlength\n";
print OUTFILE3P "GENE\tREAD LENGTH(nt)\\POSITIONS\t";
for ($i = 0; $i < $totalpos; $i++) {
	my $pval = $i-$lpos;
	print OUTFILE3P "$pval\t";
}
print OUTFILE3P "\n";


#(1) Calculate the minimum size of CDS based on the input left & right positions
#------------------------------------------------------------------------------
$minlen = $rpos+18; #For detail, see the OBJECTIVE section above


#(2) Go through every mRNA in the list, keep # read count in @readcount
#       this should include all the transcript variants!!?!?!
#------------------------------------------------------------------------------

#CDS list file is already opened

while ($line = <CDSFILE>) { #read line, each line = 1 mRNA

#chromosome	gene	strand	attributes #CDS	CDS1start,stop;CDS2start,stop
#9	mRNA	+	ID=GRMZM2G581216;Name=GRMZM2G581216;biotype=transposable_element	1	19970,20092
#file format --Version 3
#{chromosome}	{gene}	{attributes}	{#CDS}	{CDS1start,stop;CDS2start,stop}	{strand info}

	chomp ($line); #remove end-of-line char
	@temp = split (/\t/, $line);
	
	#calculate the length of CDS#################
	my $numCDS = $temp[3]; #number of CDS
	my @CDSinfo = split (/\;/, $temp[4]); #coordinate of the CDS
	my @strandinfo = split (/\;/, $temp[5]); #strand info of the CDS
	
	my $CDSlen = 0;
	for (my $i = 0; $i < $numCDS; $i++) {
		my @pos = split(/\,/, $CDSinfo[$i]);
		$CDSlen = $CDSlen+$pos[1]-$pos[0]+1; ##CAUTION!! CDS stop position always > start
	}
	
	#only continue if the length is at least $minlen
	if ($CDSlen >= $minlen) {

#9	mRNA	+	ID=GRMZM2G163722_T01;Parent=GRMZM2G163722;Name=GRMZM2G163722_T01;biotype=protein_coding	9	23314,23495;23940,24060;24284,24471;24541,24754;24842,24919;25127,25294;25364,25399;25845,25934;26231,26371
#9	mRNA	-	ID=GRMZM2G354611_T01;Parent=GRMZM2G354611;Name=GRMZM2G354611_T01;biotype=protein_coding	5	68562,68582;67887,68432;67067,67141;66607,66670;66347,66534
		
		#keep a list of need to check position
		my @plist = (); #keep the genome coordinate of the postion
		my %pcount5p = (); #keep the read count at that position. key = read length
		my %pcount3p = (); #keep the read count at that position. key = read length

		#check strand
		if ($strandinfo[0] eq '+') { #forward. unfortunately use the first CDS strand info only
		
			#make the list of position to check
			
			my @firstCDS = split(/\,/, $CDSinfo[0]); #get first position of CDS
			#write out the negative side position till before the first CDS position
			for (my $i = ($firstCDS[0]-$lpos); $i < $firstCDS[0]; $i++) { #neg region
				push (@plist, $i);
				for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
					$pcount5p{$k}{scalar(@plist)-1}= 0;
					$pcount3p{$k}{scalar(@plist)-1}= 0;
				}
			}
			#Write out the positive side position, tricky from differ exons
			my $neednt = $rpos +1;
			my $added = 0; #count how many pushed position
			my $exon = 0;
			my @pospair = split (/\,/,$CDSinfo[$exon]);
			my $current = $pospair[0];
			while ($added < $neednt){
				if ($current <= $pospair[1]) { #within exon
					push (@plist, $current);
					for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
						$pcount5p{$k}{scalar(@plist)-1}= 0;
						$pcount3p{$k}{scalar(@plist)-1}= 0;
					}
					$added++; #count
					$current++; #move next
				} else { #outside exon
					$exon++; #move to next exon
					@pospair = split (/\,/,$CDSinfo[$exon]); #new pos pair
					$current = $pospair[0]; #reset current
					push (@plist, $current);
					for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
						$pcount5p{$k}{scalar(@plist)-1}= 0;
						$pcount3p{$k}{scalar(@plist)-1}= 0;
					}
					$added++; #count
					$current++; #move next
				}
			}
			
			#RUN SAMTOOLS to check each position 5' counts of specific read length

			
			my $samcoor = "+ $temp[0]\:$plist[0]-".$plist[$totalpos-1];
			if (!exists($USECOOR{$samcoor})) { #new coordinate -> do analysis
				$USECOOR{$samcoor} = 1;
			
				#get all read cover this whole region
				open (SAMOUT, "samtools view ".
              				  "$bamfname ".
              				  "$temp[0]\:$plist[0]-".
              				  $plist[$totalpos-1]." |")
           		or die "Can't run samtools: $!\n";
			
				#for every read in SAMOUT (each line)
				while ( $samline = <SAMOUT> ) {
					chomp ($samline); #remove endline char
					my @col = (); #clear array
					@col = split (/\t/, $samline);
					my $cnum = scalar (@col); #number of column, can be differ, good to track
				
					my $sign = "+";
					my $bitcheck = $col[1] & 16; #check strand, - only if contain flag = 16
					if ($bitcheck == 16) { $sign = "-"; }
				
				
					#re-calculate POS and READ length based on the CIGAR operation
					my ($newpos, $matchlen, $newlen) = CIGARMOVE($sign, $col[3], $col[5], $col[9] );
				
					#newpos = 5' pos of either - or + strand read
					my $pos5end = $newpos;
					my $pos3end =0;
					if ($sign eq "-") { 
						$pos3end = $col[3];
					} else {
						$pos3end = $newpos+$matchlen-1; #use match length here b/c it is depict position on reference genome
						#$pos3end = $newpos+$newlen-1; 
					}
					
					#check if this read have with read length within the defined size
					if (($newlen >= $minreadlength) && ($newlen <= $maxreadlength)) {
						#count the end based on the new POS
						for (my $ii= 0; $ii < scalar(@plist); $ii++) {
							if ($plist[$ii] == $pos5end) {
								$nummatch5p++;
								$pcount5p{$newlen}{$ii} = $pcount5p{$newlen}{$ii}+1; #count this gene
								$readcount5p{$newlen}{$ii} = $readcount5p{$newlen}{$ii}+1; #count genome
							}
							if ($plist[$ii] == $pos3end) {
								$nummatch3p++;
								$pcount3p{$newlen}{$ii} = $pcount3p{$newlen}{$ii}+1; #count this gene
								$readcount3p{$newlen}{$ii} = $readcount3p{$newlen}{$ii}+1; #count genome
							}
						}#end count
					} #end check read length
				
				} #end SAMOUT
				close (SAMOUT);
			
				#print the count
				for (my $l = $minreadlength; $l <= $maxreadlength; $l++) {
					print OUTFILE5P "$temp[1]\t$l\t";
					print OUTFILE3P "$temp[1]\t$l\t";
					for (my $cc= 0; $cc < $totalpos; $cc++) {
						print OUTFILE5P "\t$pcount5p{$l}{$cc}";
						print OUTFILE3P "\t$pcount3p{$l}{$cc}";
					}
					print OUTFILE5P "\n";
					print OUTFILE3P "\n";
				}
				
			} else { #repeat coordinate, found in %USECOOR 
				print OUTFILE5P "$temp[1]\tREPEAT COORDINATE--SKIP\n";
				print OUTFILE3P "$temp[1]\tREPEAT COORDINATE--SKIP\n";
			}
			
		} else { #$temp[2] eq '-' # reverse strand

			#make the list of position to check
			
			my @firstCDS = split(/\,/, $CDSinfo[0]); #get first position of CDS
			#write out the negative side position till before the first CDS position
			for (my $i = ($firstCDS[1]+$lpos); $i > $firstCDS[1]; $i--) { #neg region
				push (@plist, $i);
				for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
					$pcount5p{$k}{scalar(@plist)-1}= 0;
					$pcount3p{$k}{scalar(@plist)-1}= 0;
				}
			}
			#Write out the positive side position, tricky from differ exons
			my $neednt = $rpos +1;
			my $added = 0; #count how many pushed position
			my $exon = 0;
			my @pospair = split (/\,/,$CDSinfo[$exon]);
			my $current = $pospair[1];
			while ($added < $neednt){
				if ($current >= $pospair[0]) { #within exon
					push (@plist, $current);
					for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
						$pcount5p{$k}{scalar(@plist)-1}= 0;
						$pcount3p{$k}{scalar(@plist)-1}= 0;
					}
					$added++; #count
					$current--; #move next
				} else { #outside exon
					$exon++; #move to next exon
					@pospair = split (/\,/,$CDSinfo[$exon]); #new pos pair
					$current = $pospair[1]; #reset current
					push (@plist, $current);
					for (my $k=$minreadlength; $k<=$maxreadlength; $k++){
						$pcount5p{$k}{scalar(@plist)-1}= 0;
						$pcount3p{$k}{scalar(@plist)-1}= 0;
					}
					$added++; #count
					$current--; #move next
				}
			}
			
			#RUN SAMTOOLS to check each position 5' counts of specific read length


			my $samcoor = "- $temp[0]\:".$plist[$totalpos-1]."-".$plist[0]; #####version 4##### 
			if (!exists($USECOOR{$samcoor})) { #new coordinate -> do analysis #####version 4##### 
				$USECOOR{$samcoor} = 1;

				#get all read cover this whole region
				open (SAMOUT, "samtools view ".
              				  "$bamfname ".
              				  "$temp[0]\:".$plist[$totalpos-1].
              				  "-$plist[0] |")
            	or die "Can't run samtools: $!\n";
			
				#for every read in SAMOUT (each line)
				while ( $samline = <SAMOUT> ) {
					chomp ($samline); #remove endline char
					my @col = (); #clear array
					@col = split (/\t/, $samline);
					my $cnum = scalar (@col); #number of column, can be differ, good to track
				
					my $sign = "+";
					my $bitcheck = $col[1] & 16; #check strand, - only if contain flag = 16
					if ($bitcheck == 16) { $sign = "-"; }
				
					#re-calculate POS and READ length based on the CIGAR operation
					my ($newpos, $matchlen, $newlen) = CIGARMOVE($sign,$col[3],$col[5],$col[9]);
				
					#newpos = 5' pos of either - or + strand read
					my $pos5end = $newpos;
					my $pos3end =0;
					if ($sign eq "-") { 
						$pos3end = $col[3];
					} else {
						$pos3end = $newpos+$matchlen-1; #use match length here b/c it is depict position on reference genome
						#$pos3end = $newpos+$newlen-1; 
					}
				
					#check if this read have with read length within the defined size
					if (($newlen >= $minreadlength) && ($newlen <= $maxreadlength)) {						#print "Found matched length read";
						#count the end based on the new POS
						for (my $ii= 0; $ii < scalar(@plist); $ii++) {
							if ($plist[$ii] == $pos5end) {
								$nummatch5p++;
								$pcount5p{$newlen}{$ii} = $pcount5p{$newlen}{$ii}+1; #count this gene
								$readcount5p{$newlen}{$ii} = $readcount5p{$newlen}{$ii}+1; #count genome
							}
							if ($plist[$ii] == $pos3end) {
								$nummatch3p++;
								$pcount3p{$newlen}{$ii} = $pcount3p{$newlen}{$ii}+1; #count this gene
								$readcount3p{$newlen}{$ii} = $readcount3p{$newlen}{$ii}+1; #count genome
							}
						}#end count
					} #end check read length
				
				} #end SAMOUT
				close (SAMOUT);
				
				#print the count
				for (my $l = $minreadlength; $l <= $maxreadlength; $l++) {
					print OUTFILE5P "$temp[1]\t$l\t";
					print OUTFILE3P "$temp[1]\t$l\t";
					for (my $cc= 0; $cc < $totalpos; $cc++) {
						print OUTFILE5P "\t$pcount5p{$l}{$cc}";
						print OUTFILE3P "\t$pcount3p{$l}{$cc}";
					}
					print OUTFILE5P "\n";
					print OUTFILE3P "\n";
				}
				
			} else { #repeat coordinate, found in %USECOOR 
				print OUTFILE5P "$temp[1]\tREPEAT COORDINATE--SKIP\n";
				print OUTFILE3P "$temp[1]\tREPEAT COORDINATE--SKIP\n";
			}

		}#end each strand as separate treatment
	
	}
	#if the CDS is too short, skip it
	$gcount++;
	#print "done with $gcount genes\n";
	
}# end while read line CDS file


#(3) Write output
#------------------------------------------------------------------------------
#print "\n\nFINAL COUNT: $gcount genes, $nummatch5p matched reads in 5'end, $nummatch3p matched reads in 3'end\n";
print OUTFILE5P "\n\n======================SUMMARY=============================\n";
print OUTFILE3P "\n\n======================SUMMARY=============================\n";
print OUTFILE5P "READ LENGTH(nt)\\POSITIONS";
print OUTFILE3P "READ LENGTH(nt)\\POSITIONS";
print OUTSUMFILE5P "READ LENGTH(nt)\\POSITIONS";
print OUTSUMFILE3P "READ LENGTH(nt)\\POSITIONS";
for ($i = 0; $i <= $lpos+$rpos; $i++) {
	my $pval = $i-$lpos;
	print OUTSUMFILE5P "\t$pval";
	print OUTSUMFILE3P "\t$pval";
}
print OUTSUMFILE5P "\n";
print OUTSUMFILE3P "\n";

for (my $l = $minreadlength; $l <= $maxreadlength; $l++) {
	print OUTFILE5P "$l";
	print OUTFILE3P "$l";
	print OUTSUMFILE5P "$l";
	print OUTSUMFILE3P "$l";
	
	for (my $c= 0; $c < $totalpos; $c++) {
		print OUTFILE5P "\t$readcount5p{$l}{$c}";
		print OUTFILE3P "\t$readcount3p{$l}{$c}";
		print OUTSUMFILE5P "\t$readcount5p{$l}{$c}";
		print OUTSUMFILE3P "\t$readcount3p{$l}{$c}";
	}
	
	print OUTFILE5P "\n";
	print OUTFILE3P "\n";
	print OUTSUMFILE5P "\n";
	print OUTSUMFILE3P "\n";
}

close (CDSFILE);
close (OUTFILE5P);
close (OUTFILE3P);
close (OUTSUMFILE5P);
close (OUTSUMFILE3P);

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

