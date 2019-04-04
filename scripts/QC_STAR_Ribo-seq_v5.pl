#!/usr/bin/perl
############################################################################
# Created by Prakitchai Chotewutmontri, 19 June 2018
# Postdoc, Barkan Lab, U of Oregon
#
# QC_STAR_Ribo-seq_v1.pl
#
# OBJECTIVES
#  Quality control of Ribo-seq results
#    1) Ribosome footprint length distribution (cp, mt, nuc)
#    2) Coverage of 5'end of RF around START codon (cp, mt, nuc)
#    3) Coverage of 3'end of RF around START codon (cp, mt, nuc)
#    4) Coverage of 5'end of RF around STOP codon (cp, mt, nuc)
#    5) Coverage of 3'end of RF around STOP codon (cp, mt, nuc)
#    6) Coverage of 3'end of RF for the whole genome (cp)
#    7) Read coverage in the CDS and Pearson correlation (only for MAIZE cp genome for now -- need standard coverage -- only made the maize standard...will make arabidopsis at some point)
#    8) 3-nt periodicity (cp, mt, nuc)
#
#Version 2 -- Add extra in-line argument for selection if it maize or Arabidopsis
#Version 3 -- generalize to take in 1 bam file at a time -- no longer
#Version 4 -- Nick's .bai file is in this format "file.bam.bai" but Non's script need "file.bai" format --fix this by copy .bam.bai to .bai
#Version 5 -- create a new .bai file in this 
############################################################################

############################################################################
#DECLARATIONS
############################################################################
use strict;
use warnings;

#global var
my $sysline = "";
my $timestring = "";

my $i = 0;

#usage
my $USAGE = "usage: $0 <QC folder path> <Zm or At> <cp/mt/nuc> <.bam file>\n\n";

#check argument
unless (@ARGV == 4) { 
	print $USAGE;
	exit -1;
}

#store argv
my $QCdirPath = $ARGV[0];
my $Species = $ARGV[1]; #have to be Zm or At
my $genomechoices = $ARGV[2]; #have to be Zm or At
my $bamfname = $ARGV[3];
#my $baifname = "$bamfname".".bai"; ##-- VERSION 4


#my $shortbai = substr($bamfname,0,length($bamfname)-4); ##-- VERSION 4
#$shortbai = "$shortbai".".bai"; ##-- VERSION 4
#copy Nick's .bai format filename to Non's filename -- VERSION 4
#system('cp', $baifname, $shortbai); ##-- VERSION 4

#version 5
my $shortbai = substr($bamfname,0,length($bamfname)-4); ##-- VERSION 5
$shortbai = "$shortbai".".bai"; ##-- VERSION 5
$sysline = "samtools index $bamfname $shortbai";  ##-- VERSION 5
system($sysline); ##-- VERSION 5

## Make selection for species specific annotation files
###PREDEFINED ANNOTATION FILES in the QC folder
my $cpAnnotation = "";
my $cpAnnotationEXTRA = "";
my $mtAnnotation = "";
my $nucAnnotation = "";
my $mtAnnotationALL = "";
my $nucAnnotationALL = "";
my $cp_pos_cover = "";
my $cp_neg_cover = "";

if ( $Species eq 'Zm' ) { #maize v4 pipeline
	###PREDEFINED ANNOTATION FILES in the QC folder
	$cpAnnotation = "$QCdirPath/Alice_082714_CDS_Pt_withNum_v9.trucate_117737.txt";
	$cpAnnotationEXTRA = "$QCdirPath/Alice_082714_CDS_Pt_withNum_v9_last-exon.trucate_117737.txt";
	$mtAnnotation = "$QCdirPath/Zea_mays.AGPv4.38.chromosome.Mt.UnionCDS.txt";
	$nucAnnotation = "$QCdirPath/Zea_mays.AGPv4.38.NUConly.UnionCDS.txt";
	$mtAnnotationALL = "$QCdirPath/Zea_mays.AGPv4.38.chromosome.Mt.All_transcript.txt";
	$nucAnnotationALL = "$QCdirPath/Zea_mays.AGPv4.38.NUConly.All_transcript.txt";
	$cp_pos_cover = "$QCdirPath/X86563_avg_coverage_per_million_forward.txt";
	$cp_neg_cover = "$QCdirPath/X86563_avg_coverage_per_million_reverse.txt";
} elsif ( $Species eq 'At' ) { #arabidopsis TAIR10
	###PREDEFINED ANNOTATION FILES in the QC folder
	$cpAnnotation = "$QCdirPath/TAIR10_AllCDS_ChrC_Non_v3_truncate129736.txt";
	$cpAnnotationEXTRA = "$QCdirPath/TAIR10_AllCDS_ChrC_plus_last_exon_Non_v3_truncate129736.txt";
	$mtAnnotation = "$QCdirPath/TAIR10_UnionCDS_ChrM.txt";
	$nucAnnotation = "$QCdirPath/TAIR10_UnionCDS_Chr1-5.txt";
	$mtAnnotationALL = "$QCdirPath/TAIR10_AllCDS_ChrM.txt";
	$nucAnnotationALL = "$QCdirPath/TAIR10_AllCDS_Chr1-5.txt";
	##$cp_pos_cover = "$QCdirPath/"; #don't have the standard coverage distribution for Arabidopsis yet 
	##$cp_neg_cover = "$QCdirPath/"; #don't have the standard coverage distribution for Arabidopsis yet 
} else { #not correctly specify "Zm" or "At"
	print "Specify species choice incorrectly! Use Zm or At only";
	exit -1;	
}

if ( $genomechoices ne 'cp' ) {
	if ( $genomechoices ne 'mt' ) {
		if ( $genomechoices ne 'nuc' ) {
			print "Specify genome choice incorrectly! Use cp, mt or nuc only";
			exit -1;	
		}
	}
}

##Perl script files
my $perl_RFlength = "$QCdirPath/CountReadLengthInCDS_v3.pl";
my $perl_start = "$QCdirPath/Count5and3endAroundStart_v1.pl";
my $perl_stop = "$QCdirPath/Count5and3endAroundStop_v1.pl";
my $perl_3pcover = "$QCdirPath/Count3primeCoverage_v1.pl";
my $perl_Pearson = "$QCdirPath/CDS_coverage_and_correlation_v2.pl";
my $perl_3ntcp = "$QCdirPath/CalPeriodicity_v4_cp.pl";
my $perl_3ntmt = "$QCdirPath/CalPeriodicity_v4_mt_26-30nt.pl";
my $perl_3ntnuc = "$QCdirPath/CalPeriodicity_v5_nuc_28-33nt.pl";

##output files
my $cpRFoutfname = "$bamfname.QC_RF_length.txt";
my $mtRFoutfname = "$bamfname.QC_RF_length.txt";
my $nucRFoutfname = "$bamfname.QC_RF_length.txt";

my $cpout5endstart = "$bamfname.QC_5pSTART_detail.txt";
my $cpout5endstartsum = "$bamfname.QC_5pSTART_sum.txt";
my $cpout3endstart = "$bamfname.QC_3pSTART_detail.txt";
my $cpout3endstartsum = "$bamfname.QC_3pSTART_sum.txt";

my $mtout5endstart = "$bamfname.QC_5pSTART_detail.txt";
my $mtout5endstartsum = "$bamfname.QC_5pSTART_sum.txt";
my $mtout3endstart = "$bamfname.QC_3pSTART_detail.txt";
my $mtout3endstartsum = "$bamfname.QC_3pSTART_sum.txt";

my $nucout5endstart = "$bamfname.QC_5pSTART_detail.txt";
my $nucout5endstartsum = "$bamfname.QC_5pSTART_sum.txt";
my $nucout3endstart = "$bamfname.QC_3pSTART_detail.txt";
my $nucout3endstartsum = "$bamfname.QC_3pSTART_sum.txt";

my $cpout5endstop = "$bamfname.QC_5pSTOP_detail.txt";
my $cpout5endstopsum = "$bamfname.QC_5pSTOP_sum.txt";
my $cpout3endstop = "$bamfname.QC_3pSTOP_detail.txt";
my $cpout3endstopsum = "$bamfname.QC_3pSTOP_sum.txt";

my $mtout5endstop = "$bamfname.QC_5pSTOP_detail.txt";
my $mtout5endstopsum = "$bamfname.QC_5pSTOP_sum.txt";
my $mtout3endstop = "$bamfname.QC_3pSTOP_detail.txt";
my $mtout3endstopsum = "$bamfname.QC_3pSTOP_sum.txt";

my $nucout5endstop = "$bamfname.QC_5pSTOP_detail.txt";
my $nucout5endstopsum = "$bamfname.QC_5pSTOP_sum.txt";
my $nucout3endstop = "$bamfname.QC_3pSTOP_detail.txt";
my $nucout3endstopsum = "$bamfname.QC_3pSTOP_sum.txt";

my $cpout3endcover = "$bamfname.QC_3pCOVER.txt";
my $cpout3endcoverwig = "$bamfname.QC_3pCOVER";

my $cp3ntperi = "$bamfname.QC_3ntPeriodicity";
my $mt3ntperi = "$bamfname.QC_3ntPeriodicity";
my $nuc3ntperi = "$bamfname.QC_3ntPeriodicity";


############################################################################
#MAIN PROGRAM
############################################################################
#    1) RF length
#    usage: CountReadLengthInCDS.pl <input .BAM file> <input CDS list> <output filename>
# CP
if ( $genomechoices eq 'cp' ) {
	$sysline = "perl $perl_RFlength $bamfname $cpAnnotation $cpRFoutfname";
	system ("$sysline");
}
# MT
if ( $genomechoices eq 'mt' ) {
	$sysline = "perl $perl_RFlength $bamfname $mtAnnotation $mtRFoutfname";
	system ("$sysline");
}
# NUC
if ( $genomechoices eq 'nuc' ) {
	$sysline = "perl $perl_RFlength $bamfname $nucAnnotation $nucRFoutfname";
	system ("$sysline");
}

############################################################################
#    2) to 3) Coverage of 5'end and 3' End of RF around START codon
#    usage: Count5and3endAroundStart_v1.pl <input .BAM file> <input CDS list> <5' outfile> <5' out sum file> <3' outfile> <3' out sum file> <negative pos, no sign> <positive pos, from zero> <min read length> <max read length>
#RF length of 16 to 40
# CP
if ( $genomechoices eq 'cp' ) {
	$sysline = "perl $perl_start $bamfname $cpAnnotation $cpout5endstart $cpout5endstartsum $cpout3endstart $cpout3endstartsum 24 62 16 40";
	system ("$sysline");
}
# MT
if ( $genomechoices eq 'mt' ) {
	$sysline = "perl $perl_start $bamfname $mtAnnotationALL $mtout5endstart $mtout5endstartsum $mtout3endstart $mtout3endstartsum 24 62 16 40";
	system ("$sysline");
}
# NUC
if ( $genomechoices eq 'nuc' ) {
	$sysline = "perl $perl_start $bamfname $nucAnnotationALL $nucout5endstart $nucout5endstartsum $nucout3endstart $nucout3endstartsum 24 62 16 40";
	system ("$sysline");
}


############################################################################
#    4) to 5) Coverage of 5'end and 3' End of RF around STOP codon
#    usage: Count5and3endAroundStop_v1.pl <input .BAM file> <input CDS list> <5' outfile> <5' out sum file> <3' outfile> <3' out sum file> <negative pos, no sign> <positive pos, from zero> <min read length> <max read length>
#RF length of 16 to 40
# CP
if ( $genomechoices eq 'cp' ) {
	$sysline = "perl $perl_stop $bamfname $cpAnnotation $cpout5endstop $cpout5endstopsum $cpout3endstop $cpout3endstopsum 48 50 16 40";
	system ("$sysline");
}
# MT
if ( $genomechoices eq 'mt' ) {
	$sysline = "perl $perl_stop $bamfname $mtAnnotationALL $mtout5endstop $mtout5endstopsum $mtout3endstop $mtout3endstopsum 48 50 16 40";
	system ("$sysline");
}
# NUC
if ( $genomechoices eq 'nuc' ) {
	$sysline = "perl $perl_stop $bamfname $nucAnnotationALL $nucout5endstop $nucout5endstopsum $nucout3endstop $nucout3endstopsum 48 50 16 40";
	system ("$sysline");
}


############################################################################
#    6) Coverage of 3'end of RF for the whole cp genome
if ( $genomechoices eq 'cp' ) {
	$sysline = "perl $perl_3pcover $bamfname $cpout3endcover $cpout3endcoverwig";
	system ("$sysline");
}

############################################################################
#    7) CDS read coverage and Pearson correlation
## only run this for Maize for now!! - 6/19/2018 (Non)
if ( $genomechoices eq 'cp' ) {
	open (BAMIDXSTATS, "samtools idxstats $bamfname |")
	or die "Can't run samtools: $!\n";
	#Mt	569630	1085155	0
	#*	0	0	13555311
	#
	#X86563	140384	3477745	0
	#*	0	0	16121217
	my $idxstats = <BAMIDXSTATS>;
	my @tmp = split (/\t/, $idxstats);
	my $chr = $tmp[0]; ## maize chloroplast genome that Nick use is "Pt"
	my $clength = $tmp[1]; ##maize truncated cp genome = 117737 bp
	close (BAMIDXSTATS);

	if ($chr eq 'Pt' && $clength == 117737) { ##only run the maize genome 
		$sysline = "perl $perl_Pearson $bamfname $cpAnnotationEXTRA $cp_pos_cover $cp_neg_cover";
		system ("$sysline");
	}
}

############################################################################
#    8) 3nt-periodicity
# CP
if ( $genomechoices eq 'cp' ) {
	$sysline = "perl $perl_3ntcp $bamfname $cpAnnotation $cp3ntperi";
	system ("$sysline");
}
# MT
if ( $genomechoices eq 'mt' ) {
	$sysline = "perl $perl_3ntmt $bamfname $mtAnnotationALL $mt3ntperi";
	system ("$sysline");
}
# NUC
if ( $genomechoices eq 'nuc' ) {
	$sysline = "perl $perl_3ntnuc $bamfname $nucAnnotationALL $nuc3ntperi";
	system ("$sysline");
}


############################################################################
#    9) remove the detail file to save space
$sysline = "rm *.QC*detail.*";
system ("$sysline");

############################################################################


#VERSION 4 remove Non's format bai file
#system('rm', $shortbai); ##-- VERSION 4
            
exit 0;