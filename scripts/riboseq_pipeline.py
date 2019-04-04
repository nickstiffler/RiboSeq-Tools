#!/usr/bin/env python3

import subprocess, os, argparse, re

reference_path = "/projects/barkanlab/shared/genomes/"
gtf_path = "/projects/barkanlab/shared/gff/"
run_folder = os.getcwd()
star_cmd = "/projects/barkanlab/shared/bin/STAR"
zm = ("Zmays_tRNA_rRNA.star", "Zea_mays.Pt.star", "Zea_mays.Mt.star", "Zea_mays.AGPv4.star2")
at = ("Arabidopsis.tRNA_rRNA.star", "At_TAIR10_cp_truncate", "Arabidopsis.Mt.star", "Arabidopsis.TAIR10.star")
qc_path = "/projects/barkanlab/shared/QCfiles"

threads = 28

parser = argparse.ArgumentParser()
parser.add_argument("--trimRPF", action="store_true")
parser.add_argument("--trimRNA", action="store_true")
parser.add_argument("--trimRIP", action="store_true")
parser.add_argument("--filter", action="store_true")
parser.add_argument("--alignCP", action="store_true")
parser.add_argument("--countsCP", action="store_true")
parser.add_argument("--alignMT", action="store_true")
parser.add_argument("--countsMT", action="store_true")
parser.add_argument("--alignNuc", action="store_true")
parser.add_argument("--countsNuc", action="store_true")
parser.add_argument("--at", action="store_true")
parser.add_argument("--errorRate", default=0.1)
parser.add_argument("fastq")
args = parser.parse_args()

genome_files = zm
if args.at:
	genome_files = at

gtf_prefix = "maize"
if args.at:
	gtf_prefix = "arab"

def pipeline():
	# Remove extension from fastq
	if "gz" in args.fastq:
		current = os.path.splitext(os.path.basename(os.path.splitext(os.path.basename(args.fastq))[0]))[0]
	else:
		current = os.path.splitext(os.path.basename(args.fastq))[0]
	if args.trimRPF:
		current = runCutAdaptRPF(args.fastq, current)
	if args.trimRNA:
		current = runCutAdaptRNA(args.fastq, current)
	if args.trimRIP:
		current = runCutAdaptRIP(args.fastq, current)
	if args.filter:
		current = runFilter(current, genome_files[0])
	if args.alignCP:
		current = runSTAR(current, genome_files[1], "--alignIntronMax 2500 --alignEndsType EndToEnd --outFilterMismatchNoverLmax " + args.errorRate)
		if args.trimRPF:
			runQC(current[1], "cp")
	if args.countsCP:
		if args.trimRNA:
			counts = runFeatureCountsRNA(current[1], gtf_prefix + "_pt_rnaseq_ripseq.gff")
		elif args.trimRIP:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_pt_ripseq.gff")
			calcRPKM(counts[0], counts[1])
		else:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_pt_riboseq.gff")

	if args.alignMT:
		current = runSTAR(current[0], genome_files[2], "--alignEndsType EndToEnd")
		if args.trimRPF:
			runQC(current[1], "mt")
	if args.countsMT:
		if args.trimRNA:
			counts = runFeatureCountsRNA(current[1], gtf_prefix + "_mt_rnaseq_ripseq.gff")
		elif args.trimRIP:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_mt_rnaseq_ripseq.gff")
			calcRPKM(counts[0], counts[1])
		else:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_mt_riboseq.gff")

	if args.alignNuc:
		current = runSTAR(current[0], genome_files[3], "--alignEndsType EndToEnd --outFilterMismatchNoverLmax " + args.errorRate)
		if args.trimRPF:
			runQC(current[1], "nuc")
	if args.countsNuc:
		if args.trimRNA:
			counts = runFeatureCountsRNA(current[1], gtf_prefix + "_rnaseq_ripseq.gff")
		elif args.trimRIP:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_rnaseq_ripseq.gff")
			calcRPKM(counts[0], counts[1])
		else:
			counts = runFeatureCountsRPF(current[1], gtf_prefix + "_riboseq.gff")

		
def runCutAdaptRPF(filename, sample):
	output = os.path.join(run_folder, sample + ".fq")
	command = "{} {} {} {} {} {} {}".format("cutadapt -j", threads, "-a NNNNTGGAATTCTCGGGTGCCAAGG", "-u 4", "-o " + output, "-m 18 -M 40", filename)
	print(command)
	subprocess.call(command, shell=True)
	return os.path.splitext(output)[0]

def runCutAdaptRIP(filename, sample):
	output = os.path.join(run_folder, sample + ".fq")
	command = "{} {} {} {} {} {} {}".format("cutadapt -j", threads, "-a NNNNTGGAATTCTCGGGTGCCAAGG", "-u 4", "-o " + output, "-m 2", filename)
	print(command)
	subprocess.call(command, shell=True)
	return os.path.splitext(output)[0]

def runCutAdaptRNA(filename, sample):
	output = os.path.join(run_folder, sample + ".fq")
	command = "{} {} {} {} {} {} {}".format("cutadapt -j", threads, "-a NNNNNNNNNTGGAATTCTCGGGTGCCAAGG", "-u 9", "-o " + output, "-m 2", filename)
	print(command)
	subprocess.call(command, shell=True)
	return os.path.splitext(output)[0]

def runSamtools(sample):
	command = "{} {} {} {} {} {}{}".format("samtools", "view --threads ", threads, "-b -S", os.path.join(run_folder, sample) + " >", os.path.join(run_folder, sample), "_unsorted.bam")
	subprocess.call(command, shell=True)
	
	command = "{} {} {} {}{} {}".format("samtools", "sort --threads ", threads, os.path.join(run_folder, sample), "_unsorted.bam", " -o " + os.path.join(run_folder, sample) + ".bam")
	print(command)
	subprocess.call(command, shell=True)
	
	command = "{} {} {}{}".format("samtools", "index", os.path.join(run_folder, sample), ".bam")
	subprocess.call(command, shell=True)

def runPileup(sample, ref):
	command = "{} {} {} {}{} {}{}".format("samtools", "mpileup -d 10000000 -f ", ref, sample, ".bam | cut -f 1-4 > ", sample, "_rRNA_coverage.xls")
	subprocess.call(command, shell=True)

def runSTAR(sample, ref, params=""):
	command = "{} {} {} {} {} {} {} {} {} {}".format(star_cmd, "--runThreadN", threads, "--runMode alignReads", "--genomeDir", os.path.join(reference_path,ref), "--readFilesIn", sample + ".fq", "--outReadsUnmapped Fastx", params)
	print(command)
	subprocess.call(command, shell=True)
	sam_output = sample + "_" + ref
	os.rename("Aligned.out.sam", sam_output)
	runSamtools(sam_output)	
	output = sample + "_no_" + ref + ".fq"
	os.rename("Unmapped.out.mate1", output)
	return (os.path.splitext(output)[0], sam_output + ".bam")

def runFilter(sample, ref):
	command = "{} {} {} {} {} {} {} {} {}".format(star_cmd, "--runThreadN", threads, "--runMode alignReads", "--genomeDir", os.path.join(reference_path,ref), "--readFilesIn", sample + ".fq", "--outReadsUnmapped Fastx")
	print(command)
	subprocess.call(command, shell=True)
	output = sample + "_no_" + ref + ".fq"
	os.rename("Unmapped.out.mate1", output)
	return os.path.splitext(output)[0]

def runFeatureCountsRNA(sample, gtf):
	command = "{} {} {} {} {} {} {} {} {}".format("featureCounts", "-O -T", threads, "-s 2 --primary --largestOverlap -t CDS -g gene_id", "-a", os.path.join(gtf_path, gtf), "-o", sample + ".counts", sample)

	print(command)
	counts_out = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	print(counts_out.stderr.decode("utf-8"))
	total_reads = re.search(r'Successfully assigned reads : (\d+)', counts_out.stderr.decode("utf-8")).group(1)
	return (sample + ".counts", int(total_reads))


def runFeatureCountsRPF(sample, gtf):
	command = "{} {} {} {} {} {} {} {} {} {}".format("featureCounts", "-O -T", threads, "-s 1", "--primary --largestOverlap -t CDS -g gene_id", "-a", os.path.join(gtf_path, gtf), "-o", sample + ".counts", sample)

	print(command)
	counts_out = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	print(counts_out.stderr.decode("utf-8"))
	total_reads = re.search(r'Successfully assigned reads : (\d+)', counts_out.stderr.decode("utf-8")).group(1)
	return (sample + ".counts", int(total_reads))


def calcTPM(featureCount):
	rpks = []
	total_rpk = 0
	with open(featureCount) as fc_file:
		next(fc_file)
		next(fc_file)
		for line in fc_file:
			cols = line.strip().split()
			rpks.append(float(cols[6]) / float(cols[5]) / 1000.0)
			total_rpk += (float(cols[6]) / float(cols[5]) / 1000.0)
	scaling_factor = total_rpk / 1000000.0
	i = 0
	with open(featureCount) as fc_file:
		next(fc_file)
		next(fc_file)
		print("Gene name\tCounts\tTPM")
		for line in fc_file:
			cols = line.strip().split()
			print("{}\t{}\t{}".format(cols[0], cols[6], rpks[i] / scaling_factor))
			i += 1

def calcRPKM(featureCount, total_counts):
	counts_file = open(featureCount + ".rpkm", 'w')
#	total_counts = 0
#	with open(featureCount) as fc_file:
#		next(fc_file)
#		next(fc_file)
#		for line in fc_file:
#			cols = line.strip().split()
#			total_counts += float(cols[6])
	scaling_factor = total_counts / 1000000.0
	i = 0
	length_fixes = {"ycf1_5UTR": 100, "ycf1_coding": 5360, "ycf1_all": 5560}
	with open(featureCount) as fc_file:
		next(fc_file)
		next(fc_file)
		counts_file.write("Gene name\tLength\tCounts\tRPKM\n")
		for line in fc_file:
			cols = line.strip().split()
			if cols[0] in length_fixes:
				cols[5] = length_fixes[cols[0]]
			counts_file.write("{}\t{}\t{}\t{}\n".format(cols[0], cols[5], cols[6], (float(cols[6]) / scaling_factor / (float(cols[5]) / 1000.0))))
			i += 1

def runQC(bamFile, reference):
	organism = "Zm"
	if args.at:
		organism = "At"
	command = "{} {} {} {} {}".format("QC_STAR_Ribo-seq_v5.pl", qc_path, organism, reference, bamFile)

	print(command)
	subprocess.call(command, shell=True)

	

pipeline()

