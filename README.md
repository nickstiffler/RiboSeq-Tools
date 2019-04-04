# RiboSeq-Tools

Set of scripts for routine data analysis of plant transcriptomic Illumina sequences. Although originally written to handle Ribo-seq data, can be used for RNA-seq and Rip-seq. 

# Dependencies

These packages are required by the pipeline and need to be in the path of executing user.

- Python 3
- Cutadapt (https://cutadapt.readthedocs.io/en/stable/installation.html)
- STAR (https://github.com/alexdobin/STAR)
- Samtools (http://www.htslib.org/)
- FeatureCounts (http://subread.sourceforge.net/)

# Pipeline

The main data analysis is completed using the `riboseq_pipeline.py` script. Before execution, the script must be configured. Open the script a text editor and set the `reference_path` to the location of the reference genome files indexed using STAR. There needs to be a separate file for chloroplast, mitochondria and nuclear genomes. Then set the `gtf_path` to the location of the gff files. There also needs to be a separate gff file for each organelle genome. The `threads` parameter can be set to indicate the number of cores that can be used for multi-threaded applications.

`usage: riboseq_pipeline.py [-h] [--trimRPF] [--trimRNA] [--trimRIP] [--filter]
                           [--alignCP] [--countsCP] [--alignMT] [--countsMT]
                           [--alignNuc] [--countsNuc] [--at]
                           [--errorRate ERRORRATE]
                           fastq
`
Command-line parameters:

`--trimRPF` - Run the trimming step (cutadapt) and the data is Ribo-seq.

`--trimRNA` - Run the trimming step (cutadapt) and the data is RNA-seq.

`--trimRIP` - Run the trimming step (cutadapt) and the data is Rip-seq.

`--filter` - Run the filter step to remove contaminate by aligning reads to rRNA and tRNA reference using STAR.

`--alignCP` - Align reads to the chloroplast genome.

`--countCP` - Produce table of counts from the chloroplast alignment using `featureCounts`.

`--alignMT` - Align reads to the mitochondria genome.

`--countMT` - Produce table of counts from the mitochondria alignment using `featureCounts`.

`--alignNuc` - Align reads to the nuclear genome.

`--countNuc` - Produce table of counts from the nuclear alignment using `featureCounts`.

`--at` - Use Aribdopsis reference genome instead of maize

`--errorRate` - Specify error rate to use for STAR alignments. The default is 0.1.

# GUI

The script `submit_riboseq_analysis.py` is a graphical user interface designed to run using X11 forwarding over SSH so users and configure and submit the pipeline analysis on a compute cluster. It will produce a script to run the on the `SLURM` queueing system on the Talapas cluster at the University of Oregon.

# QC

An additional script is provided to run the quality control steps called `QC_STAR_Ribo-seq_v5.pl`. 

# Bootstrap

A script called `barkan_bootstrap.sh` is included to setup the user accounts on Talapas to make it easier for users to run the GUI. It updates the profile to include the location of these scripts in the user path, loads the module for python3, automatically put the user in the correct project folder on login, updates permissions to the user's project folder so othe r members of the lab can access their files and creates a symbolic link to the shared Illumina data folder into the user's project folder. It must be run only once.
