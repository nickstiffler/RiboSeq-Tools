#!/usr/bin/env python3
#
# Create and submit batch script for
# running ribo-seq data analysis
#
# Nicholas Stiffler
# Barkan Lab
# February 15, 2018
#

from tkinter import *
from tkinter.filedialog import askopenfilename
import datetime, os, subprocess
from os import path

class SubmitRiboseq(Frame):

	def __init__(self, master=None):
		super().__init__(master)
		self.grid()
		self.create_widgets()

	def create_widgets(self):
	
		self.fileName = StringVar()
		self.fileName.set("Select FASTQ")
		#self.file_label = Label(self, textvariable=self.fileName).grid(row=0, column=0)

		self.file_button = Button(self, textvariable=self.fileName, command=self.get_file)
		self.file_button.grid(row=0, column=0)

		self.genome = StringVar()
		self.genome.set("Zea mays")
		self.genome_select = OptionMenu(self, self.genome, "Zea mays", "Arabidopsis thaliana")
		self.genome_select.grid(row=0, column=1)

		Label(self, text="Email address").grid(row=1, column=0)
		self.email = StringVar()
		try:
			self.email.set(os.getlogin() + "@uoregon.edu")
		except FileNotFoundError:
			print("Can't find the username.")

		Entry(self, textvariable=self.email).grid(row=1, column=1)

#		Label(self, text="DESeq2 Differential Expression").grid(row=2, column=2)
#		Label(self, text="Xtail Translational Efficiency").grid(row=2, column=3)

		Label(self, text="Trim adapter").grid(row=2, column=0)

		self.cutAdapt = StringVar()
		self.cutAdapt.set("Ribo-seq")
		self.cutAdapt_select = OptionMenu(self, self.cutAdapt, "Ribo-seq", "RNA-seq", "RIP-seq")
		self.cutAdapt_select.grid(row=2, column=1)

		Label(self, text="Alignment error rate").grid(row=3, column=0)
		self.error_rate = StringVar()
		self.error_rate.set("0.1")
		Spinbox(self, from_=0.01, to=0.3, increment=0.01, textvariable=self.error_rate).grid(row=3, column=1)
#		Checkbutton(self, text="Trim adapter small RNA kit", variable=self.cutAdapt).grid(row=3, column=0, sticky=W)

#		self.cutAdaptRNA = IntVar()
#		Checkbutton(self, text="Trim adapter RNA-seq kit", variable=self.cutAdaptRNA).grid(row=3, column=1, sticky=W)

		Label(self, text="Star Alignments").grid(row=4, column=0, sticky=W)
		Label(self, text="Feature Counts").grid(row=4, column=1, sticky=W)

		self.alignFilter = IntVar()
		self.alignFilter.set(1)
		Checkbutton(self, text="tRNA/rRNA filter", variable=self.alignFilter).grid(row=5, column=0, sticky=W)
		
		self.alignCP = IntVar()
		self.alignCP.set(1)
		Checkbutton(self, text="Align CP", variable=self.alignCP).grid(row=6, column=0, sticky=W)
		self.countsCP = IntVar()
		Checkbutton(self, text="CP Counts", variable=self.countsCP).grid(row=6, column=1, sticky=W)
		self.deseqCP = IntVar()
#		Checkbutton(self, text="CP DESeq2", variable=self.deseqCP).grid(row=4, column=2, sticky=W)
		self.xtailCP = IntVar()
#		Checkbutton(self, text="CP Xtail", variable=self.xtailCP).grid(row=4, column=3, sticky=W)
		
		self.alignMito = IntVar()
		self.alignMito.set(1)
		Checkbutton(self, text="Align Mito", variable=self.alignMito).grid(row=7, column=0, sticky=W)
		self.countsMito = IntVar()
#		Checkbutton(self, text="Mito Counts", variable=self.countsMito).grid(row=7, column=1, sticky=W)
		self.deseqMito = IntVar()
#		Checkbutton(self, text="Mito DESeq2", variable=self.deseqMito).grid(row=5, column=2)
		self.xtailMito = IntVar()
#		Checkbutton(self, text="Mito Xtail", variable=self.xtailMito).grid(row=5, column=3)
	
		self.alignNuc = IntVar()
		self.alignNuc.set(1)
		Checkbutton(self, text="Align Nuc", variable=self.alignNuc).grid(row=8, column=0, sticky=W)
		self.countsNuc = IntVar()
		self.countsNuc.set(1)
		Checkbutton(self, text="Nuc Counts", variable=self.countsNuc).grid(row=8, column=1, sticky=W)
		self.deseqNuc = IntVar()
#		Checkbutton(self, text="Nuc DESeq2", variable=self.deseqNuc).grid(row=6, column=2)
		self.xtailNuc = IntVar()
#		self.xtailNuc.set(1)
#		Checkbutton(self, text="Nuc Xtail", variable=self.xtailNuc).grid(row=6, column=3)
	
		Button(self, text="Create job", command=self.write_script).grid(row=9, column=0)
		Button(self, text="Submit job", command=self.submit_job).grid(row=9, column=1)

		self.console = Text()
		self.console.grid(row=10, columnspan=3)

#		self.quit = Button(self, text="QUIT", fg="red", command=root.destroy)
#		self.quit.pack(side="bottom")

	def get_file(self):
		self.fastq = askopenfilename()
		self.fileName.set(path.basename(self.fastq))

	def write_script(self):
		srun_content = """#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --time=0-10:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=""" + self.email.get()

		srun_content += """

module load python3
module load samtools
module load subread

/projects/barkanlab/shared/bin/riboseq_pipeline.py"""
		if self.cutAdapt.get() == "Ribo-seq":
			srun_content += " --trimRPF"
		if self.cutAdapt.get() == "RNA-seq":
			srun_content += " --trimRNA"
		if self.cutAdapt.get() == "RIP-seq":
			srun_content += " --trimRIP"
		if self.alignFilter.get():
			srun_content += " --filter"
		if self.alignCP.get():
			srun_content += " --alignCP"
		if self.countsCP.get():
			srun_content += " --countsCP"
		if self.deseqCP.get():
			srun_content += " --deseqCP"
		if self.xtailCP.get():
			srun_content += " --xtailCP"
		if self.alignMito.get():
			srun_content += " --alignMT"
		if self.countsMito.get():
			srun_content += " --countsMT"
		if self.deseqMito.get():
			srun_content += " --deseqMito"
		if self.xtailMito.get():
			srun_content += " --xtailMito"
		if self.alignMito.get():
			srun_content += " --alignNuc"
		if self.countsNuc.get():
			srun_content += " --countsNuc"
		if self.deseqNuc.get():
			srun_content += " --deseqNuc"
		if self.xtailNuc.get():
			srun_content += " --xtailNuc"
		if self.genome.get() == "Arabidopsis thaliana":
			srun_content += " --at"
		srun_content += " --errorRate " + self.error_rate.get()
		srun_content += " " + self.fastq + "\n"
		self.console.delete(1.0, END)	
		self.console.insert(INSERT, srun_content)

	def submit_job(self):
		srun_script = path.basename(self.fastq) + "_" + str(datetime.datetime.now())
		srun_script = srun_script.replace(' ', '_')
		srun_script = srun_script.replace(':', '-')
		os.mkdir(path.join(os.getcwd(), srun_script))
		os.chdir(path.join(os.getcwd(), srun_script))
		srun = open(srun_script, "w")
		srun.write(self.console.get(1.0, END))
		subprocess.Popen("sbatch " + path.join(os.getcwd(),srun_script), shell=True)
		os.chdir("..")
		self.console.insert(INSERT, "Batch script " + srun_script + " has been submitted")

root = Tk()
root.title("Maize Ribo-seq Analysis")
app = SubmitRiboseq(master=root)
app.mainloop()
