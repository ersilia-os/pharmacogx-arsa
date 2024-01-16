# Import libraries
import argparse
import os
import subprocess
import glob
import time
import sys


# Input
parser = argparse.ArgumentParser(prog='pipeline_1kGPhg38_GENCODEv43_intersect',
								 description='Pipeline to intersect 1kGPhg38 variants with GENCODEv43')


parser.add_argument('-i', '--inDir',
					dest = "inDir",
					action = "store",
					help = "Path to input directory containing vcf.gz files") 


parser.add_argument('-bed', '--bedFile',
					dest = "bedFile",
					required=True,
					action = "store",
					help = "Path and name of GENCODEv43 regions file in BED format" )

parser.add_argument('-w', '--workDir',
					dest = "workDir",
					required=True,
					action = "store",
					help = "Working directory")   

parser.add_argument('-t', '--time',
					dest = "time",
					action = "store",
					default = "72",
					help = "Maximum execution time (in hours). Default = 72.")


parser.add_argument('-@', '--cpus',
					dest = "cpus",
					action = "store",
					default = "4",
					choices=['1', '2', '4', '8', '16','32'],
					help = "Number of CPUs [1, 2, 4, 8, 16, 32]. Default = 4.") 


options = parser.parse_args()


# Starting time
date_str = time.strftime("%Y/%m/%d_%H/%M/%S").replace("/","")


# Prepare commands
listVcfs = [fq for fq in glob.glob(options.inDir+"/*.vcf.gz")]
if options.tasks != "":
	nTasks = options.tasks
else: 
	nTasks = str(len(listVcfs)) if len(listVcfs) <= 20 else "20"


# Define variables
inDir = options.inDir+"/" if not options.inDir.endswith("/") else options.inDir
target_files_folder = options.bedFile
workDir = options.workDir+"/" if not options.workDir.endswith("/") else options.workDir
intermediate_data = workDir+"intermediate_data/"
outDir = workDir+"intragenic_variants/"
os.makedirs(outDir)

currentDir = os.getcwd()


# Intersect files
chrList=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]

for chr in chrList: 

	# Run one job per chromosome
	numJobs=0

	O = open("job_script.sh", "w")
	O.write("#!/bin/bash\n")
	O.write("#SBATCH --job-name=varSubset\n")
	O.write("#SBATCH -D '.'\n")
	O.write("#SBATCH --output=%j.out\n")
	O.write("#SBATCH --error=%j.err\n")
	O.write("#SBATCH --ntasks=1\n") 
	O.write("#SBATCH --time=00:30:00\n")
	O.write("#SBATCH --cpus-per-task=8\n")
	O.write("\nmodule load bedtools/2.29.1 htslib/1.10.2 bcftools/1.10.2\n\n")


	bashArguments="bedtools intersect -u -a "+inDir+"ALL."+chr+".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz  -b "+target_files_folder+"adapted_knownGenes_gencodeV43_hg38.bed | cut -f1,2,3,4,5,6,7,8,9  > "+intermediate_data+"uniq_coding_"+chr+"_1kGPh38.vcf"
	O.write(bashArguments+"\n")

	bashArguments="cat <(cut -f1,2,3,4,5,6,7,8,9 <(bcftools view -h "+inDir+"ALL."+chr+".shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz)) "+intermediate_data+"uniq_coding_"+chr+"_1kGPh38.vcf > "+intermediate_data+"/reheaded_uniq_coding_"+chr+"_1kGPhg38.vcf"
	O.write(bashArguments+"\n")

	bashArguments="bgzip -c "+intermediate_data+"reheaded_uniq_coding_"+chr+"_1kGPhg38.vcf >  "+outDir+"reheaded_uniq_coding_"+chr+"_1kGPhg38.vcf.gz" 
	O.write(bashArguments+"\n")
	
	bashArguments="tabix "+outDir+"reheaded_uniq_coding_"+chr+"_1kGPhg38.vcf.gz"
	O.write(bashArguments+"\n")

	O.close()

	# Launch script

	subprocess.call("sbatch job_script.sh", shell=True)
	numJobs += 1






