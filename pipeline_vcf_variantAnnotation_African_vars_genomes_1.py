# Libraries:
import subprocess
import argparse
import time
import os
import glob


# Inputs:
parser = argparse.ArgumentParser(prog='pipeline_vcf_variantAnnotation_African_vars_genomes_1',
                                 description='''Pipeline to run variant annotation''')

parser.add_argument('-s', '--pathToReferencePipelines',
                    dest="pathToReferencePipelines",
                    action="store",
                    default="/slgpfs/projects/cli79/reference_pipelines",
                    help="Path to folder reference_pipelines")

parser.add_argument('-p', '--pathToReferenceData',
					dest = "pathToReferenceData",
					action = "store",
					default = "/slgpfs/projects/cli79/reference_data",
					help = "Path to folder reference_data")

parser.add_argument('-pr', '--project',
                    dest="project",
                    action="store",
                    choices= ['1kGPhg38', 'PharmGKB', 'Fedorova'],
                    help="Project name")

parser.add_argument('-i', '--inDir',
                    dest="inDir",
                    action="store",
                    help="Path to input directory containing vcf files")

parser.add_argument('-w', '--workDir',
                    dest="workDir",
                    action="store",
                    help="Working directory")

parser.add_argument('-tBed', '--targetBedFile',
                    dest="targetBedFile",
                    required=True,
                    action="store",
                    help="Target bed file")

parser.add_argument('-trans', '--transcriptsToUse',
                    dest="transcriptsToUse",
                    default="canonical",
                    action="store",
                    help="Txt file with transcripts to use for annotations or use 'canonical' to annotate canonical files")

parser.add_argument('-intron', '--intronMutations',
                    dest="intronMutations",
                    default="yes",
                    action="store",
                    choices=['yes', 'no'],
                    help="Should intron mutations be annotated? Defaults is 'yes'")

parser.add_argument('-intergenic', '--intergenicMutations',
                    dest="intergenicMutations",
                    default="no",
                    action="store",
                    choices=['yes', 'no'],
                    help="Should intergenic mutations be annotated? Defaults is 'no'")

parser.add_argument('-downstream', '--downstreamMutations',
                    dest="downstreamMutations",
                    default="no",
                    action="store",
                    choices=['yes', 'no'],
                    help="Should downstream mutations be annotated? Defaults is 'no'")

parser.add_argument('-upstream', '--upstreamMutations',
                    dest="upstreamMutations",
                    default="no",
                    action="store",
                    choices=['yes', 'no'],
                    help="Should upstream mutations be annotated? Defaults is 'no'")

parser.add_argument('-t', '--time',
                    dest="time",
                    action="store",
                    default="72",
                    help="Maximum execution time (in hours)")

parser.add_argument('-@', '--cpus',
                    dest="cpus",
                    action="store",
                    default="1",
                    help="Number of cpus per task")

parser.add_argument('-tasks', '--tasks',
                    dest="tasks",
                    action="store",
                    default="",
                    help="Number of tasks in greasy")

parser.add_argument('-n', '--analysisName',
                    dest="analysisName",
                    action="store",
                    default="",
                    help="Analysis name")

options = parser.parse_args()


# Starting time
date_str = time.strftime("%Y/%m/%d").replace("/", "")


# define variables:
refGenome = options.pathToReferenceData + "/data/genome_GRCh38.p13_GCA_000001405.28/BWA_and_PICARD/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
refGenomeFolder = options.pathToReferenceData + "/data/genome_GRCh38.p13_GCA_000001405.28/BWA_and_PICARD/"
dbSNP = options.pathToReferenceData+"/data/snpEff/db/GRCh38/dbSnp/00-All.vcf.gz"
dbNSFP = options.pathToReferenceData + "/data/snpEff/db/GRCh38/dbSNFP/dbNSFP4.3a.txt.gz"
dbCOSMIC = options.pathToReferenceData + "/data/snpEff/db/GRCh38/cosmic/sorted_subs_reheaded_joined_Cosmic_Coding_And_Noncoding.vcf.bgz"
db1kGPhg38 = options.pathToReferenceData + "/data/snpEff/db/GRCh38/1kGPhg38/ERZ822766_merged.vcf.gz"
dbgnomAD_genomes =  options.pathToReferenceData + "/data/snpEff/db/GRCh38/gnomad/release/4.0/vcf/genomes/merged_gnomad.genomes.v4.0.vcf.bgz"

projectName = options.project

# Create output directory
aName = "_"+options.analysisName if options.analysisName != "" else ""
outDir = options.workDir+"/"+date_str + "_pipeline_vcf_variantAnnotation_African_vars_genomes"+aName
os.makedirs(outDir)
os.makedirs(outDir+"/tmp")
os.makedirs(outDir+"/log")
os.makedirs(outDir+"/snpeff")


# Prepare commands
listVcfs = [fq for fq in glob.glob(options.inDir+"/*.vcf")]

if options.tasks != "":
    nTasks = options.tasks
else:
    nTasks = str(len(listVcfs)) if len(listVcfs) <= 20 else "20"


# Create bedList
comm1 = "java -Xmx8G -jar /apps/PICARD/2.24.0/bin/picard.jar BedToIntervalList -I " + options.targetBedFile+" -O "+outDir + "/tmp/targetsCoordinatesPicard.txt -SD "+refGenome+".dict"


# prepare commands
O = open(outDir+"/commands.txt", "w")
for vcfFile in listVcfs:
    
    O.write("python "+options.pathToReferencePipelines+"/pipeline_vcf_variantAnnotation_African_vars_genomes_2.py -s "+options.pathToReferencePipelines+" -p "+options.pathToReferenceData+" -pr "+options.project+" -vcf "+vcfFile+" -R "+refGenome+" -RF "+refGenomeFolder+" -dbSNPdatabase "+dbSNP+" -db1kGPhg38database "+db1kGPhg38+" -dbgnomADgenDatabase "+dbgnomAD_genomes+" -tBed "+options.targetBedFile+" -tBedIL "+outDir+"/tmp/targetsCoordinatesPicard.txt -trans "+options.transcriptsToUse+" -o "+outDir+" -@ "+options.cpus+" -intron "+options.intronMutations+" -intergenic "+options.intergenicMutations+" -downstream "+options.downstreamMutations+" -upstream "+options.upstreamMutations+"\n")
O.close()


# Prepare script
O = open(outDir+"/job_script.sh", "w")
O.write("#!/bin/bash\n")
O.write("#SBATCH --job-name=gen_variantAnnot\n")
O.write("#SBATCH -D '.'\n")
O.write("#SBATCH --output=%j.out\n")
O.write("#SBATCH --error=%j.err\n")
O.write("#SBATCH --ntasks="+nTasks+"\n")
O.write("#SBATCH --time="+options.time+":00:00\n")
O.write("#SBATCH --cpus-per-task="+options.cpus+"\n")
if int(options.time) <= 2:
    O.write("#SBATCH --qos=debug\n")
O.write("\nmodule load python/3.6.5 greasy/2.2 samtools/1.9 dotnet/2.2.402 pisces/5.3.0.0-FORK java/12.0.2 gatk/4.1.8.1 htslib/1.10.2 bcftools/1.10.2 snpeff/5.1 picard/2.24.0\n\n")
O.write(comm1+"\n\n")
O.write("greasy commands.txt\n\n")
O.write("python "+options.pathToReferencePipelines+"/pipeline_vcf_variantAnnotation_African_vars_genomes_3.py -p " +outDir+"/snpeff/"+" -pr "+options.project+" -n "+aName+"\n\n")
O.close()


# Run job
os.chdir(outDir)
subprocess.call("sbatch job_script.sh", shell=True)
