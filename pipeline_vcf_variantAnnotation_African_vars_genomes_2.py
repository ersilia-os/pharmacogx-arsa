# Libraries:
import subprocess 
import argparse
import time
import sys
import os


# Inputs:
parser  = argparse.ArgumentParser(prog='pipeline_vcf_variantAnnotation_African_vars_genomes_2', description='''Pipeline to run variant annotation''')


parser.add_argument('-s', '--pathToReferencePipelines',
					dest = "pathToReferencePipelines",
					action = "store",
					default = "/slgpfs/projects/cli79/reference_pipelines",
					help = "Path to folder reference_pipelines")

parser.add_argument('-p', '--pathToReferenceData',
					dest = "pathToReferenceData",
					action = "store",
					default = "/slgpfs/projects/cli79/reference_data",
					help = "Path to folder reference_data")

parser.add_argument('-vcf', '--vcfFile',
					dest = "vcfFile",
					action = "store",
					help = "Path to vcf file") 

parser.add_argument('-v', '--genomeVersion', 
					dest = "genomeVersion",
					action = "store",
					default = "hg38",
					choices=['hg19', 'hg38'],
					help = "Genome version")

parser.add_argument('-R', '--refGenome',
					dest = "refGenome",
					action = "store",
					help = "Path to reference genome")

parser.add_argument('-RF', '--refGenomeFolder',
					dest = "refGenomeFolder",
					action = "store",
					help = "Path to reference genome folder")

parser.add_argument('-tBed', '--targetBedFile',
					dest = "targetBedFile",
					required=True,
					action = "store",
					help = "Target bed file with the regions of interest for variant calling")

parser.add_argument('-tBedIL', '--targetBedFileIntervalList',
					dest = "targetBedFileIntervalList",
					required=True,
					action = "store",
					help = "Target bed interval list file with the regions of interest for variant calling")

parser.add_argument('-trans', '--transcriptsToUse',
					dest = "transcriptsToUse",
					default = "canonical",
					action = "store",
					help = "Txt file with transcripts to use for annotations")

parser.add_argument('-o', '--outDir',
					dest = "outDir",
					action = "store",
					help = "Working directory")

parser.add_argument('-@', '--cpus',
					dest = "cpus",
					action = "store",
					default="1",
					help = "Number of cpus")

parser.add_argument('-pr', '--project',
					dest = "project",
					action = "store",
					choices= ['1kGPhg38', 'PharmGKB', 'Fedorova'],
					help = "Project name")

parser.add_argument('-intron', '--intronMutations',
					dest = "intronMutations",
					default = "yes",
					action = "store",
					choices=['yes', 'no'],
					help = "Should intron mutations be annotated?")

parser.add_argument('-intergenic', '--intergenicMutations',
					dest = "intergenicMutations",
					default = "no",
					action = "store",
					choices=['yes', 'no'],
					help = "Should intergenic mutations be annotated? Defaults is 'no'")

parser.add_argument('-downstream', '--downstreamMutations',
					dest = "downstreamMutations",
					default = "no",
					action = "store",
					choices=['yes', 'no'],
					help = "Should downstream mutations be annotated? Defaults is 'no'")

parser.add_argument('-upstream', '--upstreamMutations',
					dest = "upstreamMutations",
					default = "no",
					action = "store",
					choices=['yes', 'no'],
					help = "Should upstream mutations be annotated? Defaults is 'no'")

parser.add_argument('-dbSNPdatabase', '--dbSNP',
 					dest = "dbSNP",
 					action = "store",
 					help = "Path to dbSNP database")

parser.add_argument('-dbgnomADgenDatabase', '--dbgnomAD_genomes',
 					dest = "dbgnomAD_genomes",
 					action = "store",
 					help = "Path to gnomAD genomes database")

parser.add_argument('-db1kGPhg38database', '--db1kGPhg38',
					dest = "db1kGPhg38",
					action = "store",
					help = "Path to db1kGPhg38 database")

parser.add_argument('-snpEffdataGenome', '--snpEffdata',
					dest = "snpEffdata",
					action = "store",
					help = "Path to snpEff genome data")


options = parser.parse_args()


##
## 1) DEFINE VARIABLES AND FILES
##
pathToReferencePipelines = options.pathToReferencePipelines
pathToReferenceData = options.pathToReferenceData
vcfFile = options.vcfFile
genomeVersion = options.genomeVersion
refGenome = options.refGenome
refGenomeFolder = options.refGenomeFolder
targetBedFile = options.targetBedFile
targetBedFileIntervalList = options.targetBedFileIntervalList
transcriptsToUse = options.transcriptsToUse
outDir = options.outDir
cpus = options.cpus
intronMutations = options.intronMutations
intergenicMutations = options.intergenicMutations
downstreamMutations = options.downstreamMutations
upstreamMutations = options.upstreamMutations
dbSNP = options.dbSNP
db1kGPhg38 = options.db1kGPhg38
dbgnomAD_genomes = options.dbgnomAD_genomes
snpEffdata = options.snpEffdata
projectName = options.project

start_time = time.time()
sample = vcfFile.split("/")[-1].split(".")[0]
LOG_FILE_OUT = open(outDir+"/log/log_out_"+sample+".txt", 'w')
LOG_FILE_ERR = open(outDir+"/log/log_err_"+sample+".txt", 'w')



##
## 5) ANNOTATION (snpeff + snpsift)
### -canon => to use only canonical genes
### -onlyTr my_transcripts.txt => to annote only transcripts of interest (txt file with one transcript ID per line)
### 5.1) snpsift + snpeff

vcfAnnot =  outDir+"/snpeff/"+sample+"_ann.vcf"

fields1kGPhg38 = "AF,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF"

fields_gnomAD_genomes = "AF_afr,AF_amr,AF_asj,AF_eas,AF_fin,AF_mid,AF_nfe,AF_sas,AF_remaining,AF_non_ukb_afr,AF_non_ukb_amr,AF_non_ukb_asj,AF_non_ukb_eas,AF_non_ukb_fin,AF_non_ukb_mid,AF_non_ukb_nfe,AF_non_ukb_sas,AF_non_ukb_remaining,AF_joint_afr,AF_joint_amr,AF_joint_asj,AF_joint_eas,AF_joint_fin,AF_joint_mid,AF_joint_nfe,AF_joint_sas,AF_joint_remaining,AF_joint_ami"

dataDir = pathToReferenceData+"/data/snpEff/data/"

## SnpSift annotate
# genomes
bashArguments =  "snpsift annotate -noLog -a -tabix -id -noInfo "+dbSNP+" "+vcfFile+" | snpsift annotate -noLog -a -tabix -noId -info "+fields1kGPhg38+" "+db1kGPhg38+" /dev/stdin | snpsift annotate -noLog -a -tabix -noId -info "+fields_gnomAD_genomes+" "+dbgnomAD_genomes+" /dev/stdin | snpEff -noStats -noLog "+("" if intronMutations == "yes" else "-no-intron")+(" " if intergenicMutations == "yes" else " -no-intergenic")+(" " if downstreamMutations == "yes" else " -no-downstream")+(" " if upstreamMutations == "yes" else " -no-upstream")+(" -canon" if transcriptsToUse == "canonical" else " -onlyTr "+transcriptsToUse)+" -dataDir "+dataDir+" GRCh38.105 > "+vcfAnnot

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
	LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
	LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
	sys.exit(1)


### 5.2) snpsift extractFields
fieldsToExtract_1 = "CHROM POS REF ALT" 

fieldsToExtract_2 = "ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].EFFECT ANN[*].IMPACT ANN[*].RANK ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ALLELE ANN[*].ERRORS"

fields_1kGPhg38 = "ID EAS_AF EUR_AF AFR_AF AMR_AF SAS_AF"

fields_gnomAD_genomes = "AF_afr AF_amr AF_asj AF_eas AF_fin AF_mid AF_nfe AF_sas AF_remaining AF_joint_afr AF_joint_amr AF_joint_asj AF_joint_eas AF_joint_fin AF_joint_mid AF_joint_nfe AF_joint_sas AF_joint_remaining AF_joint_ami"

fieldsToExtract_3 = "fedorova_label adme_label abundantAFR_label maxFreqs AFR_overrepresentation AFR_overrepresentation_mod AFR_overrepresentation_label specificAFR_label"

# Extract field for genomes
tsvAnnot_tmp = vcfAnnot.replace(".vcf", "_all_annotations.tsv")
tsvAnnot = vcfAnnot.replace(".vcf", ".tsv")

bashArguments = "cat "+vcfAnnot+" | "+("" if transcriptsToUse == "canonical" else "/apps/SNPEFF/5.1/scripts/vcfEffOnePerLine.pl | ")+"snpsift extractFields -s ',' -e '.' - '"+fieldsToExtract_1+" "+fieldsToExtract_2+" "+fields_1kGPhg38+" "+fields_gnomAD_genomes+" "+fieldsToExtract_3+"' > "+tsvAnnot_tmp

e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
	LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
	LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
	sys.exit(1)



### 5.3) keep only lines with transcripts 
# For exomes
I = open(tsvAnnot_tmp, "r")
O = open(tsvAnnot, "w")
for iLine in I:
	if iLine.startswith("CHROM"): 
		list = iLine.split("\t")
		for i in range(24,29):
			list[i] = "_".join(['1kGPhg38', list[i]])
		for i in range(29,48):
			list[i] = "_".join(['gnomAD_genomes', list[i]])
		
		newLine = '\t'.join(list)
		O.write(newLine)

	else:
		iList = iLine.rstrip("\n").split("\t")
		list = iList[23].split(";")
		# Substituir la columna ID conjunt per l'ID RS
		rs = (",".join([idx for idx in list if idx.startswith('rs') or idx.startswith('"rs')])).strip('"')  
		iList[23] = rs

		newLine = '\t'.join(iList)
		
		O.write(newLine+"\n")

I.close()
O.close()
os.remove(tsvAnnot_tmp)



# 6) CLOSE LOGS
#
LOG_FILE_OUT.write("DONE. EXECUTION TIME IN MINUTES: %s\n" %(round((time.time() - start_time)/60, 4))) 
LOG_FILE_OUT.close()
LOG_FILE_ERR.close()