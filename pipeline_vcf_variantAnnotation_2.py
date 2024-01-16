# Libraries:
import subprocess 
import argparse
import time
import sys
import os


# Inputs:
parser  = argparse.ArgumentParser(prog='pipeline_bam_targetVariantCallingRed', description='''Pipeline to run variant calling in target NGS''')


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

parser.add_argument('-dbNSFPdatabase', '--dbNSFP',
					dest = "dbNSFP",
					action = "store",
					help = "Path to dbNSFP database")

parser.add_argument('-dbCOSMICdatabase', '--dbCOSMIC',
					dest = "dbCOSMIC",
					action = "store",
					help = "Path to dbCOSMIC database")

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
dbNSFP = options.dbNSFP
dbCOSMIC = options.dbCOSMIC
db1kGPhg38 = options.db1kGPhg38
snpEffdata = options.snpEffdata
projectName = options.project

start_time = time.time()
sample = vcfFile.split("/")[-1].split(".")[0]
LOG_FILE_OUT = open(outDir+"/log/log_out_"+sample+".txt", 'w')
LOG_FILE_ERR = open(outDir+"/log/log_err_"+sample+".txt", 'w')

"""

##
## 5) ANNOTATION (snpeff + snpsift)
##
### (-no-intron) -no-intergenic -no-downstream -no-upstream
### -canon => to use only canonical genes
### -onlyTr my_transcripts.txt => to annote only transcripts of interest (txt file with one transcript ID per line)
### 5.1) snpsift + snpeff

vcfAnnot =  outDir+"/snpeff/"+sample+"_ann.vcf"

fieldsOfInterest = "gnomAD_exomes_AC,gnomAD_exomes_AN,gnomAD_exomes_AF,gnomAD_exomes_POPMAX_AC,gnomAD_exomes_POPMAX_AN,gnomAD_exomes_POPMAX_AF,gnomAD_exomes_AFR_AC,gnomAD_exomes_AFR_AN,gnomAD_exomes_AFR_AF,gnomAD_exomes_NFE_AC,gnomAD_exomes_NFE_AN,gnomAD_exomes_NFE_AF,gnomAD_exomes_AMR_AC,gnomAD_exomes_AMR_AN,gnomAD_exomes_AMR_AF,gnomAD_exomes_ASJ_AC,gnomAD_exomes_ASJ_AN,gnomAD_exomes_ASJ_AF,gnomAD_exomes_EAS_AC,gnomAD_exomes_EAS_AN,gnomAD_exomes_EAS_AF,gnomAD_exomes_FIN_AC,gnomAD_exomes_FIN_AN,gnomAD_exomes_FIN_AF,gnomAD_exomes_SAS_AC,gnomAD_exomes_SAS_AN,gnomAD_exomes_SAS_AF,gnomAD_genomes_AC,gnomAD_genomes_AN,gnomAD_genomes_AF,gnomAD_genomes_POPMAX_AC,gnomAD_genomes_POPMAX_AN,gnomAD_genomes_POPMAX_AF,gnomAD_genomes_AFR_AC,gnomAD_genomes_AFR_AN,gnomAD_genomes_AFR_AF,gnomAD_genomes_NFE_AN,gnomAD_genomes_NFE_AC,gnomAD_genomes_NFE_AF,gnomAD_genomes_AMI_AC,gnomAD_genomes_AMI_AN,gnomAD_genomes_AMI_AF,gnomAD_genomes_AMR_AC,gnomAD_genomes_AMR_AN,gnomAD_genomes_AMR_AF,gnomAD_genomes_ASJ_AC,gnomAD_genomes_ASJ_AN,gnomAD_genomes_ASJ_AF,gnomAD_genomes_EAS_AC,gnomAD_genomes_EAS_AN,gnomAD_genomes_EAS_AF,gnomAD_genomes_FIN_AN,gnomAD_genomes_FIN_AF,gnomAD_genomes_MID_AC,gnomAD_genomes_MID_AN,gnomAD_genomes_MID_AF,gnomAD_genomes_SAS_AC,gnomAD_genomes_SAS_AN,gnomAD_genomes_SAS_AF,1000Gp3_AC,1000Gp3_AF,1000Gp3_AFR_AC,1000Gp3_AFR_AF,1000Gp3_AMR_AC,1000Gp3_AMR_AF,1000Gp3_EAS_AC,1000Gp3_EAS_AF,1000Gp3_EUR_AC,1000Gp3_EUR_AF,1000Gp3_SAS_AC,1000Gp3_SAS_AF,ESP6500_AA_AC,ESP6500_AA_AF,ESP6500_EA_AC,ESP6500_EA_AF,ExAC_AC,ExAC_AF,ExAC_Adj_AC,ExAC_Adj_AF,ExAC_AFR_AC,ExAC_AFR_AF,ExAC_AMR_AC,ExAC_AMR_AF,ExAC_EAS_AC,ExAC_EAS_AF,ExAC_FIN_AC,ExAC_FIN_AF,ExAC_NFE_AC,ExAC_NFE_AF,ExAC_SAS_AC,ExAC_SAS_AF,CADD_phred,Polyphen2_HDIV_score,Polyphen2_HDIV_pred,SIFT_score,SIFT_pred,MutationAssessor_score,MutationAssessor_pred,MutationTaster_score,MutationTaster_pred,PROVEAN_pred,VEST4_score,clinvar_id,clinvar_clnsig,clinvar_trait,clinvar_review,clinvar_MedGen_id,clinvar_OMIM_id,clinvar_Orphanet_id,Interpro_domain"

fields1kGPhg38 = "AF,AC,NS,AN,EAS_AF,EUR_AF,AFR_AF,AMR_AF,SAS_AF"

dataDir = pathToReferenceData+"/data/snpEff/data/"

## SnpSift annotate
bashArguments =  "snpsift annotate -noLog -a -tabix -id -noInfo "+dbSNP+" "+vcfFile+" | snpsift annotate -noLog -a -tabix -noId -info "+fields1kGPhg38+" "+db1kGPhg38+" /dev/stdin | snpsift dbnsfp -f '"+fieldsOfInterest+"' -noLog -collapse -a -m -db "+dbNSFP+" /dev/stdin | snpEff -noStats -noLog "+("" if intronMutations == "yes" else "-no-intron")+(" " if intergenicMutations == "yes" else " -no-intergenic")+(" " if downstreamMutations == "yes" else " -no-downstream")+(" " if upstreamMutations == "yes" else " -no-upstream")+(" -canon" if transcriptsToUse == "canonical" else " -onlyTr "+transcriptsToUse)+" -dataDir "+dataDir+" GRCh38.105 > "+vcfAnnot


e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
	LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
	LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
	sys.exit(1)


"""
#vcfAnnot =  outDir+"/snpeff/"+sample+"_ann.vcf"

vcfAnnot = '/slgpfs/scratch/cli79/cli79334/projects/other/1kGPhg38/results/20231129_210444_pipeline_bam_targetVariantCallingRed_coding_1kGPhg38_test2/snpeff/'+sample+'_ann.vcf'
#vcfAnnot = '/slgpfs/scratch/cli79/cli79334/projects/other/AFR_ARSA_FedorovaL/02_snpEff_annotation/results/20231127_131226_pipeline_bam_targetVariantCallingRed_Fedorova_test3/snpeff/ARSAtableAFR_hg38_sorted_ann.vcf'

### 5.2) snpsift extractFields
tsvAnnot_tmp = vcfAnnot.replace(".vcf", "_all_annotations.tsv")
#tsvAnnot_tmp = outDir+"/snpeff/"+vcfAnnot.split('/')[-1].replace(".vcf", "_all_annotations.tsv")
tsvAnnot = vcfAnnot.replace(".vcf", ".tsv")
#tsvAnnot = outDir+"/snpeff/"+vcfAnnot.split('/')[-1].replace(".vcf", ".tsv")
fieldsToExtract_1 = "CHROM POS REF ALT" 
fieldsToExtract_2 = "ANN[*].GENE ANN[*].GENEID ANN[*].FEATURE ANN[*].FEATUREID ANN[*].BIOTYPE ANN[*].EFFECT ANN[*].IMPACT ANN[*].RANK ANN[*].HGVS_C ANN[*].HGVS_P ANN[*].CDNA_POS ANN[*].CDNA_LEN ANN[*].CDS_POS ANN[*].CDS_LEN ANN[*].AA_POS ANN[*].AA_LEN ANN[*].DISTANCE ANN[*].ALLELE ANN[*].ERRORS"


# For 1kGPhg38
if projectName == '1kGPhg38':
	fieldsToExtract_3 = "ID AF AC NS AN EAS_AF EUR_AF AFR_AF AMR_AF SAS_AF dbNSFP_gnomAD_exomes_AC dbNSFP_gnomAD_exomes_AN dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_exomes_POPMAX_AC dbNSFP_gnomAD_exomes_POPMAX_AN dbNSFP_gnomAD_exomes_POPMAX_AF dbNSFP_gnomAD_exomes_AFR_AC dbNSFP_gnomAD_exomes_AFR_AN dbNSFP_gnomAD_exomes_AFR_AF dbNSFP_gnomAD_exomes_NFE_AC dbNSFP_gnomAD_exomes_NFE_AN dbNSFP_gnomAD_exomes_NFE_AF dbNSFP_gnomAD_exomes_AMR_AC dbNSFP_gnomAD_exomes_AMR_AN dbNSFP_gnomAD_exomes_AMR_AF dbNSFP_gnomAD_exomes_ASJ_AC dbNSFP_gnomAD_exomes_ASJ_AN dbNSFP_gnomAD_exomes_ASJ_AF dbNSFP_gnomAD_exomes_EAS_AC dbNSFP_gnomAD_exomes_EAS_AN dbNSFP_gnomAD_exomes_EAS_AF dbNSFP_gnomAD_exomes_FIN_AC dbNSFP_gnomAD_exomes_FIN_AN dbNSFP_gnomAD_exomes_FIN_AF dbNSFP_gnomAD_exomes_SAS_AC dbNSFP_gnomAD_exomes_SAS_AN dbNSFP_gnomAD_exomes_SAS_AF dbNSFP_gnomAD_genomes_AC dbNSFP_gnomAD_genomes_AN dbNSFP_gnomAD_genomes_AF dbNSFP_gnomAD_genomes_POPMAX_AC dbNSFP_gnomAD_genomes_POPMAX_AN dbNSFP_gnomAD_genomes_POPMAX_AF dbNSFP_gnomAD_genomes_AFR_AC dbNSFP_gnomAD_genomes_AFR_AN dbNSFP_gnomAD_genomes_AFR_AF dbNSFP_gnomAD_genomes_NFE_AN dbNSFP_gnomAD_genomes_NFE_AC dbNSFP_gnomAD_genomes_NFE_AF dbNSFP_gnomAD_genomes_AMI_AC dbNSFP_gnomAD_genomes_AMI_AN dbNSFP_gnomAD_genomes_AMI_AF dbNSFP_gnomAD_genomes_AMR_AC dbNSFP_gnomAD_genomes_AMR_AN dbNSFP_gnomAD_genomes_AMR_AF dbNSFP_gnomAD_genomes_ASJ_AC dbNSFP_gnomAD_genomes_ASJ_AN dbNSFP_gnomAD_genomes_ASJ_AF dbNSFP_gnomAD_genomes_EAS_AC dbNSFP_gnomAD_genomes_EAS_AN dbNSFP_gnomAD_genomes_EAS_AF dbNSFP_gnomAD_genomes_FIN_AN dbNSFP_gnomAD_genomes_FIN_AF dbNSFP_gnomAD_genomes_MID_AC dbNSFP_gnomAD_genomes_MID_AN dbNSFP_gnomAD_genomes_MID_AF dbNSFP_gnomAD_genomes_SAS_AC dbNSFP_gnomAD_genomes_SAS_AN dbNSFP_gnomAD_genomes_SAS_AF dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_1000Gp3_AFR_AC dbNSFP_1000Gp3_AFR_AF dbNSFP_1000Gp3_AMR_AC dbNSFP_1000Gp3_AMR_AF dbNSFP_1000Gp3_EAS_AC dbNSFP_1000Gp3_EAS_AF dbNSFP_1000Gp3_EUR_AC dbNSFP_1000Gp3_EUR_AF dbNSFP_1000Gp3_SAS_AC dbNSFP_1000Gp3_SAS_AF dbNSFP_ESP6500_AA_AC dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AC dbNSFP_ESP6500_EA_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF dbNSFP_ExAC_Adj_AC dbNSFP_ExAC_Adj_AF dbNSFP_ExAC_AFR_AC dbNSFP_ExAC_AFR_AF dbNSFP_ExAC_AMR_AC dbNSFP_ExAC_AMR_AF dbNSFP_ExAC_EAS_AC dbNSFP_ExAC_EAS_AF dbNSFP_ExAC_FIN_AC dbNSFP_ExAC_FIN_AF dbNSFP_ExAC_NFE_AC dbNSFP_ExAC_NFE_AF dbNSFP_ExAC_SAS_AC dbNSFP_ExAC_SAS_AF dbNSFP_CADD_phred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_PROVEAN_pred dbNSFP_VEST4_score dbNSFP_clinvar_id dbNSFP_clinvar_clnsig dbNSFP_clinvar_trait dbNSFP_clinvar_review dbNSFP_clinvar_MedGen_id dbNSFP_clinvar_OMIM_id dbNSFP_clinvar_Orphanet_id dbNSFP_Interpro_domain"


# For PharmGkB
elif projectName == 'PharmGKB': 
	fieldsToExtract_3 = "ID AF AC NS AN EAS_AF EUR_AF AFR_AF AMR_AF SAS_AF PGKBVID RSID GENE PGKBGID HGVS dbNSFP_gnomAD_exomes_AC dbNSFP_gnomAD_exomes_AN dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_exomes_POPMAX_AC dbNSFP_gnomAD_exomes_POPMAX_AN dbNSFP_gnomAD_exomes_POPMAX_AF dbNSFP_gnomAD_exomes_AFR_AC dbNSFP_gnomAD_exomes_AFR_AN dbNSFP_gnomAD_exomes_AFR_AF dbNSFP_gnomAD_exomes_NFE_AC dbNSFP_gnomAD_exomes_NFE_AN dbNSFP_gnomAD_exomes_NFE_AF dbNSFP_gnomAD_exomes_AMR_AC dbNSFP_gnomAD_exomes_AMR_AN dbNSFP_gnomAD_exomes_AMR_AF dbNSFP_gnomAD_exomes_ASJ_AC dbNSFP_gnomAD_exomes_ASJ_AN dbNSFP_gnomAD_exomes_ASJ_AF dbNSFP_gnomAD_exomes_EAS_AC dbNSFP_gnomAD_exomes_EAS_AN dbNSFP_gnomAD_exomes_EAS_AF dbNSFP_gnomAD_exomes_FIN_AC dbNSFP_gnomAD_exomes_FIN_AN dbNSFP_gnomAD_exomes_FIN_AF dbNSFP_gnomAD_exomes_SAS_AC dbNSFP_gnomAD_exomes_SAS_AN dbNSFP_gnomAD_exomes_SAS_AF dbNSFP_gnomAD_genomes_AC dbNSFP_gnomAD_genomes_AN dbNSFP_gnomAD_genomes_AF dbNSFP_gnomAD_genomes_POPMAX_AC dbNSFP_gnomAD_genomes_POPMAX_AN dbNSFP_gnomAD_genomes_POPMAX_AF dbNSFP_gnomAD_genomes_AFR_AC dbNSFP_gnomAD_genomes_AFR_AN dbNSFP_gnomAD_genomes_AFR_AF dbNSFP_gnomAD_genomes_NFE_AN dbNSFP_gnomAD_genomes_NFE_AC dbNSFP_gnomAD_genomes_NFE_AF dbNSFP_gnomAD_genomes_AMI_AC dbNSFP_gnomAD_genomes_AMI_AN dbNSFP_gnomAD_genomes_AMI_AF dbNSFP_gnomAD_genomes_AMR_AC dbNSFP_gnomAD_genomes_AMR_AN dbNSFP_gnomAD_genomes_AMR_AF dbNSFP_gnomAD_genomes_ASJ_AC dbNSFP_gnomAD_genomes_ASJ_AN dbNSFP_gnomAD_genomes_ASJ_AF dbNSFP_gnomAD_genomes_EAS_AC dbNSFP_gnomAD_genomes_EAS_AN dbNSFP_gnomAD_genomes_EAS_AF dbNSFP_gnomAD_genomes_FIN_AN dbNSFP_gnomAD_genomes_FIN_AF dbNSFP_gnomAD_genomes_MID_AC dbNSFP_gnomAD_genomes_MID_AN dbNSFP_gnomAD_genomes_MID_AF dbNSFP_gnomAD_genomes_SAS_AC dbNSFP_gnomAD_genomes_SAS_AN dbNSFP_gnomAD_genomes_SAS_AF dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_1000Gp3_AFR_AC dbNSFP_1000Gp3_AFR_AF dbNSFP_1000Gp3_AMR_AC dbNSFP_1000Gp3_AMR_AF dbNSFP_1000Gp3_EAS_AC dbNSFP_1000Gp3_EAS_AF dbNSFP_1000Gp3_EUR_AC dbNSFP_1000Gp3_EUR_AF dbNSFP_1000Gp3_SAS_AC dbNSFP_1000Gp3_SAS_AF dbNSFP_ESP6500_AA_AC dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AC dbNSFP_ESP6500_EA_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF dbNSFP_ExAC_Adj_AC dbNSFP_ExAC_Adj_AF dbNSFP_ExAC_AFR_AC dbNSFP_ExAC_AFR_AF dbNSFP_ExAC_AMR_AC dbNSFP_ExAC_AMR_AF dbNSFP_ExAC_EAS_AC dbNSFP_ExAC_EAS_AF dbNSFP_ExAC_FIN_AC dbNSFP_ExAC_FIN_AF dbNSFP_ExAC_NFE_AC dbNSFP_ExAC_NFE_AF dbNSFP_ExAC_SAS_AC dbNSFP_ExAC_SAS_AF dbNSFP_CADD_phred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_PROVEAN_pred dbNSFP_VEST4_score dbNSFP_clinvar_id dbNSFP_clinvar_clnsig dbNSFP_clinvar_trait dbNSFP_clinvar_review dbNSFP_clinvar_MedGen_id dbNSFP_clinvar_OMIM_id dbNSFP_clinvar_Orphanet_id dbNSFP_Interpro_domain"

elif projectName == 'Fedorova': 
	fieldsToExtract_3 = "ID AF AC NS AN EAS_AF EUR_AF AFR_AF AMR_AF SAS_AF Simons_AFR Simons_OTHER1 Simons_OTHER2 Simons_OTHER3 1000Gp3_AFR 1000Gp3_EUR 1000Gp3_EAS 1000Gp3_AMR 1000Gp3_SAS dbNSFP_gnomAD_exomes_AC dbNSFP_gnomAD_exomes_AN dbNSFP_gnomAD_exomes_AF dbNSFP_gnomAD_exomes_POPMAX_AC dbNSFP_gnomAD_exomes_POPMAX_AN dbNSFP_gnomAD_exomes_POPMAX_AF dbNSFP_gnomAD_exomes_AFR_AC dbNSFP_gnomAD_exomes_AFR_AN dbNSFP_gnomAD_exomes_AFR_AF dbNSFP_gnomAD_exomes_NFE_AC dbNSFP_gnomAD_exomes_NFE_AN dbNSFP_gnomAD_exomes_NFE_AF dbNSFP_gnomAD_exomes_AMR_AC dbNSFP_gnomAD_exomes_AMR_AN dbNSFP_gnomAD_exomes_AMR_AF dbNSFP_gnomAD_exomes_ASJ_AC dbNSFP_gnomAD_exomes_ASJ_AN dbNSFP_gnomAD_exomes_ASJ_AF dbNSFP_gnomAD_exomes_EAS_AC dbNSFP_gnomAD_exomes_EAS_AN dbNSFP_gnomAD_exomes_EAS_AF dbNSFP_gnomAD_exomes_FIN_AC dbNSFP_gnomAD_exomes_FIN_AN dbNSFP_gnomAD_exomes_FIN_AF dbNSFP_gnomAD_exomes_SAS_AC dbNSFP_gnomAD_exomes_SAS_AN dbNSFP_gnomAD_exomes_SAS_AF dbNSFP_gnomAD_genomes_AC dbNSFP_gnomAD_genomes_AN dbNSFP_gnomAD_genomes_AF dbNSFP_gnomAD_genomes_POPMAX_AC dbNSFP_gnomAD_genomes_POPMAX_AN dbNSFP_gnomAD_genomes_POPMAX_AF dbNSFP_gnomAD_genomes_AFR_AC dbNSFP_gnomAD_genomes_AFR_AN dbNSFP_gnomAD_genomes_AFR_AF dbNSFP_gnomAD_genomes_NFE_AN dbNSFP_gnomAD_genomes_NFE_AC dbNSFP_gnomAD_genomes_NFE_AF dbNSFP_gnomAD_genomes_AMI_AC dbNSFP_gnomAD_genomes_AMI_AN dbNSFP_gnomAD_genomes_AMI_AF dbNSFP_gnomAD_genomes_AMR_AC dbNSFP_gnomAD_genomes_AMR_AN dbNSFP_gnomAD_genomes_AMR_AF dbNSFP_gnomAD_genomes_ASJ_AC dbNSFP_gnomAD_genomes_ASJ_AN dbNSFP_gnomAD_genomes_ASJ_AF dbNSFP_gnomAD_genomes_EAS_AC dbNSFP_gnomAD_genomes_EAS_AN dbNSFP_gnomAD_genomes_EAS_AF dbNSFP_gnomAD_genomes_FIN_AN dbNSFP_gnomAD_genomes_FIN_AF dbNSFP_gnomAD_genomes_MID_AC dbNSFP_gnomAD_genomes_MID_AN dbNSFP_gnomAD_genomes_MID_AF dbNSFP_gnomAD_genomes_SAS_AC dbNSFP_gnomAD_genomes_SAS_AN dbNSFP_gnomAD_genomes_SAS_AF dbNSFP_1000Gp3_AC dbNSFP_1000Gp3_AF dbNSFP_1000Gp3_AFR_AC dbNSFP_1000Gp3_AFR_AF dbNSFP_1000Gp3_AMR_AC dbNSFP_1000Gp3_AMR_AF dbNSFP_1000Gp3_EAS_AC dbNSFP_1000Gp3_EAS_AF dbNSFP_1000Gp3_EUR_AC dbNSFP_1000Gp3_EUR_AF dbNSFP_1000Gp3_SAS_AC dbNSFP_1000Gp3_SAS_AF dbNSFP_ESP6500_AA_AC dbNSFP_ESP6500_AA_AF dbNSFP_ESP6500_EA_AC dbNSFP_ESP6500_EA_AF dbNSFP_ExAC_AC dbNSFP_ExAC_AF dbNSFP_ExAC_Adj_AC dbNSFP_ExAC_Adj_AF dbNSFP_ExAC_AFR_AC dbNSFP_ExAC_AFR_AF dbNSFP_ExAC_AMR_AC dbNSFP_ExAC_AMR_AF dbNSFP_ExAC_EAS_AC dbNSFP_ExAC_EAS_AF dbNSFP_ExAC_FIN_AC dbNSFP_ExAC_FIN_AF dbNSFP_ExAC_NFE_AC dbNSFP_ExAC_NFE_AF dbNSFP_ExAC_SAS_AC dbNSFP_ExAC_SAS_AF dbNSFP_CADD_phred dbNSFP_Polyphen2_HDIV_score dbNSFP_Polyphen2_HDIV_pred dbNSFP_SIFT_score dbNSFP_SIFT_pred dbNSFP_MutationAssessor_score dbNSFP_MutationAssessor_pred dbNSFP_MutationTaster_score dbNSFP_MutationTaster_pred dbNSFP_PROVEAN_pred dbNSFP_VEST4_score dbNSFP_clinvar_id dbNSFP_clinvar_clnsig dbNSFP_clinvar_trait dbNSFP_clinvar_review dbNSFP_clinvar_MedGen_id dbNSFP_clinvar_OMIM_id dbNSFP_clinvar_Orphanet_id dbNSFP_Interpro_domain"

else:
	print('Project not available.')


# Extract fields
bashArguments = "cat "+vcfAnnot+" | "+("" if transcriptsToUse == "canonical" else "/apps/SNPEFF/5.1/scripts/vcfEffOnePerLine.pl | ")+"snpsift extractFields -s ',' -e '.' - '"+fieldsToExtract_1+" "+fieldsToExtract_2+" "+fieldsToExtract_3+"' > "+tsvAnnot_tmp
e = subprocess.call(bashArguments, shell=True, stdout=LOG_FILE_OUT, stderr=LOG_FILE_ERR)
if e != 0:
	LOG_FILE_ERR.write("ERROR IN COMMAND: "+bashArguments+"\n")
	LOG_FILE_OUT.write("ERROR IN COMMAND: "+bashArguments+"\n")
	sys.exit(1)


### 5.3) keep only lines with transcripts 
# For 1kGPhg38
if projectName == '1kGPhg38':
	I = open(tsvAnnot_tmp, "r")
	O = open(tsvAnnot, "w")
	for iLine in I:
		if iLine.startswith("CHROM"): 
			list = iLine.split("\t")
			newLine = '\t'.join(list)
			O.write(newLine)

		else:
			iList = iLine.rstrip("\n").split("\t")
			list = iList[23].split(";")
			# Substituir la columna ID conjunt per l'ID RS
			rs = (",".join([idx for idx in list if idx.startswith('rs') or idx.startswith('"rs')])).strip('"')  
			iList[23] = rs

			# Keep only the most damaging score in case more than one is present (Not necessary for CADD scores)
			for i in range(125,133,2):
				scores = iList[i].split(",")
				if i != 127:
					# Polyphen2_HDIV, MutAssessor, MutTaster
					iList[i] = max(scores)
				else:
					# SIFT
					iList[i] = min(scores)
				
			# categories...
			iList[126] = "D" if "D" in iList[126] else "P" if "P" in iList[126] else "B" if "B" in iList[126] else "." # Polyphen2_HDIV
			iList[128] = "D" if "D" in iList[128] else "T" if "T" in iList[128] else "." # SIFT
			iList[130] = "H" if "H" in iList[130] else "M" if "M" in iList[130] else "L" if "L" in iList[130] else "N" if "N" in iList[130] else "." # MutAssessor
			iList[132] = "A" if "A" in iList[132] else "D" if "D" in iList[132] else "N" if "N" in iList[132] else "P" if "P" in iList[132] else "." # MutTaster
			iList[133] = "D" if "D" in iList[133] else "N" if "N" in iList[133] else "." # Provean

			newLine = '\t'.join(iList)
			
			O.write(newLine+"\n")

	I.close()
	O.close()
	#os.remove(tsvAnnot_tmp)

elif projectName == 'PharmGKB':

	#PharmGKB
	I = open(tsvAnnot_tmp, "r")
	O = open(tsvAnnot, "w")
	for iLine in I:
		if iLine.startswith("CHROM"): 
			list = iLine.split("\t")
			newLine = '\t'.join(list)
			O.write(newLine)

		else:
			# Separar RS i COsv
			iList = iLine.rstrip("\n").split("\t")
			list = iList[23].split(";")
			# Substituir la columna ID conjunt per l'ID RS
			rs = (",".join([idx for idx in list if idx.startswith('rs') or idx.startswith('"rs')])).strip('"')  
			iList[23] = rs

			# Mantenir nomes els scores mes damaging (CADD no cal tocar-lo)
			for i in range(130,138,2):
				scores = iList[i].split(",")
				if i != 132:
					# Polyphen2_HDIV, MutAssessor, MutTaster
					iList[i] = max(scores)
				else:
					# SIFT
					iList[i] = min(scores)
					#iList[i] = min(scores)
				
			# categories...
			iList[131] = "D" if "D" in iList[131] else "P" if "P" in iList[131] else "B" if "B" in iList[131] else "." # Polyphen2_HDIV_pred
			iList[133] = "D" if "D" in iList[133] else "T" if "T" in iList[133] else "." # SIFT
			iList[135] = "H" if "H" in iList[135] else "M" if "M" in iList[135] else "L" if "L" in iList[135] else "N" if "N" in iList[135] else "." # MutAssessor
			iList[137] = "A" if "A" in iList[137] else "D" if "D" in iList[137] else "N" if "N" in iList[137] else "P" if "P" in iList[137] else "." # MutTaster
			iList[138] = "D" if "D" in iList[138] else "N" if "N" in iList[138] else "." # Provean

			newLine = '\t'.join(iList)
			
			O.write(newLine+"\n")

	I.close()
	O.close()
	os.remove(tsvAnnot_tmp)

else:
	# Fedorova
	I = open(tsvAnnot_tmp, "r")
	O = open(tsvAnnot, "w")
	for iLine in I:
		if iLine.startswith("CHROM"): 
			list = iLine.split("\t")
			newLine = '\t'.join(list)
			O.write(newLine)

		else:
			# Separar RS i COsv
			iList = iLine.rstrip("\n").split("\t")
			list = iList[23].split(";")
			# Substituir la columna ID conjunt per l'ID RS
			rs = (",".join([idx for idx in list if idx.startswith('rs') or idx.startswith('"rs')])).strip('"')  
			iList[23] = rs

			# Mantenir nomes els scores mes damaging (CADD no cal tocar-lo)
			for i in range(134,142,2):
				scores = iList[i].split(",")
				if i != 136:
					# Polyphen2_HDIV, MutAssessor, MutTaster
					iList[i] = max(scores)
				else:
					# SIFT
					iList[i] = min(scores)

				
			# categories...
			iList[135] = "D" if "D" in iList[135] else "P" if "P" in iList[135] else "B" if "B" in iList[135] else "." # Polyphen2_HDIV
			iList[137] = "D" if "D" in iList[137] else "T" if "T" in iList[137] else "." # SIFT
			iList[139] = "H" if "H" in iList[139] else "M" if "M" in iList[139] else "L" if "L" in iList[139] else "N" if "N" in iList[139] else "." # MutAssessor
			iList[141] = "A" if "A" in iList[141] else "D" if "D" in iList[141] else "N" if "N" in iList[141] else "P" if "P" in iList[141] else "." # MutTaster
			iList[142] = "D" if "D" in iList[142] else "N" if "N" in iList[142] else "." # Provean

			newLine = '\t'.join(iList)
			
			O.write(newLine+"\n")

	I.close()
	O.close()
	#os.remove(tsvAnnot_tmp)



# 6) CLOSE LOGS
#
LOG_FILE_OUT.write("DONE. EXECUTION TIME IN MINUTES: %s\n" %(round((time.time() - start_time)/60, 4))) 
LOG_FILE_OUT.close()
LOG_FILE_ERR.close()