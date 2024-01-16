import glob
import argparse

parser  = argparse.ArgumentParser(prog='pipeline_bam_targetVariantCallingRed', description='''Pipeline to run variant calling in target NGS (script 3-merge)''')

parser.add_argument('-p', '--path',
                    dest = "path",
                    action = "store",
                    help = "Path to folder containing snpeff tsv files")

parser.add_argument('-n', '--name',
                    dest = "name",
                    action = "store",
                    help = "Run name") 

parser.add_argument('-pr', '--project',
                    dest = "project",
                    action = "store",
                    help = "Project name")

options = parser.parse_args()

projectName = options.project

O = open(options.path.replace("snpeff/", "mutations"+options.name+".tsv"), "w")
headLine = "no"

for tsv in glob.glob(options.path+"*.tsv"):

    if projectName == '1kGPhg38':
        sampleName = tsv.split("/")[-1].split("_")[3].replace(".tsv", "")
    else:
        sampleName = tsv.split("/")[-1].replace(".tsv", "") 

    TSV = open(tsv, "r")
    for tLine in TSV:
        if tLine.startswith("CHROM"):
            if headLine == "no":
                O.write("SAMPLE\t%s" %tLine)
                headLine = "yes"
        else:
            O.write("%s\t%s" %(sampleName, tLine))
    TSV.close()

O.close()