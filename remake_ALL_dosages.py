#Remake the ALL chunk dosages because the ALL does not have header in some of the file
import pandas as pd

#read the chrom1 and take the header, cause some chroms do not have header
#i read in only two rows to make it fast
#we will use the header to complete the other chroms
#this is to avoid key erro when i start accessing the df by col name

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("chr", action="store", help="put chromosome no")
parser.add_argument("chunk", action="store", help="put chromosome chunk no")
args = parser.parse_args()
chrom = str(args.chr)
chunk = str(args.chunk)

for_header = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t", nrows=5)
header=list(for_header.columns)
#because the headers in the chroms change between header and no header
#I use condition to capture it

if header[0] == "id":
    all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", header=0, sep="\t")
    all1.to_csv("/home/pokoro/data/mesa_models/all/chunks/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", header=True, index=False, sep="\t")

else:
    for_header = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr1_genotype_chunk1.txt.gz", sep="\t", nrows=5) #take header from chrom1    
    header=list(for_header.columns)
    #read the chrom without header
    all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t", header=None)
    all1.columns = header #put header
    all1.to_csv("/home/pokoro/data/mesa_models/all/chunks/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", header=True, index=False, sep="\t")
