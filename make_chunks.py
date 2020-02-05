#break the cau dosages into five chunks like the ALL dosages
#this is to speed up imputation

import pandas as pd

#read the chrom1 and take the header, cause some chroms do not have header
#i read in only two rows to make it fast
#we will use the header to complete the other chroms
#this is to avoid key erro when i start accessing the df by col name

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("chr", action="store", help="put chromosome no")
parser.add_argument("chunk", action="store", help="put chromosome chunk no")
parser.add_argument("pops", action="store", help="put populations name")
args = parser.parse_args()
chrom = str(args.chr)
chunk = str(args.chunk)
pop = str(args.pops)

#read in the all chunks and use it to break the new dosages
all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/chunks/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t")

tr_chr1 = pd.read_csv("/home/pokoro/data/lauren_mesa/ml_dosages/"+pop+"/chr"+chrom+".txt",sep="\t")

snps = list(all1["id"])

tr_chr1_chunk = tr_chr1[(tr_chr1["id"].isin(snps))]

#this is slow
tr_chr1_chunk.to_csv("/home/pokoro/data/lauren_mesa/ml_dosages/chunks/"+pop+"/chr"+chrom+"_chunk"+chunk+".txt.gz", header=True, index=False, sep="\t")
#index=False is rownames=False
