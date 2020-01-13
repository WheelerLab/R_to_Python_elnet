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
args = parser.parse_args()
chrom = str(args.chr)
chunk = str(args.chunk)

for_header = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk1.txt.gz", sep="\t", nrows=2)
header=list(for_header.columns)
#because the headers in the chroms change between header and no header
#I use condition to capture it

if header[0] == "id":
    #read without header
    all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t")
    all1.columns = header #put header

    tr_chr1 = pd.read_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr"+chrom+".txt",sep="\t")

    snps = list(all1["id"])

    tr_chr1_chunk = tr_chr1[(tr_chr1["id"].isin(snps))]

    #this is slow
    tr_chr1_chunk.to_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr"+chrom+"_chunk"+chunk+".txt", header=True, index=False, sep="\t")
    #index=False is rownames=False

else:
    #read in chrom1 and take its header
    for_header = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr1_genotype_chunk1.txt.gz", sep="\t", nrows=2)
    header=list(for_header.columns)
    #read without header
    all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t", header=None)
    all1.columns = header #put header

    tr_chr1 = pd.read_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr"+chrom+".txt",sep="\t")

    snps = list(all1["id"])

    tr_chr1_chunk = tr_chr1[(tr_chr1["id"].isin(snps))]

    #this is slow
    tr_chr1_chunk.to_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr"+chrom+"_chunk"+chunk+".txt", header=True, index=False, sep="\t")
    #index=False is rownames=False


"""
#read without header
all1 = pd.read_csv("/home/pokoro/data/mesa_models/all/ALL_chr"+chrom+"_genotype_chunk"+chunk+".txt.gz", sep="\t", header=None)
all1.columns = header #put header

tr_chr1 = pd.read_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_dosage_chr"+chrom+".txt",sep="\t")

snps = list(all1["id"])

tr_chr1_chunk = tr_chr1[(tr_chr1["id"].isin(snps))]

#this is slow
tr_chr1_chunk.to_csv("/home/pokoro/data/mesa_models/mesa_pheno/thrombotic/cau_imputation_chr"+chrom+"_chunk"+chunk+".txt", header=True, index=False, sep="\t")
#index=False is rownames=False
"""

#This is very fast
#tr_chr1_chunk.to_hdf("tr_chr_chunk3.txt", header=True, index=False, key="stage")

#However, it writes the object as a binary file, and can only be read by
#pd.read_hdf()
#this makes it hard to deal the df for model building.
#I will figure it out later

#read_tr2 = pd.read_hdf("tr_chr_chunk3.txt")

#That will require me to make the few changes below in get_maf_filtered_genotype

"""
def get_maf_filtered_genotype(genotype_file_name,  maf):
	gt_df = pd.read_hdf(genotype_file_name)
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

def get_maf_filtered_genotype_csv(genotype_file_name,  maf):
	gt_df = pd.read_csv(genotype_file_name, 'r', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

test_gt_df1.drop("id", inplace=True) #Just to drop the id column
test_gt_df1 = test_gt_df1.rename(columns=test_gt_df1.iloc[0]) #to remove the row
#that was used as the column
"""
