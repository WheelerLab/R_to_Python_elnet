import numpy as np
#from sklearn.svm import SVR
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean
import math
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame
#import pickle
#from sklearn.ensemble import RandomForestRegressor
#from sklearn.neighbors import KNeighborsRegressor
#import time
from scipy import stats
import argparse
from sklearn.linear_model import ElasticNet
#from sklearn.model_selection import GridSearchCV
#from sklearn.metrics import mean_squared_error
#from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score
from sklearn.metrics import make_scorer#use to convert metrics to scoring callables
#from hyperopt import fmin, tpe, hp, Trials, STATUS_OK
#import hyperopt

r2 = make_scorer(r2_score, greater_is_better=True)

parser = argparse.ArgumentParser()
parser.add_argument("chr", action="store", help="put chromosome no")
args = parser.parse_args()
chrom = str(args.chr)
pop = "AFA"

#important functions needed
def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")  
     snpanot = snpanot[(((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) |
                        ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) |
                        ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) |
                        ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) |
                        ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) |
                        ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) |
                        ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) |
                        ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))) &
                       (snpanot["rsid"].notna())]
     snpanot = snpanot.drop_duplicates(["varID"])
     return snpanot


def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) &
                           (gene_anot["gene_type"].isin(gene_types))]
     return gene_anot

def get_gene_type (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_type = gene_type.iloc[0,5]
     return gene_type

def get_gene_name (gene_anot, gene):
     gene_name = gene_anot[gene_anot["gene_id"]==gene]
     gene_name = gene_name.iloc[0,2]
     return gene_name

def get_gene_coords (gene_anot, gene):
     gene_type = gene_anot[gene_anot["gene_id"]==gene]
     gene_coord = [gene_type.iloc[0,3], gene_type.iloc[0,4]]
     return gene_coord

def get_covariates (cov_filepath):
      cov = pd.read_csv(cov_filepath, sep=" ")
      cov = cov.set_index("IID") #make IID to be the row names
      cov.index.names = [None] # remove the iid name from the row
      pc = ["PC1", "PC2", "PC3"] #a list of the PCs to retain
      cov = cov[pc]
      return cov

def get_gene_expression(gene_expression_file_name, gene_annot):
	expr_df = pd.read_csv(gene_expression_file_name, header = 0, index_col = 0, delimiter='\t')
	expr_df = expr_df.T
	inter = list(set(gene_annot['gene_id']).intersection(set(expr_df.columns)))
	#print(len(inter))
	expr_df = expr_df.loc[:, inter ]
	return expr_df

def adjust_for_covariates (expr_vec, cov_df):   
      reg = LinearRegression().fit(cov_df, expr_vec)
      ypred = reg.predict(cov_df)
      residuals = expr_vec - ypred
      residuals = scale(residuals)
      return residuals

def get_maf_filtered_genotype(genotype_file_name,  maf=0.01):
	gt_df = pd.read_csv(genotype_file_name, 'r', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

def get_cis_genotype (gt_df, snp_annot, coords, cis_window=1000000):
      snp_info = snpannot[(snpannot['pos'] >= (coords[0] - cis_window)) &
                          (snpannot['rsid'].notna()) &
                          (snpannot['pos'] <= (coords[1] + cis_window))]
      if len(snp_info) == 0:
          return 0
      else:
           gtdf_col = list(gt_df.columns)
           snpinfo_col = list(snp_info["varID"])
           intersect = snps_intersect(gtdf_col, snpinfo_col) #this function was defined earlier
           cis_gt = gt_df[intersect]
           return cis_gt

def snps_intersect(list1, list2):
     return list(set(list1) & set(list2))

# Set file paths

snp_dosage_file = "/home/pokoro/data/mesa_models/"+pop.lower()+"/"+pop.upper()+"_"+chrom+"_snp.txt"
gene_expression_file = "/home/pokoro/data/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
pc_file = "/home/pokoro/data/mesa_models/"+pop.lower()+"/"+pop.upper()+"_3_PCs.txt"
gene_annotation_file = "/home/pokoro/data/mesa_models/gencode.v18.annotation.parsed.txt"
snp_annotation_file = "/home/pokoro/data/mesa_models/"+pop.lower()+"/"+pop.upper()+"_"+chrom+"_annot.txt"


# parse the files

snpannot = get_filtered_snp_annot(snp_annotation_file)
geneannot = get_gene_annotation(gene_annotation_file, chrom)
cov = get_covariates(pc_file)
expr_df = get_gene_expression(gene_expression_file, geneannot)
genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(snp_dosage_file)


en = ElasticNet(max_iter=10000, random_state=1234)


#where to write out result
open("/home/pokoro/data/mesa_models/en_R_Python_compare/"+pop+"_en_py_chr"+chrom+
     ".txt", "w").write("gene_id"+"\t"+"gene_name"+"\t"+"chr"+"\t"+"cvr2")

#Go through all protein coding genes

for gene in genes:
    coords = get_gene_coords(geneannot, gene)
    gene_name = get_gene_name(geneannot, gene)
    expr_vec = expr_df[gene]
    
    adj_exp = adjust_for_covariates(list(expr_vec), cov)
    cis_gt = get_cis_genotype(gt_df, snpannot, coords)

    #build the model
    
    if (type(cis_gt) != int) & (cis_gt.shape[1] > 0):
         

         x = cis_gt.values
         y = adj_exp.ravel()



         #Elastic Net
         cvr2 = cross_val_score(en, x, y, scoring=r2, cv=5).mean()

         open("/home/pokoro/data/mesa_models/en_R_Python_compare/"+pop+"_en_py_chr"+
              chrom+".txt", "a").write("\n"+gene+"\t"+gene_name+"\t"+chrom+
                                       "\t"+str(cvr2))
         
