import numpy as np
from sklearn.svm import SVR
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean
import math
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame
import pickle
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import KNeighborsRegressor
import time
from scipy import stats
import argparse
from sklearn.linear_model import ElasticNet
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.metrics import explained_variance_score
from sklearn.metrics import r2_score
from sklearn.metrics import make_scorer#use to convert metrics to scoring callables

mse = make_scorer(mean_squared_error, greater_is_better=False)
r2 = make_scorer(r2_score, greater_is_better=True)
evs = make_scorer(explained_variance_score, greater_is_better=True)


#important functions needed
def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")  
     snpanot = snpanot[(((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))) & (snpanot["rsid"].notna())]
     snpanot = snpanot.drop_duplicates(["varID"])
     return snpanot


def get_gene_annotation (gene_anot_filepath, chrom, gene_types=["protein_coding"]):
     gene_anot = pd.read_csv(gene_anot_filepath, sep="\t")
     gene_anot = gene_anot[(gene_anot["chr"]==str(chrom)) & (gene_anot["gene_type"].isin(gene_types))]
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

def get_maf_filtered_genotype(genotype_file_name,  maf):
	gt_df = pd.read_csv(genotype_file_name, 'r', header = 0, index_col = 0,delimiter='\t')
	effect_allele_freqs = gt_df.mean(axis=1)
	effect_allele_freqs = [ x / 2 for x in effect_allele_freqs ]
	effect_allele_boolean = pd.Series([ ((x >= maf) & (x <= (1 - maf))) for x in effect_allele_freqs ]).values
	gt_df = gt_df.loc[ effect_allele_boolean ]
	gt_df = gt_df.T
	return gt_df

def get_cis_genotype (gt_df, snp_annot, coords, cis_window=1000000):
      snp_info = snpannot[(snpannot['pos'] >= (coords[0] - cis_window)) & (snpannot['rsid'].notna()) & (snpannot['pos'] <= (coords[1] + cis_window))]
      if len(snp_info) == 0:
          return 0
      else:
           gtdf_col = list(gt_df.columns)
           snpinfo_col = list(snp_info["varID"])
           intersect = snps_intersect(gtdf_col, snpinfo_col) #this function was defined earlier
           cis_gt = gt_df[intersect]
           return cis_gt

def calc_R2 (y, y_pred):
    tss = 0
    rss = 0
    for i in range(len(y)):
        tss = tss + (y[i])**2
        rss = rss + (((y[i]) - (y_pred[i]))**2)
    tss = float(tss)
    rss = float(rss)
    r2 = 1 - (rss/tss)
    return r2


def calc_corr (y, y_pred):
    num = 0
    dem1 = 0
    dem2 = 0
    for i in range(len(y)):
        num = num + ((y[i]) * (y_pred[i]))
        dem1 = dem1 + (y[i])**2
        dem2 = dem2 + (y_pred[i])**2
    num = float(num)
    dem1 = math.sqrt(float(dem1))
    dem2 = math.sqrt(float(dem2))
    rho = num/(dem1*dem2)
    return rho

def snps_intersect(list1, list2):
     return list(set(list1) & set(list2))

chrom = 1
"""
afa_snp = "Z:/data/mesa_models/split_mesa/HIS_chr1_genotype_chunk1.txt"
gex = "Z:/data/mesa_models/split_mesa/HIS_chr1_gex_chunk1.txt"
cov_file = "Z:/data/mesa_models/his/HIS_3_PCs.txt"
geneanotfile = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt"
snpfilepath = "Z:/data/mesa_models/his/HIS_1_annot.txt"

"""
pop = "AFHI"
afa_snp = "Z:/data/mesa_models/afhi/AFHI_chr1_genotype_chunk1.txt.gz"
gex = "Z:/data/mesa_models/afhi/AFHI_chr1_gex_chunk1.txt.gz"
cov_file = "Z:/data/mesa_models/afhi/PC3_AFHI_PCs_sorted.txt"
geneanotfile = "Z:/data/mesa_models/gencode.v18.annotation.parsed.txt"
snpfilepath = "Z:/data/mesa_models/afhi/AFHI.chr1.anno.txt.gz"


snpannot = get_filtered_snp_annot(snpfilepath)
geneannot = get_gene_annotation(geneanotfile, chrom)
cov = get_covariates(cov_file)
expr_df = get_gene_expression(gex, geneannot)
genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(afa_snp, 0.01)
