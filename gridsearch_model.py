import numpy as np
from sklearn.svm import SVR
#from sklearn.model_selection import train_test_split
import pandas as pd
from sklearn.model_selection import cross_val_score
from statistics import mean
import math
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import scale
from pandas import DataFrame
import pickle
from sklearn.ensemble import RandomForestRegressor
#from sklearn.svm import SVR
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


parser = argparse.ArgumentParser()
parser.add_argument("chr", action="store", help="put chromosome no")
args = parser.parse_args()
chrom = args.chr
chrom = str(chrom)

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


afa_snp = "/home/paul/mesa_models/AFA_"+chrom+"_snp.txt"
gex = "/home/paul/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
cov_file = "/home/paul/mesa_models/AFA_3_PCs.txt"
geneanotfile = "/home/paul/mesa_models/gencode.v18.annotation.parsed.txt"
snpfilepath = "/home/paul/mesa_models/AFA_"+chrom+"_annot.txt"


snpannot = get_filtered_snp_annot(snpfilepath)
geneannot = get_gene_annotation(geneanotfile, chrom)
cov = get_covariates(cov_file)
expr_df = get_gene_expression(gex, geneannot)
genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(afa_snp, 0.01)

#algorithms to use
rf = RandomForestRegressor()
n_estimators = [int(i) for i in range(50,301,50)]
rf_grid = {"n_estimators": n_estimators}
rfgs = GridSearchCV(rf, rf_grid, cv=5, iid=False, scoring=r2)

svr = SVR(gamma="scale")
kernel = ["linear", "poly", "rbf", "sigmoid"]
degree = [2, 3, 4, 5, 6, 7]
C = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0, 1.5]
svr_grid = {"kernel": kernel,
            "degree": degree, "C": C}
svrgs = GridSearchCV(svr, svr_grid, cv=5, iid=False, scoring=r2)

knn = KNeighborsRegressor()
n_neighbors = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
weights = ["uniform", "distance"]
p = [1, 2, 3]
knn_grid = {"n_neighbors": n_neighbors,
            "weights": weights, "p": p}
knngs = GridSearchCV(knn, knn_grid, cv=5, iid=False, scoring=r2)

#text file where to write out the cv results
open("/home/paul/mesa_models/python_ml_models/results/grid_rf_cv_chr"+chrom+".txt", "w").write("Gene_ID"+"\t"+"Gene_Name"+"\t"+"CV_R2"+"\t"+"n_estimators"+"\n")
open("/home/paul/mesa_models/python_ml_models/results/grid_knn_cv_chr"+chrom+".txt", "w").write("Gene_ID"+"\t"+"Gene_Name"+"\t"+"CV_R2"+"\t"+"n_neigbors"+"\t"+"weights"+"\t"+"p"+"\n")
open("/home/paul/mesa_models/python_ml_models/results/grid_svr_cv_chr"+chrom+".txt", "w").write("Gene_ID"+"\t"+"Gene_Name"+"\t"+"CV_R2"+"\t"+"kernel"+"\t"+"degree"+"\t"+"C"+"\n")

for gene in genes:
    coords = get_gene_coords(geneannot, gene)
    gene_name = get_gene_name(geneannot, gene)
    expr_vec = expr_df[gene]
    
    adj_exp = adjust_for_covariates(list(expr_vec), cov)
    cis_gt = get_cis_genotype(gt_df, snpannot, coords)

    #build the model
    
    if (type(cis_gt) != int) & (cis_gt.shape[1] > 0):
         

         cis_gt = cis_gt.values
         
         #Random Forest
         
         rfgs.fit(cis_gt, adj_exp.ravel())
         rf_cv = str(rfgs.best_score_)
         n = str(rfgs.best_params_["n_estimators"])
         open("/home/paul/mesa_models/python_ml_models/results/grid_rf_cv_chr"+chrom+".txt", "w").write(gene+"\t"+gene_name+"\t"+rf_cv+"\t"+n+"\n")
         
         #SVR
         svrgs.fit(cis_gt, adj_exp.ravel())
         svr_cv = str(svrgs.best_score_)
         svr_kernel = str(svrgs.best_params_["kernel"])
         svr_degree = str(svrgs.best_params_["degree"])
         svr_c = str(svrgs.best_params_["C"])
         open("/home/paul/mesa_models/python_ml_models/results/grid_svr_cv_chr"+chrom+".txt", "w").write(gene+"\t"+gene_name+"\t"+svr_cv+"\t"+svr_kernel+"\t"+svr_degree+"\t"+svr_c+"\n")

         #KNN
         knngs.fit(cis_gt, adj_exp.ravel())
         knn_cv = str(knngs.best_score_)
         knn_n = str(knngs.best_params_["n_neighbors"])
         knn_w = str(knngs.best_params_["weights"])
         knn_p = str(knngs.best_params_["p"])
         open("/home/paul/mesa_models/python_ml_models/results/grid_knn_cv_chr"+chrom+".txt", "w").write(gene+"\t"+gene_name+"\t"+knn_cv+"\t"+knn_n+"\t"+knn_w+"\t"+knn_p+"\n")
         

