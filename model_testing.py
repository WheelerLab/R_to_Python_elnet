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
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
import time
from scipy import stats


#time the whole script per chromosome
chrom = 20
pop = "CAU"
#open("/home/paul/mesa_models/python_ml_models/whole_script_chr"+str(chrom)+"_timer.txt", "w").write("Chrom"+"\t"+"Time(s)"+"\n")
#t0 = time.time()

#important functions needed
def get_filtered_snp_annot (snpfilepath):
     snpanot = pd.read_csv(snpfilepath, sep="\t")
     #snpanot = snpanot[((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="A") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="A")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="G")) | ((snpanot["refAllele"]=="G") & (snpanot["effectAllele"]=="T")) | ((snpanot["refAllele"]=="T") & (snpanot["effectAllele"]=="C")) | ((snpanot["refAllele"]=="C") & (snpanot["effectAllele"]=="T"))]
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
          return NaN
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

#chrom = 21 #chromosome number. #this is removed. and initialized early at the top

afa_snp = "/home/paul/mesa_models/AFA_"+str(chrom)+"_snp.txt"
gex = "/home/paul/mesa_models/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno_Nk-10.txt"
cov_file = "/home/paul/mesa_models/AFA_3_PCs.txt"
geneanotfile = "/home/paul/mesa_models/gencode.v18.annotation.parsed.txt"
snpfilepath = "/home/paul/mesa_models/AFA_"+str(chrom)+"_annot.txt"

#test data files
test_snp = "/home/paul/mesa_models/CAU_"+str(chrom)+"_snp.txt"
test_gex = "/home/paul/mesa_models/meqtl_sorted_CAU_MESA_Epi_GEX_data_sidno_Nk-10.txt"
test_covfile = "/home/paul/mesa_models/CAU_3_PCs.txt"
test_snpfile = "/home/paul/mesa_models/CAU_"+str(chrom)+"_annot.txt"

#train functioning
snpannot = get_filtered_snp_annot(snpfilepath)
geneannot = get_gene_annotation(geneanotfile, chrom)
cov = get_covariates(cov_file)
expr_df = get_gene_expression(gex, geneannot)
genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(afa_snp, 0.01)
train_ids = list(gt_df.index)

#test functioning
test_snpannot = get_filtered_snp_annot(test_snpfile)
test_cov = get_covariates(test_covfile)
test_expr_df = get_gene_expression(test_gex, geneannot)
test_genes = list(test_expr_df.columns)
test_gt_df = get_maf_filtered_genotype(test_snp, 0.01)
test_ids = list(test_gt_df.index)
ypred_frame = pd.DataFrame() #frame to store the ypred
test_adj_exp_frame = pd.DataFrame() #frame to store test adjusted expression

#algorithms to use
rf = RandomForestRegressor(max_depth=None, random_state=1234, n_estimators=100)
svrl = SVR(kernel="linear", gamma="auto")
svr = SVR(kernel="rbf", gamma="auto")
knn = KNeighborsRegressor(n_neighbors=10, weights = "distance")
#models = [rf,svrl,svr,knn]

#text file where to write out the cv and test results
open("/home/paul/mesa_models/python_ml_models/results/"+pop+"_rf_test_chr"+str(chrom)+".txt", "w").write("gene_id"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+"pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+"spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")
open("/home/paul/mesa_models/python_ml_models/results/"+pop+"_knn_test_chr"+str(chrom)+".txt", "w").write("gene_id"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+"pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+"spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")
open("/home/paul/mesa_models/python_ml_models/results/"+pop+"_svr_linear_cv_chr"+str(chrom)+".txt", "w").write("gene_id"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+"pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+"spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")
open("/home/paul/mesa_models/python_ml_models/results/"+pop+"_svr_rbf_cv_chr"+str(chrom)+".txt", "w").write("gene_id"+"\t"+"pearson_yadj_vs_ypred (a)"+"\t"+"a_pval"+"\t"+"pearson_yobs_vs_ypred (b)"+"\t"+"b_pval"+"\t"+"spearman_yadj_vs_ypred (c)"+"\t"+"c_pval"+"\t"+"spearman_yobs_vs_ypred (d)"+"\t"+"d_pval"+"\n")

for gene in genes:
    if gene in test_genes:
        coords = get_gene_coords(geneannot, gene)
        #print(gene)
        expr_vec = expr_df[gene]#observed exp
        test_expr_vec = test_expr_df[gene]#observed exp
        #print(expr_vec)
        adj_exp = adjust_for_covariates(list(expr_vec), cov)#adjusted exp
        test_adj_exp = adjust_for_covariates(list(test_expr_vec), test_cov)#adjusted exp
        #break
        #expr_vec = expr_df[gene]
        #adj_exp = adjust_for_covariates(expr_vec, cov)
        cis_gt = get_cis_genotype(gt_df, snpannot, coords)
        test_cis_gt = get_cis_genotype(test_gt_df, test_snpannot, coords)
        gg = [gene] #just to cast the gene id to list because pandas need it to be in list before it can be used as col name

        #take the snps
        train_snps = list(cis_gt.columns)
        test_snps = list(test_cis_gt.columns)
        snp_intersect = snps_intersect(train_snps, test_snps)

        cis_gt = cis_gt[snp_intersect]
        test_cis_gt = test_cis_gt[snp_intersect]
        
        #build the model
        #adj_exp = adj_exp.values #not needed after making adj_exp a numpy array above
        cis_gt = cis_gt.values
        test_cis_gt = test_cis_gt.values
        test_yobs = test_expr_vec.values

        #prepare test_adj_exp for writing out to a file
        test_adj_exp_pd = pd.DataFrame(test_adj_exp)
        test_adj_exp_pd.columns = gg
        test_adj_exp_pd.index = test_ids
        test_adj_exp_frame = pd.concat([test_adj_exp_frame,test_adj_exp_pd], axis=1)
        
        #these steps can be shortened with a loop where the models are in a list or dictionary
        #Random Forest
        #rf_t0 = time.time()#do rf and time it
        #rf_cv = str(float(mean(cross_val_score(rf, cis_gt, adj_exp.ravel(), cv=5))))
        #rf_t1 = time.time()
        #rf_tt = str(float(rf_t1 - rf_t0))
        rf.fit(cis_gt, adj_exp.ravel())
        ypred = rf.predict(test_cis_gt)

        #prepare ypred for writing out to a file
        yrep_pd = pd.DataFrame(ypred)
        
        ypred_pd.columns = gg
        ypred_pd.index = test_ids
        ypred_frame_rf = pd.concat([ypred_frame_rf, ypred_pd], axis=1)
        
        pa = stats.pearsonr(test_adj_exp, ypred)
        pacoef = str(float(pa[0]))
        papval = str(float(pa[1]))
        pb = stats.pearsonr(test_yobs, ypred)
        pbcoef = str(float(pb[0]))
        pbpval = str(float(pb[1]))
        pc = stats.spearmanr(test_adj_exp, ypred)
        pccoef = str(float(pc[0]))
        pcpval = str(float(pc[1]))
        pd = stats.spearmanr(test_yobs, ypred)
        pdcoef = str(float(pd[0]))
        pdpval = str(float(pd[1]))
        open("/home/paul/mesa_models/python_ml_models/results/"pop+"_rf_test_chr"+str(chrom)+".txt", "a").write(gene+"\t"+rf_cv+"\t"+rf_tt+"\n")

        #SVR Linear
        svrl_t0 = time.time()#time it
        svrl_cv = str(float(mean(cross_val_score(svrl, cis_gt, adj_exp.ravel(), cv=5))))
        svrl_t1 = time.time()
        svrl_tt = str(float(svrl_t1 - svrl_t0))
        open("/home/paul/mesa_models/python_ml_models/results/svr_linear_cv_chr"+str(chrom)+".txt", "a").write(gene+"\t"+svrl_cv+"\t"+svrl_tt+"\n")
        #SVR RBF
        svr_t0 = time.time()#time it
        svr_cv = str(float(mean(cross_val_score(svr, cis_gt, adj_exp.ravel(), cv=5))))
        svr_t1 = time.time()
        svr_tt = str(float(svr_t1 - svr_t0))
        open("/home/paul/mesa_models/python_ml_models/results/svr_rbf_cv_chr"+str(chrom)+".txt", "a").write(gene+"\t"+svr_cv+"\t"+svr_tt+"\n")
        #KNN
        knn_t0 = time.time()#time it
        knn_cv = str(float(mean(cross_val_score(knn, cis_gt, adj_exp.ravel(), cv=5))))
        knn_t1 = time.time()
        knn_tt = str(float(knn_t1 - knn_t0))
        open("/home/paul/mesa_models/python_ml_models/results/knn_cv_chr"+str(chrom)+".txt", "a").write(gene+"\t"+knn_cv+"\t"+knn_tt+"\n")
        

ypred_frame_rf.to_csv("/home/paul/mesa_models/python_ml_models/results/AFA_2_"+pop+"_rf_predicted_gene_expr_chr"+str(chrom)+".txt", header=True, index=True, sep="\t")
test_adj_exp_frame.to_csv("/home/paul/mesa_models/python_ml_models/results/"+pop+"_pc_adjusted_gene_expr_chr"+str(chrom)+".txt", header=True, index=True, sep="\t")
#t1 = time.time()
#total = str(float(t1-t0))
#open("/home/paul/mesa_models/python_ml_models/whole_script_chr"+str(chrom)+"_timer.txt", "a").write(str(chrom)+"\t"+total+"\n")
#coords = get_gene_coords(geneannot, "geneID")#this is where to loop for gene id
#adj_exp = adjust_for_covariates(expr_vec, cov) #this is loop side
#cis_gt = get_cis_genotype(gt_df, snpannot, coords) #this is loop side

