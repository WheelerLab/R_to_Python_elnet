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
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("chr", action="store", help="put chromosome no")
args = parser.parse_args() #1-8
chrom = args.chr
chrom = str(chrom)
pop = "afa"

#time the whole script per chromosome

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

def get_maf_filtered_genotype(genotype_file_name,  maf): #the input file must have column names
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

#chrom = 21 #chromosome number. #this is removed. and initialized early at the top
folder = "all"
tr_pop = "ALL"
#train data files
afa_snp = "/home/pokoro/data/mesa_models/"+folder+"/whole_genotypes/"+tr_pop+".chr"+chrom+".genotype.txt.gz"
gex = "/home/pokoro/data/mesa_models/"+folder+"/"+tr_pop+"_PF10.txt.gz"
cov_file = "/home/pokoro/data/mesa_models/"+folder+"/PC3_"+tr_pop+"_PCs_sorted.txt"
geneanotfile = "/home/pokoro/data/mesa_models/gencode.v18.annotation.parsed.txt"
snpfilepath = "/home/pokoro/data/mesa_models/"+folder+"/"+tr_pop+".chr"+chrom+".anno.txt.gz"

#test data files
test_snp = "/home/pokoro/data/lauren_mesa/ml_dosages/"+pop+"/chr"+chrom+".gz"
test_annot = "/home/pokoro/data/lauren_mesa/snp_annotation/"+pop+"/chr"+chrom+"_annot.txt"

#train functioning
snpannot = get_filtered_snp_annot(snpfilepath)
geneannot = get_gene_annotation(geneanotfile, chrom)
expr_df = get_gene_expression(gex, geneannot) #this had to created early to avoid empty df downstream
annot_geneid = geneannot["gene_id"]#remove decimal from gene_id
annot_geneid = list(annot_geneid)
agid = []
for i in annot_geneid:
	agid.append(i[0:(i.find("."))])
geneannot["gene_id"] = agid #replace with non decimal gene_id

cov = get_covariates(cov_file)

genes = list(expr_df.columns)
gt_df = get_maf_filtered_genotype(afa_snp, 0.01)
train_ids = list(gt_df.index)
train_g = [] #where to store the non decimal gene_id
for i in genes:
	train_g.append(i[0:(i.find("."))])
expr_df.columns = train_g #use the non decimal gene_id to rename the expr_df col
genes = list(expr_df.columns) #take out the new non decimal gene_id
#adj_exp_frame = pd.DataFrame()


#test functioning
test_snpannot = get_filtered_snp_annot(test_annot)
test_gt_df = get_maf_filtered_genotype(test_snp, 0.01)
test_ids = list(test_gt_df.index)


#frame to store the ypred and test adjusted expression
ypred_frame_rf = pd.DataFrame()
ypred_frame_svr = pd.DataFrame()
ypred_frame_knn = pd.DataFrame()

#test_adj_exp_frame = pd.DataFrame()

#read in the grid search best result files and take the params to fit the model
rf_grid = pd.read_csv("/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/"+tr_pop+"_best_grid_rf_chr"+chrom+"_full.txt", sep="\t")
knn_grid = pd.read_csv("/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/"+tr_pop+"_best_grid_knn_chr"+chrom+"_full.txt", sep="\t")
svr_grid = pd.read_csv("/home/pokoro/data/mesa_models/python_ml_models/merged_chunk_results/"+tr_pop+"_best_grid_svr_chr"+chrom+"_full.txt", sep="\t")

#create file with header to write out expression to file immediately
open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf.txt", "a").write("gene_id")
for i in range(len(test_ids)):
     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf.txt", "a").write("\t" + str(test_ids[i]))

open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn.txt", "a").write("gene_id")
for i in range(len(test_ids)):
     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn.txt", "a").write("\t" + str(test_ids[i]))

open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr.txt", "a").write("gene_id")
for i in range(len(test_ids)):
     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr.txt", "a").write("\t" + str(test_ids[i]))


#algorithms to use

for gene in genes:
    coords = get_gene_coords(geneannot, gene)
    gene_name = get_gene_name(geneannot, gene)
        
    #print(gene)
    expr_vec = expr_df[gene]#observed exp
        
    #print(expr_vec)
    adj_exp = adjust_for_covariates(list(expr_vec), cov)#adjusted exp
        
        
        
    cis_gt = get_cis_genotype(gt_df, snpannot, coords)
    test_cis_gt = get_cis_genotype(test_gt_df, test_snpannot, coords)

    if (type(cis_gt) != int) & (type(test_cis_gt) != int):#just to be sure the cis genotype is not empty
        gg = [gene] #just to cast the gene id to list because pandas need it to be in list before it can be used as col name

        #take the snps
        train_snps = list(cis_gt.columns)
        test_snps = list(test_cis_gt.columns)
        snp_intersect = snps_intersect(train_snps, test_snps)

        cis_gt = cis_gt[snp_intersect]
        test_cis_gt = test_cis_gt[snp_intersect]

        if (cis_gt.shape[1] > 0) & (test_cis_gt.shape[1] > 0): #make sure that the cis_gt is not empty
            #build the model
                  
            cis_gt = cis_gt.values
            test_cis_gt = test_cis_gt.values
                  
                  
            #Random Forest
            gnlist = list(rf_grid['Gene_Name'])
            f = gene_name in gnlist
            if f != False: #Just to be sure that the gene exist in the RF best grid dataframe
                n_tree = rf_grid[rf_grid['Gene_Name']==gene_name].iloc[0,3]
                rf = RandomForestRegressor(random_state=1234, n_estimators=n_tree)
                rf.fit(cis_gt, adj_exp.ravel())
                ypred = rf.predict(test_cis_gt)

                #write out ypred quickly
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf.txt", "a").write("\n")
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf.txt", "a").write(str(gene))
                for j in range(len(ypred)):
                     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf.txt", "a").write("\t"+str(ypred[j]))

                #prepare ypred for writing out to a file
                ypred_pd = pd.DataFrame(ypred)
                       
                ypred_pd.columns = gg
                ypred_pd.index = test_ids
                ypred_frame_rf = pd.concat([ypred_frame_rf, ypred_pd], axis=1, sort=True)
                       
                      
            #Support Vector Machine
            gnlist = list(svr_grid['Gene_Name'])
            f = gene_name in gnlist
            if f != False: #Just to be sure that the gene exist in the SVR best grid dataframe
                kernel = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,3]
                degree = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,4]
                c = svr_grid[svr_grid['Gene_Name']==gene_name].iloc[0,5]
                svr = SVR(gamma="scale", kernel=kernel, degree=degree, C=c)
                svr.fit(cis_gt, adj_exp.ravel())
                ypred = svr.predict(test_cis_gt)

                #write out ypred quickly
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr.txt", "a").write("\n")
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr.txt", "a").write(str(gene))
                for j in range(len(ypred)):
                     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr.txt", "a").write("\t"+str(ypred[j]))
       
                #prepare ypred for writing out to a file
                yprep_pd = pd.DataFrame(ypred)
                       
                ypred_pd.columns = gg
                ypred_pd.index = test_ids
                ypred_frame_svr = pd.concat([ypred_frame_svr, ypred_pd], axis=1, sort=True)
                       
                       
                  
            #K Nearest Neighbour
            gnlist = list(knn_grid['Gene_Name'])
            f = gene_name in gnlist
            if f != False: #Just to be sure that the gene exist in the KNN best grid dataframe
                k = knn_grid[knn_grid['Gene_Name']==gene_name].iloc[0,3]
                weight = knn_grid[knn_grid['Gene_Name']==gene_name].iloc[0,4]
                knn = KNeighborsRegressor(n_neighbors=k, weights = weight)
                knn.fit(cis_gt, adj_exp.ravel())
                ypred = knn.predict(test_cis_gt)

                #write out ypred quickly
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn.txt", "a").write("\n")
                open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn.txt", "a").write(str(gene))
                for j in range(len(ypred)):
                     open("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn.txt", "a").write("\t"+str(ypred[j]))
                       
                #prepare ypred for writing out to a file
                yprep_pd = pd.DataFrame(ypred)
                       
                ypred_pd.columns = gg
                ypred_pd.index = test_ids
                ypred_frame_knn = pd.concat([ypred_frame_knn, ypred_pd], axis=1, sort=True)
                       
                       
        

ypred_frame_rf.to_csv("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_rf_bestorder.txt", header=True, index=True, sep="\t")
ypred_frame_svr.to_csv("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_svr_bestorder.txt", header=True, index=True, sep="\t")
ypred_frame_knn.to_csv("/home/pokoro/data/twas_mesa/ml_pred/"+pop+"_chr"+chrom+"_knn_bestorder.txt", header=True, index=True, sep="\t")
