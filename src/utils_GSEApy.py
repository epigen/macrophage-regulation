#!/bin/env python

### libraries
# general
import os
import pandas as pd
import numpy as np

# visualization
import seaborn as sns
import matplotlib.pyplot as plt

# enrichment analysis
import gseapy as gp

def do_enrichment_all(gene_lists, background, databases, databases_strict, adj_pvalue, adj_pvalue_strict, dir_results, top_n=5):
    """perform enrichment analysis on given gene lists and make summary tables and summary heatmaps

    parameters:
    gene_lists = dictionary consisting of lists of genes in the correct format (recommendation: upper case and gene IDs)
    background = list of background genes in correct format (recommendation: upper case and gene IDs)
    databases = list of Enrichr database names
    databases_strict = list of Enrichr database names, where the more strict adj.p-value filter should be applied
    adj_pvalue =  adj.p-value filter
    adj_pvalue_strict = strict adj.p-value filter
    dir_result = os.path. where the results should be saved
    top_n = number of top terms of each gene_list to be included in the summary plots
    """

    # download all databases to local dict 
    db_dict = dict()
    for db in databases:
        db_dict[db]=gp.parser.gsea_gmt_parser(db, min_size=0, max_size=100000,)
    db_dict.keys()

    # perform enrichment analyses
    bg_n = len(background)

    res = dict()
    for db in db_dict.keys():
        print(db)
        res[db] = dict()

        for key in gene_lists.keys():
#             print("{} \n".format(key))

            res[db][key] = gp.enrichr(gene_list=gene_lists[key],
                             gene_sets=db_dict[db],
                             background=background,
                             organism='mouse',
                             outdir=os.path.join(dir_results, key, db),
                             top_term=25,
                             cutoff=0.05,
                             format='svg',
                             verbose=False,
                            )
            # move on if result is empty
            if res[db][key].results.shape[0]==0:
                continue

            # annotate used gene set
            res[db][key].results['Gene_set'] = db

            # separate export
            res[db][key].results.to_csv(os.path.join(dir_results, key, db, "Enrichr_{}_{}.csv".format(db, key)))

    # get all enrichment terms that are stat. significant in at least one gene list per database and top_n for plotting
    sig_terms = dict()
    plot_terms = dict()

    for db in res.keys():
        sig_terms[db]=list()
        plot_terms[db]=list()
        
        for gene_list in res[db].keys():
            # skip empty results
            if res[db][gene_list].results.shape[0]==0:
                continue

            if db in databases_strict:
                tmp_sig_terms=res[db][gene_list].results.loc[res[db][gene_list].results['Adjusted P-value']<adj_pvalue_strict,'Term'].to_list()
            else:
                tmp_sig_terms=res[db][gene_list].results.loc[res[db][gene_list].results['Adjusted P-value']<adj_pvalue,'Term'].to_list()
            
            sig_terms[db]=sig_terms[db]+tmp_sig_terms
            
            if len(tmp_sig_terms)>0:
                # select only top_n stat. sign. terms for plotting
                tmp_plot_terms = res[db][gene_list].results.loc[[x in tmp_sig_terms for x in res[db][gene_list].results['Term']],:].sort_values('P-value').loc[:,'Term'].to_list()
                plot_terms[db]=plot_terms[db] + tmp_plot_terms[:min(top_n,len(tmp_plot_terms))]
        
        
        sig_terms[db]=list(set(sig_terms[db]))
        plot_terms[db]=list(set(plot_terms[db]))
        print("{} {}".format(db, len(sig_terms[db])))

    # make summary data frames
    adj_pval=dict()
    overlap=dict()
    odds_ratio=dict()

    for db in res.keys():
        # make empty summary data frames for saving and plotting
        adj_pval[db]=pd.DataFrame(index=sig_terms[db], columns=gene_lists.keys())
        overlap[db]=pd.DataFrame(index=sig_terms[db], columns=gene_lists.keys())
        odds_ratio[db]=pd.DataFrame(index=sig_terms[db], columns=gene_lists.keys())

        for gene_list in res[db].keys():
            # skip empty results
            if res[db][gene_list].results.shape[0]==0:
                continue

            # determine sig. terms within the result
            idx_intersect = list(set(res[db][gene_list].results.set_index('Term').index).intersection(set(sig_terms[db])))

            # fill data frames
            adj_pval[db].loc[idx_intersect,gene_list]=res[db][gene_list].results.set_index('Term').loc[idx_intersect,'Adjusted P-value']
            overlap[db].loc[idx_intersect,gene_list]=res[db][gene_list].results.set_index('Term').loc[idx_intersect,'Overlap']
            odds_ratio[db].loc[idx_intersect,gene_list]=res[db][gene_list].results.set_index('Term').loc[idx_intersect,'Odds Ratio']

    # save summary data frames
    for db in res.keys():
        adj_pval[db].to_csv(os.path.join(dir_results, "summary_{}_adjpvalues.csv".format(db)))

        # convert overlaps from string to proportion by evaluation
        overlap[db][overlap[db].isna()]='0'
        overlap[db].applymap(eval).to_csv(os.path.join(dir_results, "summary_{}_overlap.csv".format(db)))

        odds_ratio[db].to_csv(os.path.join(dir_results, "summary_{}_oddsratio.csv".format(db)))

    # plot selected top_n Terms per gene_list from summary data frames as clustermaps
    for db in res.keys():
        
        # determine number of non-empty gene_list results ie columns
        gene_lists_n = sum((overlap[db].loc[plot_terms[db],:].applymap(eval)!= 0).any(axis=0))
        
        # clustermap parameters
#         width = round(len(gene_lists.keys())/4)
        width = round(gene_lists_n/3)
        height = round(len(plot_terms[db])/4)
        dendrogram_ratio = 0.01
        cbar_pos=(0.02, 0.98, 0.01, 0.05)

        ### plot adjusted p value
        # transform adj_pvalues for plotting
        mask = adj_pval[db].loc[plot_terms[db],:].isna()
        adj_pval_plot = -1*np.log10(adj_pval[db].loc[plot_terms[db],:].apply(pd.to_numeric))
        adj_pval_plot[adj_pval_plot>4]=4 # cap −log10(q)  at 4
        adj_pval_plot[mask]=0 # set NaN to 0=log10(1)  to enable clustering
        # remove columns(=gene lists) consisting only of zeros
        non_zero_cols = (adj_pval_plot != 0).any(axis=0) 
        mask=mask.loc[:,non_zero_cols]
        adj_pval_plot=adj_pval_plot.loc[:,non_zero_cols]
        # plot adj_pvalues
        tmp_map = sns.clustermap(adj_pval_plot, figsize=(width, height), mask=mask, cmap="plasma", yticklabels=True, xticklabels=True, dendrogram_ratio=dendrogram_ratio, cbar_kws={'label': '-log10(adj.p-value) \n capped at 4;  \n NaN masked & set to 0', 'shrink': 0.5}, cbar_pos=cbar_pos)
        # save fig
        tmp_map.savefig(os.path.join(dir_results,"summary_heatmap_{}_adjPvalue.svg".format(db)))

        ### plot overlap proportion
        # remove columns(=gene lists) consisting only of zeros
        non_zero_cols = (overlap[db].loc[plot_terms[db],:].applymap(eval)!= 0).any(axis=0)
        overlap_plot = overlap[db].loc[plot_terms[db],non_zero_cols].applymap(eval)
        # plot overlap
        tmp_map = sns.clustermap(overlap_plot, figsize=(width, height), cmap="plasma", yticklabels=True, xticklabels=True, dendrogram_ratio=dendrogram_ratio, cbar_kws={'label': 'Overlap', 'shrink': 0.5}, cbar_pos=cbar_pos)
        tmp_map.savefig(os.path.join(dir_results,"summary_heatmap_{}_overlap.svg".format(db)))

        ### plot odds ratio
        # transform odds ratios for plotting
        mask = odds_ratio[db].loc[plot_terms[db],:].isna()
        odds_ratio_plot = np.log10(odds_ratio[db].loc[plot_terms[db],:].apply(pd.to_numeric))
        odds_ratio_plot[mask]=0 # set NaN to 0=log10(1) to enable clustering
        # set inf to max to enable clustering
        max_or = odds_ratio_plot[odds_ratio_plot != np.inf].max()
        odds_ratio_plot.replace(np.inf,max_or,inplace=True)
        # remove columns(=gene lists) consisting only of zeros
        non_zero_cols = (odds_ratio_plot != 0).any(axis=0) 
        mask=mask.loc[:,non_zero_cols]
        odds_ratio_plot=odds_ratio_plot.loc[:,non_zero_cols]
        # plot odds ratio
        tmp_map = sns.clustermap(odds_ratio_plot, figsize=(width, height), mask=mask, cmap="plasma", yticklabels=True, xticklabels=True, dendrogram_ratio=dendrogram_ratio,  cbar_kws={'label': 'log10(Odds Ratio) \n set inf to max; NaN masked & set to 0', 'shrink': 0.5}, cbar_pos=cbar_pos)
        tmp_map.savefig(os.path.join(dir_results,"summary_heatmap_{}_oddsratio.svg".format(db)))
        
        plt.close("all")