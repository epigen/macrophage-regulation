#!/bin/env python

# libraries

# general
import os
import pandas as pd
import numpy as np
from itertools import compress
from tqdm import tqdm

# visualization
import seaborn as sns
import matplotlib.pyplot as plt

# dimensionality reduction
import umap
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# utils for enumerating string lists
from collections import defaultdict
from itertools import count
from functools import partial

# generic dim red visualization function (UMAP & PCA) for discrete and continous variables (data rows are observations)
def dimred_plot(data, annot, variables, label, results_dir=None, color_dict=None, alpha=1, centroids=True):
    
    if results_dir!=None and os.path.isfile(os.path.join(results_dir, "UMAP_{}.csv".format(label))):
        # load results
        data_umap_df=pd.read_csv(os.path.join(results_dir, "UMAP_{}.csv".format(label)), index_col=0)
    else:
        # UMAP - Dimensionality reduction via unsupervised UMAP for visualization purpose (init with PCA was horrible)
        data_umap = umap.UMAP(
            #     n_neighbors=15,
            #     min_dist=1.0,
            # init='spectral', #default
    #         densmap=True,
            n_components=2,
            random_state=42,
            metric="correlation", # tested "cosine" -> results not as good
        ).fit(data)
    #     print(data_umap.embedding_.shape)
        data_umap_df = pd.DataFrame(data_umap.embedding_, index=data.index,)
    #     data_umap_df = data_umap_df.rename_axis(("sample_name"))
        
        if results_dir!=None:
            # save results
            data_umap_df.to_csv(os.path.join(results_dir, "UMAP_{}.csv".format(label)))
    
    if results_dir!=None and os.path.isfile(os.path.join(results_dir, "PCA_{}.csv".format(label))):
        # load results
        data_pca_df=pd.read_csv(os.path.join(results_dir, "PCA_{}.csv".format(label)), index_col=0)
        explained_variance=list(pd.read_csv(os.path.join(results_dir,"PCA_{}_explained_variance.csv".format(label)), index_col=0)['0'])
    else:
        # PCA - dimensionality redcution via PCA so that 99.9% of variance is preserved
        pca_obj = PCA(0.999, random_state=42,)
        data_pca=pca_obj.fit_transform(StandardScaler().fit_transform(data))
        data_pca_df = pd.DataFrame(data_pca, index=data.index,)
        #     data_pca_df = data_pca_df.rename_axis(("sample_name"))
        explained_variance=pd.DataFrame(pca_obj.explained_variance_ratio_)
        
        if results_dir!=None:
            # save results
            data_pca_df.to_csv(os.path.join(results_dir,  "PCA_{}.csv".format(label)))
            explained_variance.to_csv(os.path.join(results_dir,  "PCA_{}_explained_variance.csv".format(label)))
        
        explained_variance=list(explained_variance.iloc[:,0])


    for dimred in ['PCA','UMAP']:
        
        # make dataframes for plotting
        if dimred=='PCA':
            data_df = data_pca_df
        if dimred=='UMAP':
            data_df = data_umap_df

        # plot data in 2D
        # sns.set(color_codes=True)
        for variable in variables:
            fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(5, 5))
            unique_vals = annot[variable].unique()
            
#             if len(unique_vals)==1:# or (sum(pd.isna(unique_vals))>0):
#                 print("{} only one value or contains NaNs".format(variable))
#                 continue

            unique_vals = unique_vals[pd.notna(unique_vals)]

            if all([isinstance(i, (str, bool, np.bool_)) for i in unique_vals]):
                print('discrete variable ', variable)
                
                if color_dict==None:
#                     colors = plt.cm.get_cmap("tab20").colors
                    cm = plt.get_cmap('gist_ncar')
                    colors=[cm(1.*i/len(unique_vals)) for i in range(len(unique_vals))]
                else:
                    colors = [color_dict[x] for x in unique_vals]

                # plot data by unique value
                for g, c in zip(unique_vals, colors):
                    c = [c]  # to remove the warnings
                    axes.scatter(
                        data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,0].values,
                        data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,1].values,
                        label=g,
                        marker=".",
                        c=c,
                        alpha=alpha,
                    )
                    
                    if centroids:
                        # centroids by mean
                        axes.scatter(
                            data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,0].mean(),
                            data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,1].mean(),
                            marker="s",
                            alpha=0.5,
                            label=str(g) + " centroid",
                            c=c,
                        )

                        axes.text(
                            data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,0].mean(),
                            data_df.loc[annot.index[annot[variable] == g].tolist(), :].iloc[:,1].mean(),
                            s=g,
                            horizontalalignment="center",
                            verticalalignment="bottom",
                            alpha=0.5,
                        )
                        
                if centroids:
                    # edit and position legend
                    handles, labels = axes.get_legend_handles_labels()
                    legend_idx = ["centroid" in label for label in labels]
                    handles = list(compress(handles, legend_idx))
                    labels = list(compress(labels, legend_idx))
                    axes.legend(handles, labels, loc="center left", bbox_to_anchor=(1.05, 0.5))

            elif all([isinstance(i, (int, float, np.int, np.float, np.int64)) for i in unique_vals]):
                print('continous variable ', variable)

                cmap = sns.cubehelix_palette(as_cmap=True)
                points = axes.scatter(
                    data_df.iloc[:, 0,],
                    data_df.iloc[:, 1,],
                    marker=".",
                    c=annot[variable],
                    s=50, 
                    cmap=cmap,
                    alpha=alpha,
                )
                fig.colorbar(points)
            else:
                print("variable type not-detected for {}".format(variable))
                continue

            # show and save figure
            x_postfix=''
            y_postfix=''
            
            if dimred=='PCA':
                x_postfix=' ({:.1f}%)'.format(100*explained_variance[0])
                y_postfix=' ({:.1f}%)'.format(100*explained_variance[1])
            
            axes.set_xlabel(dimred+" 1"+x_postfix)
            axes.set_ylabel(dimred+" 2 "+y_postfix)
            plt.title(dimred+" of " + label + " "+variable)
            plt.show()
            
            # save plot if directory is provided
            if results_dir!=None:
                fig.savefig(
                    fname=os.path.join(results_dir, dimred+"_"+label+"_"+variable+".svg"),
                    format="svg",
                    dpi=300,
                    bbox_inches="tight",
                )
