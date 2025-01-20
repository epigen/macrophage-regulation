# Integrated time series analysis and high-content CRISPR screening delineates the dynamics of macrophage immune regulation

This repository contains the code and software specifications used to create the results and figures for the manuscript **[Integrated time series analysis and high-content CRISPR screening delineates the dynamics of macrophage immune regulation](https://doi.org/XXX) by Traxler, Reichl et al. (2025) *JOURNAL***

**How to cite?** 
```
Traxler, Reichl et al. (2025) *JOURNAL* XX, XX–XX, doi: https://doi.org/XXX
```

**Website:** [http://macrophage-regulation.bocklab.org/](http://macrophage-regulation.bocklab.org/)

## Instructions
The `README`'s structure follows the generated and analyzed datasets and thereby the main figures of the manuscript. Each dataset consists of mutliple analyses, performed in order, and the used software, which are linked to the respective file within the repository.
- **Code** (`src/*`) is provided as interactive notebooks, including the last outputs, and helper scripts written in `R` or `Python`
  - Notebooks are structured using `Markdown`, start with a short description of the goal, input, output, followed by loading of libraries and helper functions, configurations and data loading steps, and subsequent code for the respective analyses
  - Input paths have to be adapted at the top of the notebooks as they  might differ after data download from GEO
  - CROP-seq analyses often contain configs at the start to decide analysis parameters or whcih data to use e.g., if only mixscape (perturbed) cells or all cells should be used for the analysis
- **Software** (`envs/*.yaml`) is documented as conda environment specification files, exported using [env_export.sh](./envs/env_export.sh), in three different flavors:
  - `fromHistory`, reflects the installtion history (linked environment file)
  - `noBuild`, includes the explicit version but not build information
  - `all`, contains all package information (name, version, build)
- **Configurations** (`config/*`) contain parameters or paths used in workflows or for visualization purposes
- **Metadata** (`metadata/*`) are bulk RNA-seq and ATAC-seq annotations and metadata used in the processing, analysis and visualization
- **External resources** used are linked at the respective analysis step.
- Results with **stochastic** elements that might not be reproducible with a seed are provided in `results_stochastic/`

## Transcriptome RNA-seq time series (RNA)
Related to Figures 1, 2, S1, S2, and Table S1.
- RNA-seq preprocessing with Snakemake workflow [rna-seq-star-deseq2 (v1.2.0)](https://github.com/snakemake-workflows/rna-seq-star-deseq2)
  - Resources are downloaded automatically by the workflow according to the provided [configuration file](config/rna_rnaseq_config.yaml)
- Genome browser track visualizations with our Snakemake workflow [genome_tracks (v2.0.1)](https://github.com/epigen/genome_tracks)
  - The 12 column BED file annotation of the `mm10` genome used for the annotation of the tracks from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables) was used. Select assembly:`mm10`, track:`NCBI RefSeq`, table:`refFlat` and output format: `BED`.
- [Quality control and processing](./src/RNA_01_processing.R.ipynb) using [limma.yaml](./envs/limma.fromHistory.yaml)
- [Unsupervised analysis](./src/RNA_02_unsupervised_analysis.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- [Differential expression analysis](./src/RNA_03_DEA.R.ipynb) using [limma.yaml](./envs/limma.fromHistory.yaml)
- [Time series analysis](./src/RNA_04_time_series.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- [Enrichment analysis](./src/RNA_05_enrichment_analysis.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)

## Epigenome ATAC-seq time series (ATAC)
Related to Figures 1, 2, S1, S2, and Table S2.
- ATAC-seq preprocessing and unsupervised analysis with our Snakemake workflow [atacseq_pipeline (v0.1.0)](https://github.com/epigen/atacseq_pipeline)
  - All required resources are provided on [Zenodo](https://zenodo.org/records/6344322)
- Genome browser track visualizations with our Snakemake workflow [genome_tracks (v2.0.1)](https://github.com/epigen/genome_tracks)
  - The 12 column BED file annotation of the `mm10` genome used for the annotation of the tracks from [UCSC](https://genome.ucsc.edu/cgi-bin/hgTables) was used. Select assembly:`mm10`, track:`NCBI RefSeq`, table:`refFlat` and output format: `BED`.
- [Differential accessibility analysis](./src/ATAC_01_DEA.R.ipynb) using [limma.yaml](./envs/limma.fromHistory.yaml)
- [Time series analysis](./src/ATAC_02_time_series.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- Enrichment analysis
  - [Preparation](./src/ATAC_03_enrichment_analysis_preparation.py.ipynb) using [atac_analysis.yaml](./envs/atac_analysis.fromHistory.yaml)
  - Genomic region enrichment analysis with our Snakemake workflow [enrichment_analysis](https://github.com/epigen/enrichment_analysis/tree/7e8425c6962290a3201f8a250dd3888b23a93d7c)
  - [Aggregation](./src/ATAC_04_enrichment_analysis_aggregation.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)
- [Get](./src/ATAC_05_get_promoters.R.ipynb), using [limma.yaml](./envs/limma.fromHistory.yaml), and [quantify](./src/ATAC_06_quantify_promoters.py.ipynb) promoter regions, using [atac_analysis.yaml](./envs/atac_analysis.fromHistory.yaml), for integrative analysis (INT) of gene-promoter pairs

## Integrative analysis of RNA-seq and ATAC-seq (INT)
Related to Figures 3, S3-5, and Table S3.
- [Processing, integration and differential analysis](./src/INT_01_processing_DEA.R.ipynb) using [limma.yaml](./envs/limma.fromHistory.yaml)
- [Unsupervised analysis](./src/INT_02_unsupervised_analysis.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- [Correlation analysis](./src/INT_03_RNA_ATAC_correlation_analysis.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- [Time series analysis](./src/INT_04_time_series.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
- [Enrichment analysis](./src/INT_05_enrichment_analysis.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)
- [Transcription factor binding site motif enrichment analysis](./src/INT_06_TFBS_motif_enrichment_analysis.R.ipynb) using [rcistarget.yaml](./envs/rcistarget.fromHistory.yaml)
  - As database we used `mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather` from [Aerts lab cistarget resources](https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/)

## Proof-of-concept CROP-seq KO15 screen (KO15)
Related to Figures 4, S6-8, Table S4, and S6.
- [Processing and unsupervised analysis of all cells](./src/KO15_01_processing_unsupervised_analysis_all.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
  - For cell cycle scoring [human-mouse homologs from Mouse Genome Informatics (MGI)](https://www.informatics.jax.org/downloads/reports/index.html#homology) were used
- [Mixscape perturbation analysis of all cells simultaneously](./src/KO15_02_mixscape_all.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Mixscape perturbation analysis of each condition separately](./src/KO15_03_mixscape_conditions.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Processing and unsupervised analysis of only Mixscape selected (perturbed) cells](./src/KO15_04_processing_unsupervised_analysis_mixscape.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- Unsupervised analysis of surface protein expression with our Snakemake workflow [unsupervised_analysis (v0.2.0)](https://github.com/epigen/unsupervised_analysis)
- [Differential expression analysis within conditions between KOs and within KOs between conditions](./src/KO15_05_DEA.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Enrichment analysis](./src/KO15_06_enrichment_analysis.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)

## Upscaled CROP-seq KO150 screen (KO150)
Related to Figures 5, 6, S9-12, Table S5, and S6.
- [Processing and unsupervised analysis of all cells](./src/KO150_01_processing_unsupervised_analysis_all.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
  - For cell cycle scoring [human-mouse homologs from Mouse Genome Informatics (MGI)](https://www.informatics.jax.org/downloads/reports/index.html#homology) were used
- [Mixscape perturbation analysis of all cells simultaneously](./src/KO150_02_mixscape_all.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Mixscape perturbation analysis of each condition separately](./src/KO150_03_mixscape_conditions.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Processing and unsupervised analysis of only Mixscape selected (perturbed) cells](./src/KO150_04_processing_unsupervised_analysis_mixscape.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- Unsupervised analysis of surface protein expression with our Snakemake workflow [unsupervised_analysis (v0.2.0)](https://github.com/epigen/unsupervised_analysis)
- [Differential expression analysis within conditions between KOs and within KOs between conditions](./src/KO150_05_DEA_all.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Enrichment analysis](./src/KO150_06_enrichment_analysis.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)
- Cross-prediction similarity graph
  - [Creation with machine learning](./src/KO150_07_KO_classifier.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)
  - [Visualization](./src/KO150_08_KO_clf_plot.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
  - [Comparison with STRING](./src/KO150_09_KO_clf_stringdb.R.ipynb) using [stringdb.yaml](./envs/stringdb.fromHistory.yaml)
  - [Interpretation](./src/KO150_10_KO_clf_interpretation_analysis.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)
- [Integrative analysis of Listeria bulk results with CROP-seq KO150 results](./src/KO150_11_DEA_INT_comparison.py.ipynb) using [enrichment_analysis.yaml](./envs/enrichment_analysis.fromHistory.yaml)

## Ep300 validation experiments
Related to Figures 5 and S10.
- qPCR data
  - [qPCR standard curves for pimer efficiency determination](./Ep300_validation/Ep300_inh_qPCR_standardcurves.csv)
  - CRISPR-based Ep300 knockout untreated [plate 1](./Ep300_validation/Ep300_ko_untreated_1_1-15_untreated_RAW.csv), [plate 2](./Ep300_validation/Ep300_ko_untreated_2_1-15_untreated_RAW.csv), and [IFN-b pre-treated](./Ep300_validation/Ep300_ko_IFNb_1-15_untreated_RAW.csv) qPCR measurements
  - [Small molecule inhibtion of Ep300 untreated and IFN-b pre-treated qPCR measurements](./Ep300_validation/Ep300_inh_1-36_RAW.csv)
- [Processing and analysis of qPCR measurements](./src/Ep300_validation_analysis.R.ipynb) using [tidyverse.yaml](./envs/tidyverse.fromHistory.yaml)
- [Visualization of qPCR results](./src/Ep300_validation_figure.R.ipynb) using [tidyverse.yaml](./envs/tidyverse.fromHistory.yaml)

## Figures & Tables

### Main Figures
- [Figure 1 - RNA-seq & ATAC-seq unsupervised and differential analyses](./src/Figure_1.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Figure 2 - RNA-seq & ATAC-seq time series and enrichment analyses](./src/Figure_2.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Figure 3 - Integrative analysis of RNA-seq and ATAC-seq](./src/Figure_3.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Figures 3, 4, and 5 volcano plots](./src/Figure_3_4_5_volcanos.R.ipynb) using [enhancedVolcano.yaml](./envs/enhancedVolcano.fromHistory.yaml)
- [Figure 4 - Proof-of-concept CROP-seq KO15 screen](./src/Figure_4.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Figure 5 - Upscaled CROP-seq KO150 screen](./src/Figure_5.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Figure 6 part 1/2 - Integrative analysis across time points](./src/Figure_6_1.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Figure 6 part 2/2 - Integrative analysis with Listeria time series](./src/Figure_6_2.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)

### Supplementary Figures & Tables
- [Supplementary Figures for RNA-seq & ATAC-seq time series](./src/Supplementary_Figures_RNA_ATAC.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Supplementary Figures for RNA-seq & ATAC-seq time series volcano plots](./src/Supplementary_Figures_RNA_ATAC_volcanos.R.ipynb) using [enhancedVolcano.yaml](./envs/enhancedVolcano.fromHistory.yaml)
- [Supplementary Figures for Integrative analysis of RNA-seq and ATAC-seq](./src/Supplementary_Figures_INT.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Supplementary Figures for Proof-of-concept CROP-seq KO15 screen](./src/Supplementary_Figures_KO15.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Supplementary Figures for Proof-of-concept CROP-seq KO15 screen volcano plots](./src/Supplementary_Figures_KO15_volcanos.R.ipynb) using [enhancedVolcano.yaml](./envs/enhancedVolcano.fromHistory.yaml)
- [Supplementary Figures for Proof-of-concept CROP-seq KO15 screen Euler diagrams](./src/Supplementary_Figures_KO15_euler.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Supplementary Figures for Upscaled CROP-seq KO150 screen](./src/Supplementary_Figures_KO150.R.ipynb) using [Seurat.yaml](./envs/Seurat.fromHistory.yaml)
- [Supplementary Figures for Upscaled CROP-seq KO150 screen volcano plots](./src/Supplementary_Figures_KO150_volcanos.R.ipynb) using [enhancedVolcano.yaml](./envs/enhancedVolcano.fromHistory.yaml)
- [Supplementary Figures for Upscaled CROP-seq KO150 screen Euler diagrams](./src/Supplementary_Figures_KO150_euler.R.ipynb) using [plot_env.yaml](./envs/plot_env.fromHistory.yaml)
- [Tables S1-5](./src/Supplementary_Tables.py.ipynb) using [python.yaml](./envs/python.fromHistory.yaml)


We generalized and expanded most of these analyses to [Snakemake](https://snakemake.github.io/) workflows in an effort to augment multi-omics research by streamlining bioinformatics analyses into modules and recipes. For more details and instructions check out the project's repository here: [MrBiomics](https://github.com/epigen/MrBiomics).

