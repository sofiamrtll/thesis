
# Predicting Interactions of m6A-Regulated RNA-Binding Proteins

**Author:** Sofia Martello  
**Master's Thesis Project | TU Munich Campus Straubing | June 30, 2023**

## Overview

This repository contains the code and data for my master's thesis project on predicting interactions of m6A-regulated RNA-binding proteins, with the aim of understanding the impact of the methylation m6A on the binding preferences of RNA-binding proteins. 

## Input Data

1. **RNA-seq for HEK293 cells (Sun et al. 2018)**
   - Dataset:  [(Sun et al. 2018)](https://pubmed.ncbi.nlm.nih.gov/30526041/) --> [GSE122425_all.counts.293_vs_293NK.edgeR_all](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122425)

2. **RNA-seq for HepG2 (Wold, ENCODE)**
   - Dataset: (Wold, ENCODE) 

3. **ENCODE eCLIP for protein-RNA interactions (Wold, ENCODE)**
   - Dataset: (Wold, ENCODE) 

4. **miCLIP for m6A modifications on HEK293 cells (Linder et al. 2015)**
   - Dataset: [(Linder et al. 2015)](https://pubmed.ncbi.nlm.nih.gov/26121403/), [GSE63753_hek293.abcam.CIMS.m6A.9536.bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63753) (processed)

## Dataset Preparation

For dataset preparation, consider three settings:

1. **Baseline:**
   - Does not contain m6A information.

2. **Setting A:**
   - Contains only sequences with at least one m6A site.

3. **Setting B:**
   - Ratio 1:1 between sequences containing and not containing m6A sites.

## Setup

1. **Duplicate the environment:**
   - Use `denbi-conda_env.yml` and `denbi-jupyterlab_env.yml`.

2. **Create a package 'src':**
   - Utilize `encoding.py`, `filtering.py`, and `model.py`.
   - Install with `pip install -e .`

3. **Organize folders as suggested:**
   ```
   |-- data
   |-- docs
   |-- results
   |-- scripts
   |-- src
   |-- tests
   ```

## Workflow

Run the Jupyter Notebooks in the following order:

1. **Preprocessing + map_gene_ids.R**
2. **Encoding**
3. **Plots**
4. **Model**
5. **Baseline_vs_Setting_A**
6. **Baseline_vs_Setting_B**
7. **Baseline_vs_Setting_A_aug**
8. **Baseline_vs_Setting_B_aug**

For the Baseline_vs_Setting notebooks, input either `negative-1` (unbound sequences) or `negative-2` (sequences bound to other RBPs) as class 0.

