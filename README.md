
# Predicting interactions of m6A-regulated RNA-binding proteins 

Sofia Martello   |   Master's Thesis Project   |   1st may 2023   |   TU Munich Campus Straubing 



### Input Data 
- RNA-seq for HEK293 cells [(Sun et al. 2018)](https://pubmed.ncbi.nlm.nih.gov/30526041/) --> [GSE122425_all.counts.293_vs_293NK.edgeR_all](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122425)
- RNA-seq for HepG2 (Wold, ENCODE) 
- ENCODE eCLIP for protein-RNA interactions (Wold, ENCODE)
- miCLIP for m6A modifications on HEK293 cells [(Linder et al. 2015)](https://pubmed.ncbi.nlm.nih.gov/26121403/) --> actually this one has been processed before importing [GSE63753_hek293.abcam.CIMS.m6A.9536.bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63753)


### Dataset Preparation 

Regarding the dataset preparation, there are three different settings to be considered: 
- baseline : does not contain m6A information
- setting A : contains only sequences which include at least one m6A site
- setting B : ratio 1:1 between the sequences containing and non-containing m6A sites 


### Setup 

- Duplicate the environment (denbi-conda_env.yml, denbi-jupyterlab_env.yml). 
- Utilize the files encoding.py, filtering.py and model.py to create a package 'src', then to be installed with `pip install -e .` 
- Organise folders as suggested :
  |-- data
  |-- docs
  |-- results
  |-- scripts
  |-- src
  |-- tests


### Workflow

The order in which the Jupyter Notebooks should be run is the following :

- Preprocessing + map_gene_ids.R
- Encoding
- Plots
- Model
- Baseline_vs_Setting_A 
- Baseline_vs_Setting_B
- Baseline_vs_Setting_A_aug
- Baseline_vs_Setting_B_aug
  
For the Baseline_vs_Setting notebooks--> possibility to input either negative-1 (unbound sequences) or negative-2 (sequences bound to ther RBPs) as class 0. 


Credits : 

LICENSE!!!!
