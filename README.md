
# Predicting interactions of m6A-regulated RNA-binding proteins 

Sofia Martello   |   Master's Thesis Project   |   1st may 2023   |   TU Munich Campus Straubing 

### Background & Motivation

RNA-binding proteins have essential functions in regulating various mRNA processing mechanisms. Since regulation is frequently performed by directly binding to particular transcript regions, identifying these locations is critical. 
RNA-binding proteins (RBPs) have a wide range of roles in human diseases like cancer, neurodegenerative and metabolic diseases. These roles include export and localization of transcripts, post-transcriptional modification, alternative splicing, and translation. To understand how RNA binding proteins work in the body and in disease, it is essential to identify their binding sites [(Horlacher et al. 2023)](https://doi.org/10.1101/2023.02.14.528560). 
Many deep learning-based models have lately been designed to explore RNA-RBPs interactions, for example DeepCLIP [(Gronning et al. 2020)](https://academic.oup.com/nar/article/48/13/7099/5859960) or DeepRAM [(Trabelsi et al. 2019)](https://academic.oup.com/bioinformatics/article/35/14/i269/5529112). However, most of them don't accurately represent high-fidelity interactions in vivo, mainly because they make predictions based solely on RNA sequence, failing to account for highly dynamic variables that affect the cell environment and can also impact binding. 


For example, RNA modifications help to this dynamism, and particularly N(6)-methyladenosine (m6A), which impacts adenosines and is the least rare of these methylations; it has been demonstrated to block or promote protein access to binding sites by changing the RNA structure. In particular, methylation or demethylation generates the so-called "m6A switch", which results in local conformational changes in the RNA strand that impact the relative accessibility of miRNA and RBP binding sites [(Das Mandal & Ray 2021)](https://doi.org/10.1016/j.ygeno.2020.12.027). 
In this project, we'll concentrate on m6A modifications and the development of a deep learning model that integrates m6A and RNA sequence data in order to better identify the binding preferences of RBPs in vivo and understand the consequences of m6A modifications on the regulatory code that controls binding. 



### Input Data 

Let’s now explore the datasets that were used :

- RNA-seq for HEK293 cells [(Sun et al. 2018)](https://pubmed.ncbi.nlm.nih.gov/30526041/) --> [GSE122425_all.counts.293_vs_293NK.edgeR_all](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122425)
- RNA-seq for HepG2 (Wold, ENCODE) 
- ENCODE eCLIP for protein-RNA interactions (Wold, ENCODE)
- miCLIP for m6A modifications on HEK293 cells [(Linder et al. 2015)](https://pubmed.ncbi.nlm.nih.gov/26121403/) --> actually this one has been processed before importing [GSE63753_hek293.abcam.CIMS.m6A.9536.bed](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63753)


RNA-seq are gene expression datasets and two of them were used: one for HEK293, is a cell line isolated from the kidney tissue of a human embryo, and one for HepG2 cell line, which holds high proliferation rates and performs many hepatic functions. 

The eCLIP dataset contains the coordinates for the binding sites of 223 proteins on the RNA; the study was conducted on HEK562 liver cells and HepG2 cells. The protocol eCLIP takes advantage of the tendency of reverse transcriptase to detach at the cross-linked nucleotide, which produces cDNAs with a 5′ end that maps to the first nucleotide next to the cross-linking site and this allows the identification of it at nucleotide-level resolution [(Hafner et al. 2021)](https://www.nature.com/articles/s43586-021-00018-1). It’s important to mention the fact that a good number of CLIP reads reflect transient RBP-RNA interactions, so false positives are commonly included in these datasets.  

The miCLIP dataset contains the coordinates for the m6A sites in the genome, on a study conducted on HEK293 cells. m6A residues are mapped by generating signature mutations using m6A-specific antibodies and UV cross-linking. The reverse transcription of RNA which is UV-cross-linked to certain m6A antibodies results in a very specific pattern of mutations or truncations in the cDNA, with precise identification of m6A residues [(Linder et al. 2015)](https://pubmed.ncbi.nlm.nih.gov/26121403/). 

## Workflow 

### 1. Setup 

Packages to be installed in the environment, ....  


### 2. Preprocessing 

The first Jupyter Notebook to be run is *THESIS - Preprocessing*. 
The RNA-seq input files are being filtered in order to retain the expressed genes: HEK293 is in a different genome version, so the filtered output is saved in 'expressed_gene_ids_HEK.xlsx'. 

Now, the script to be run is *map_gene_id.R*. This script enables the lifting of the genes included in the file 'expressed_gene_ids_HEK.xlsx', from the version hg17 to the version hg19, in order to have equivalent versions throughot all files. The output file is cleaned from empty values and saved in 'inconsistenciesENSEMBL_noNaN.xlsx', which enables the production of a dictionary to be used in the Jupyter Notebook *THESIS - Preprocessing*. 

After this step, an intersection between the two sets of expressed genes in the different cell lines is identified: this is the set on which the eCLIP and miCLIP data are then going to be filtered, as well. 

The filtered outputs miCLIP and eCLIP are respectively saved in 'miCLIP.filt.bed' and 'peaks.crosslink.anno.filt.bed', for later use. 





The first step of the preprocessing was the filtering of the RNA-seq for expressed genes , and once found the intersection between the two cell lines in exam, the datasets eCLIP and miCLIP were filtered as well on the base of this common set of  expressed genes. The intersection in red is the dataset that is going to be utilized. 

Talking about numbers, here we see how numerous the datasets are, the original and filtered  : the intersection between the two RNA-seq datsets is roughly 21 000 genes, the miCLIP is the smallest in comparison with the others and it’s 9000 sites; the eCLIp one is a bit more variable but usually each RBP-set contains 20 000 binding sites.  


### Dataset Preparation 

At this point there is the preparation of the dataset used for the baseline model, which won’t contain any m6A information. The eCLIP binding site coordinates, for every protein,  have been slopped to a standard length of 400 nucleotides and then converted into fasta sequences. Afterwards the sequences have been onehot encoded in arrays of shape (400,4) . For the integration of  the m6A sites, a new column has been added to the previous arrays, which includes the onehot encoding of the sites, producing arrays of shape (400,5). 
This is the production of the positive dataset. Regarding the negative dataset, different options were available: a set of unbound sequences and another set of sequences bound to other proteins ; these sets were produced investigating the neighbouring space from the sequences belonging to the positive set, making sure not to have any overlapping . The size of the negatives were balanced with the positives and all of the dataset was divided into 5 folders to make sure that the cross validation always contained a ratio between positives and negatives that made sense. 



In this slide we can see the different settings of datasets that are used ( so far I have tested the Baseline and Setting A ) .  For the baseline setting there is the positive, which is the bound sequence, and then the two negatives: one is the unbound sequences and one is the bound to other RBPs. So basically from this setting we can produce two datasets to be inputted into the model, one is the combination of positive and negative 1 and the other is positive and negative 2. Going ahead, setting A is similar to the basic setting but it just investigates sequences that contain m6A sites, at least one. These are the settings that I have tested so far but there is another option which would be the setting B, in which there is a ratio 1:1 between the sequences containing and non-containing m6A sites. 




## Table of contents
- `code here` 
- here some **incredibly** ***important*** *text*
>here idk
- My favorite search engine is [Duck Duck Go](https://duckduckgo.com).
- See the section on [`code`](#code).
## Installation 


Credits : 

LICENSE!!!!
