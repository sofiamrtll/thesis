
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

Now, the script to be run is *map_gene_id.R*. This script enables the lifting of the genes included in the file 'expressed_gene_ids_HEK.xlsx', from the version hg19 to the version hg38, in order to have equivalent versions throughot all files. The output file is cleaned from empty values and saved in 'inconsistenciesENSEMBL_noNaN.xlsx', which enables the production of a dictionary to be used in the Jupyter Notebook *THESIS - Preprocessing*. 

After this step, an intersection between the two sets of expressed genes in the different cell lines is identified: this is the set on which the eCLIP and miCLIP data are then going to be filtered, as well. The filtered outputs miCLIP and eCLIP are respectively saved in 'miCLIP.filt.out.bed' and 'peaks.crosslink.anno.filt.bed'. 

The next step is to produce sequences of homogenous length to be inputted in the model; the code to produce said sequences is contained in the Jupyter Notebook *THESIS - Encoding*. 
First, the file 'hg38.genome' is downloaded; it contains the lengths of human chromosomes for the genome version hg38. The RBP-binding site coordinates of the files 'peaks.crosslink.anno.filt.bed' in the respective RBP folder are now sloped to a length of 400 nucleotides and then converted into FASTA sequences using BedTools, producing respectively the files 'peaks.crosslink.anno.filt.bed.slop' and 'peaks.crosslink.anno.filt.bed.slop.fasta'. 
This fasta file contains what will later be the positive-labelled data to be inputted in the model. 
The negative-labelled sequences have been sampled exploring the binding site-neighbouring sequence space, making sure not to have any overlapping; two types of negatives were produced :

- negative-1, unbound sequences not containing RBP binding sites
- negative-2, sequences bound to other RBPs, not the one in exam 

For every RBP, a new folder is created, with organization : 
>RBP
> 
>>fold 0
>>
>>fold1
>>
>>fold2
>>
>>fold3
>>
>>fold4
>> 
>>>positive.fold-4
>>>
>>>>.bed
>>>>
>>>>.fasta
>>>>
>>>negative-1.fold-4
>>>
>>>>.bed
>>>>
>>>>.fasta
>>>>
>>>negative-2.fold-4
>>>
>>>>.bed
>>>>
>>>>.fasta
>>>>

The size of the negatives were balanced with the positives and all of the dataset was divided into 5 folds, to enable an equilibrated CrossValidation. The .fasta files are used for the dataset preparation of the baseline model, which does not contain m6A information. 

In order to include m6A data in the dataset, an additional preprocessing step is required: the filtering of the .bed files contained in the various folds for sequences containing at least one m6A site, the output is saved in the files '{file}.miclip.filt.bed.out'.  


### 3. Dataset Preparation 

Regarding the dataset preparation, there are three different settings to be considered: 
- baseline : does not contain m6A information
- setting A : contains only sequences which include at least one m6A site
- setting B : ratio 1:1 between the sequences containing and non-containing m6A sites 

To summarize, every setting has the possibility of producing two separate datasets, one containing the unbound sequences as negative-labelled elements, and the other with bound-to-other-RBPs sequences as negatives.





-------------------------

At this point there is the preparation of the dataset used for the baseline model, which won’t contain any m6A information. 

the sequences have been onehot encoded in arrays of shape (400,4) . For the integration of  the m6A sites, a new column has been added to the previous arrays, which includes the onehot encoding of the sites, producing arrays of shape (400,5). 





## Table of contents
- `code here` 
- here some **incredibly** ***important*** *text*
>here idk
- My favorite search engine is [Duck Duck Go](https://duckduckgo.com).
- See the section on [`code`](#code).
## Installation 


Credits : 

LICENSE!!!!
