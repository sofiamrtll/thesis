#module for filtering
import pandas as pd 
import numpy as np 
from pybedtools import BedTool
import pybedtools
import os , sys
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import src

def explode_annotations(ann) :
    annots = pd.Series(ann.split(';'))
    annots = annots.str.split('|', expand = True)
    annots.columns = ['gene_id', 'gene_name', 'gene_type']
    annots.gene_id= annots.gene_id.str.split('.').str[0]
    return annots

def compare_with_common_genes(feature, genes_to_compare): 
    df = explode_annotations(feature) 
    intersection = set(df.gene_id).intersection(genes_to_compare) 
    if len(intersection) > 0 : 
        return True
    return False

def filter_CLIP(input_bed, set_common_genes, output_bed):
    eCLIP = BedTool(input_bed)
    df_eCLIP = eCLIP.to_dataframe(disable_auto_names=True, header= None)
    df_eCLIP['gene'] = df_eCLIP[9]
    df_eCLIP['in intersection'] =  df_eCLIP.gene.map(lambda x: compare_with_common_genes(x, genes_to_compare=set_common_genes))
    df_eclip_ok = df_eCLIP[df_eCLIP['in intersection'] == True]
    df_eclip_ok = df_eclip_ok.drop(columns = ['gene', 'in intersection'])
    #df_eclip_ok.to_csv(output_bed, sep='\t', header=None, index=False)
    return df_eclip_ok, len(df_eCLIP) , len(df_eclip_ok)
    
    
def plot_original_modified_old(plot, cell_line) : 
    sns.set(style="darkgrid")
    labels = plot['RBP'].str.split('_').str[0]#removed the cell line because I'm going to do two plots
    bar1 = plot['original length']
    bar2 = plot['modified length']
    percentage = plot['percentage']

    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars

    fig, ax = plt.subplots(figsize=(5, 15))
    rects1 = ax.barh(x - width/2, bar1, width, color='#FA8072', label='Unfiltered set')
    rects2 = ax.barh(x + width/2, bar2, width, color='#FFA500', label='Filtered set')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_xlabel('#genes')
    ax.set_title(cell_line) #here introduce a new plot in the same figure to divide cell lines
    ax.set_yticks(x, labels)
    ax.legend()

    ax.bar_label(rects1,labels = percentage, padding=3)
    #ax.bar_label(rects2, labels = None ,padding=3)

    ax.spines[['top', 'right']].set_visible(False)
    fig.suptitle('Comparison between unfiltered and filtered set sizes', fontsize=16)
    fig.tight_layout()
    
    return fig, ax
def plot_original_modified(plot, cell_line):
    sns.set(style="darkgrid")
    plot= plot.sort_values('RBP')
    labels = plot['RBP'].str.split('_').str[0]#removed the cell line because I'm going to do two plots
    bar1 = plot['original length']
    bar2 = plot['modified length']
    percentage = plot['percentage']
    
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots(figsize=(25, 10))
    rects1 = ax.bar(x - width/2, bar1, width, color='#48AAAD', label='Total number of sequences')
    rects2 = ax.bar(x + width/2, bar2, width, color='#F9812A', label='Number of expressed sequences')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of sequences')
    ax.set_title(cell_line) #here introduce a new plot in the same figure to divide cell lines
    ax.set_xticks(x, labels, rotation=90)
    ax.set_xlabel('RBP')
    ax.legend(fontsize = 'x-small', loc = 'upper left')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.margins(0,0)
    fig.suptitle('Comparison between set size before and after filtering for expressed genes', fontsize=12)
    fig.tight_layout()
    return fig,ax