import os,sys
import importlib
import src
importlib.reload(src)
import src
import src.encoding as enc
import shlex, subprocess
import pandas as pd
from pybedtools import BedTool
import numpy as np
import seaborn as sns 
import matplotlib.pyplot as plt


def random_l_r(n, k):
    k_nums = np.random.default_rng().multinomial(n, [1/k]*k, size=1)[0]
    l = k_nums[0]
    r =  k_nums[1]
    return l,r


def get_center_pos_between_rbp_and_m6a(m6a_pos, rbp_pos):
    dist = rbp_pos - m6a_pos
    return int( rbp_pos + dist / 2)


def plot_comparison_m6a_original(plot, cell_line):
    sns.set(style="darkgrid")
    plot= plot.sort_values('percentage', ascending=False)
    labels = plot['RBP'].str.split('_').str[0]#removed the cell line because I'm going to do two plots
    bar1 = plot['original length']
    bar2 = plot['modified length']
    percentage = plot['percentage']
    
    x = np.arange(len(labels))  # the label locations
    width = 0.35  # the width of the bars
    fig, ax = plt.subplots(figsize=(25, 10))
    rects1 = ax.bar(x - width/2, bar1, width, color='#48AAAD', label='Total number of sequences')
    rects2 = ax.bar(x + width/2, bar2, width, color='#F9812A', label='Number of methylated sequences')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Number of sequences')
    ax.set_title(cell_line) #here introduce a new plot in the same figure to divide cell lines
    ax.set_xticks(x, labels, rotation=90)
    ax.set_xlabel('RBP')
    ax.legend(fontsize = 'x-small', loc = 'upper left')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.margins(0,0)
    fig.suptitle('Comparison between set size before and after filtering for m6A sites\nAbsolute number of sequences containing m6A sites', fontsize=12)
    fig.tight_layout()
    return fig,ax

def create_df_statistics(path, file1,file2):
    directory1 = os.listdir(path)
    data = pd.DataFrame(columns = ['original length', 'modified length', 'RBP', 'avgseqxsite', 'avgsite', 'avgsitesxtotal'])
    for folder in directory1: #path, folder, file1,file2 > data_out
        directory2 = os.listdir(f"{path}{folder}")   
        original = 0
        modified = 0 
        avgseqxsite = 0
        avgsite = 0
        avgsitesxtotal = 0
        percentage = 0
        for file in directory2 :
            if file == file1 : 
                clipfile = BedTool(f"{path}{folder}/{file}")
                df_CLIP = clipfile.to_dataframe(disable_auto_names=True,names = [i for i in range(27)])
                result = df_CLIP.groupby([0,1,2]).count() # number of sequences in the original 
                original = len(result)
        for file in directory2 :
            if (file2 in file) & (original != 0 ): 
                clipfile = BedTool(f"{path}{folder}/{file}")
                df_CLIP = clipfile.to_dataframe(disable_auto_names=True, names = [i for i in range(27)]) 
                if df_CLIP.columns.any() : 
                    df_CLIP_ = df_CLIP[[0,1,2,3]].copy()
                    result1 = df_CLIP_.groupby([3]).count() #grouping to get the exact number of m6A sites, not of intersections 
                    avgseqxsite = (result[3].sum())/ len(result1)
                    df_CLIPa = df_CLIP[[12,13,14,15]].copy()
                    result2 = df_CLIPa.groupby([13,14,15]).count() #number of total sequences containing m6A sites
                    modified = len(result2)
                    percentage = (modified/original)*100
                    avgsite = (result2[12].sum())/len(result2)     #number of sites per total of sequences that contain m6A sites
                    avgsitesxtotal = (result2[12].sum())/original      #number of sites on the total of sequences
                else : 
                    modified = 0 
                    avgseqxsite = 0
                    avgsite = 0
                    avgsitesxtotal = 0
                    percentage = 0 
        plot_df = pd.DataFrame({'original length':  [original],#preparing the dataset for the plots
                'modified length': [modified],
                'RBP': [folder], 
                'avgseqxsite': [avgseqxsite], #average of how many sequences contain the same m6A site
                'avgsite': [avgsite], #average of m6A sites per sequence that contains at least one m6A site
                'avgsitesxtotal': [avgsitesxtotal],
                'percentage': [percentage]}) #average of m6A sites per total of sequences
        data = pd.concat([data, plot_df])
    return data

def create_csv_meth_rates(df,lower_thresh, upper_thresh, output_path ):
    plot_ds = df
    df_low = plot_ds[plot_ds.percentage < lower_thresh]
    df_mid = plot_ds [(plot_ds.percentage > lower_thresh) & (plot_ds.percentage < upper_thresh) ]
    df_high = plot_ds [plot_ds.percentage > upper_thresh]

    df_low['methylation'] = 'low'
    df_mid['methylation'] = 'medium'
    df_high['methylation'] = 'high'
    
    df_tot = pd.concat([df_low, df_mid , df_high])
    list_exclusions = ['SLBP_K562','SERBP1_K562','SBDS_K562','RPS11_K562']
    df_tot = df_tot.loc[~df_tot["RBP"].isin(list_exclusions)]
    df_tot.to_csv(output_path)
    return df_tot

def plot_meth(df, nr_nucleotides):
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(data=df,x="percentage", y="avgsite", hue="methylation")
    ax.set_yticklabels(ax.get_yticklabels(), rotation=40, ha="right")
    plt.title(f'Methylation rate - {nr_nucleotides}')
    ax.set_xlabel('Percentage of methylated sequences')
    ax.set_ylabel('Average number m6A sites per methylated sequence')
    ax.legend(loc = 'lower right')
    #plt.gca().spines['top'].set_visible(False)
    #plt.gca().spines['right'].set_visible(False)
    plt.tight_layout()
    return plt

def plot_absxpercentage(df):
    plt.figure(figsize=(6,6))
    ax = sns.scatterplot(data=df, x="percentage", y="modified length", hue="methylation")
    ax.set_yticklabels(ax.get_yticklabels())
    plt.title(f'Absolute number m6A sites related to the percentage ')
    ax.set_xlabel('Percentage of methylated sequences')
    ax.set_ylabel('Number of m6A sites')
    plt.tight_layout()
    return plt

def plot_percentage(df,cell_line):
    df= df.sort_values('percentage', ascending=False)
    labels = df['RBP'].str.split('_').str[0]
    sns.set(style="darkgrid")
    plt.figure(figsize=(25, 10))
    ax = sns.barplot(data=df, x=labels, y="percentage", errorbar = None, color = '#48AAAD' )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    plt.title(f'Percentage of methylated sequences - {cell_line}', fontsize=12)
    ax.set_ylabel('Percentage ( % )')
    ax.set_xlabel('RBP')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.margins(0,0)
    plt.tight_layout()
    return plt


def plot_avgsitextotal(df, cell_line):
    df= df.sort_values('percentage', ascending=False)
    labels = df['RBP'].str.split('_').str[0]
    plt.figure(figsize=(25, 10))
    ax = sns.barplot(data=df, y="avgsitesxtotal", x=labels, color = '#48AAAD', errorbar = None )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    plt.title(f'Average number of m6A sites per sequence - {cell_line}')
    ax.set_ylabel('Number of sites per sequence')
    ax.set_xlabel('RBP')
    plt.margins(0,0)
    plt.tight_layout()
    return plt

def plot_avgsite(df, cell_line):
    df= df.sort_values('percentage', ascending=False)
    labels = df['RBP'].str.split('_').str[0]
    plt.figure(figsize=(25, 10))
    ax = sns.barplot(data=df, y="avgsite", x=labels, color = '#48AAAD', errorbar = None )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    plt.title(f'Average number of m6A sites per sequence on the total of sequences containing at least one m6A site - {cell_line}')
    ax.set_ylabel('Number of sites per methylated sequence')
    ax.set_xlabel('RBP')
    plt.margins(0,0)
    plt.tight_layout()
    return plt

def plot_absm6a(df, cell_line):
    df= df.sort_values('percentage', ascending=False)
    plt.figure(figsize=(25, 10))
    labels = df['RBP'].str.split('_').str[0]
    ax = sns.barplot(data=df, y="modified length", x=labels, color = '#48AAAD', errorbar = None )
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    ax.set_ylabel('Number of m6A sites')
    ax.set_xlabel('RBP')
    plt.title(f'Absolute number of m6A sites per RBP - {cell_line}')
    plt.margins(0,0)
    plt.tight_layout()
    return plt