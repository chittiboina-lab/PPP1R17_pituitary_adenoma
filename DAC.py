#!/usr/bin/env python
# coding: utf-8

# In[1]:


# pyDEseq2 requires:
# - numpy==1.23.0
# - anndata==0.8.0
# - pandas==1.4.3
# - scikit-learn==1.1.1
# - scipy==1.8.1
# - statsmodels==0.13.2

import pandas as pd
import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from subprocess import Popen, PIPE
from math import log
import pyensembl as ensembl
import io
ensembldb = ensembl.EnsemblRelease(release = 86)

genes = ensembldb.genes()
coding_genes = [gene for gene in genes if gene.biotype == 'protein_coding']



# In[20]:


# Setup Metadata
samples = {'P1_T':'../data/SN_MO/P1_T_outs/P1_T_outs_atac_fragments.tsv.gz',
           'P1_N':'../data/SN_MO/P1_N_outs/P1_N_outs_atac_fragments.tsv.gz',
           'P2_T':'../data/SN_MO/P2_T_outs/P2_T_outs_atac_fragments.tsv.gz',
           'P2_N':'../data/SN_MO/P2_N_outs/P2_N_outs_atac_fragments.tsv.gz',
           'P3_T':'../data/SN_MO/P3_T_outs/P3_T_outs_atac_fragments.tsv.gz',
           'P4_T':'../data/SN_MO/P4_T_outs/P4_T_outs_atac_fragments.tsv.gz'}

barcodes = {}

for sample in samples.keys():
    with open('../results/%s/%s_cort_barcodes_formatted.csv' % (sample, sample), 'r') as fin:
        barcodes[sample] = [line.strip() for line in fin]


# In[33]:


peaks = pd.read_csv('../results/ATAC/macs/cort_peaks.narrowPeak', sep='\t', header=None)
peaks = peaks[[0,1,2,3]]
peaks.columns = ['Chr', 'Start', 'End', 'PeakID']
peaks.index = peaks['PeakID']
peak_ids = peaks['PeakID']
queries = peak_ids.map(lambda x: '{}:{}-{}'.format(peaks.loc[x]['Chr'], peaks.loc[x]['Start'], peaks.loc[x]['End']))
queries = queries.rename('Query')
peak_key = pd.concat([peak_ids,queries], axis=1)
peak_key


# In[17]:


### Assemble Counts DF using row-wise apply function, adding columns. 
# Function goes row by row and uses the query built before and the tabix command to query the fragments file.
# Fragments are then filtered based on the corticotroph barcodes


def f(x):
    result = []
    for sample in samples.keys():
        process = Popen(['tabix', samples[sample], x['Query']], stdout=PIPE, universal_newlines=True)
        process = process.communicate()
        if process[0] == '':
            result.append(0)
        elif process[0] != '':
            frags_df = pd.read_csv(io.StringIO(process[0]), sep = '\t',header=None).drop([4], axis=1)
            frags_df = frags_df.loc[(frags_df[3].isin(barcodes[sample]))]
            result.append(len(frags_df))
    return result


# In[31]:


counts = pd.DataFrame(peak_key['Query'])
          
counts[list(samples.keys())] = counts.apply(f, axis=1, result_type='expand')

counts


# In[32]:


counts.to_csv('../results/ATAC/differential/peak_counts.csv')


# In[6]:


counts = pd.read_csv('../results/ATAC/differential/peak_counts.csv', index_col=0)
counts = counts.drop('Query', axis=1)
counts


# In[26]:


counts = counts.transpose()
# Filter out features with low fragment counts
genes_to_keep = counts.columns[counts.sum(axis=0) >= 10]
counts = counts[genes_to_keep]
counts


# In[44]:


gene_queries = {}
for gene in coding_genes:
    if gene.start < gene.end:
        gene_queries[gene.gene_name] = 'chr%s:%d-%d' % (gene.contig, gene.start-5000, gene.end+5000)
    elif gene.start > gene.end:
        gene_queries[gene.gene_name] = 'chr%s:%d-%d' % (gene.contig, gene.start+5000, gene.end-5000)    
gene_counts = pd.DataFrame(list(gene_queries.values()), index=gene_queries.keys(), columns=['Query'])
gene_counts[list(samples.keys())] = gene_counts.apply(f, axis=1, result_type='expand')
gene_counts


# In[45]:


gene_counts.to_csv('../results/ATAC/differential/gene_counts.csv')


# In[46]:


gene_counts = gene_counts.drop('Query', axis=1)
gene_counts


# In[47]:


gene_counts = gene_counts.transpose()
# Filter out features with low fragment counts
genes_to_keep = gene_counts.columns[gene_counts.sum(axis=0) >= 10]
gene_counts = gene_counts[genes_to_keep]
gene_counts


# In[48]:


clinical_df = pd.DataFrame(['Tumor', 'Normal', 'Tumor', 'Normal', 'Tumor', 'Tumor'], ['P1_T', 'P1_N', 'P2_T', 'P2_N', 'P3_T', 'P4_T'], ["Annotation"])
clinical_df


# In[50]:


## Peaks
dds = DeseqDataSet(
    counts=counts,
    clinical=clinical_df,
    design_factors="Annotation",
    refit_cooks=True,
    n_cpus=1,
)
dds


# In[51]:


dds.deseq2()


# In[31]:


stat_res = DeseqStats(dds, n_cpus=2)
stat_res.summary()
peak_results = stat_res.results_df
peak_results = peak_results.sort_values('pvalue', ascending=True)
peak_results.to_csv('../results/ATAC/differential/DAC_peaks.csv')
peak_results


# In[52]:


## Genes
dds = DeseqDataSet(
    counts=gene_counts,
    clinical=clinical_df,
    design_factors="Annotation",
    refit_cooks=True,
    n_cpus=1,
)
dds


# In[53]:


dds.deseq2()


# In[54]:


stat_res = DeseqStats(dds, n_cpus=2)
stat_res.summary()
gene_results = stat_res.results_df
gene_results = gene_results.sort_values('pvalue', ascending=True)
gene_results.to_csv('../results/ATAC/differential/DAC_gene.csv')
gene_results


# In[55]:


# peak_results[peak_results['pvalue'] < 0.05]#['log2FoldChange'] > 0
gene_results[gene_results['pvalue'] < 0.05]


# In[77]:


## generate peak files for homer analysis
peak_results = pd.read_csv('../results/ATAC/differential/DAC_peaks.csv')
# peak_results.index = peak_results['PeakID']
# peak_results = peak_results.drop(['PeakID'], axis=1)
peak_results


# In[78]:


peaks_annotated = pd.read_csv('../results/ATAC/macs/cort_peaks_annotated.narrowPeak', sep='\t')
peaks_annotated = peaks_annotated[['PeakID (cmd=annotatePeaks.pl results/ATAC/macs/cort_peaks.narrowPeak hg38)','Chr','Start','End','Strand','Peak Score','Annotation','Distance to TSS','Nearest Ensembl','Gene Name']]
peaks_annotated.columns = ['PeakID','Chr','Start','End','Strand','Peak Score','Annotation','Distance to TSS','Nearest Ensembl','Gene Name']
# peaks.index = peaks['PeakID']
# peaks = peaks.drop(['PeakID'], axis=1)
peaks_annotated.index = peaks_annotated['PeakID']
peaks_annotated = peaks_annotated.drop(['PeakID'], axis=1)
peaks_annotated


# In[79]:


peak_results_annotated = pd.DataFrame(peak_results)
peak_results_annotated.index = peak_results['PeakID']
peak_results_annotated[['Chr','Start','End','Strand','Peak Score','Annotation','Distance to TSS','Nearest Ensembl','Gene Name']] = peak_results_annotated.apply(lambda x: peaks_annotated.loc[x['PeakID']], axis=1, result_type='expand')
peak_results.index = peak_results['PeakID']
peak_results = peak_results.drop(['PeakID'], axis=1)
peak_results_annotated.to_csv('../results/ATAC/differential/DAC_peaks_annotated.csv')


# In[93]:


sig = peak_results_annotated[peak_results_annotated['pvalue']<0.05]
homer_up = sig[sig['log2FoldChange']>0.5][['Chr','Start','End','Strand']]
homer_up.to_csv('../results/ATAC/homer/peaks_up.homer', sep='\t')
homer_down = sig[sig['log2FoldChange']<-0.5][['Chr','Start','End','Strand']]
homer_down.to_csv('../results/ATAC/homer/peaks_down.homer', sep='\t')
homer_up

