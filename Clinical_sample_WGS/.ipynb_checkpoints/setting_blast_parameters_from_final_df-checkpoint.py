
# coding: utf-8

# This is a notebook to pull down all the information for figure 3 (ncbi hits). The design is to compair between each sample (barcode) instead of each flowcell as they sequenced the same thing. Therefore I just pulled all the ncbi hits for each flowcell for each barcode and merged them together according to the barcode.

# In[13]:

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import itertools
import argparse


# In[2]:

# use argparse to do this
#BASEDIR = args.BASEDIR
BASEDIR = '/home/yiheng/test'

# here we define the folder name of the dataframe it created by capturing the folder from the BASDIR
folder_name = os.path.basename(BASEDIR)
#column_name = folder_name.split('_')[-1]


# In[3]:

# first check if the analysis folder is there
folder_list = 'analysis  basecalled_data  scripts  tracking  workspace'.split(' ')
for x in range(0,folder_list.count('')):
    folder_list.remove('')
#fix this test
if not set(os.listdir(os.path.abspath(BASEDIR))) == set (folder_list):
    print("Something wrong with basefolder. check it please.")


# In[4]:

# get the dataframe there
dataframe = os.path.join(BASEDIR, 'analysis', 'summary_df_%s.tab' % folder_name)
sum_df = pd.read_csv(dataframe, sep='\t')


# In[5]:

# get the dataframe there
tax_linage_dataframe = os.path.join(BASEDIR, 'analysis', '%s_nttaxa.tab' % folder_name)
tax_df = pd.read_csv(tax_linage_dataframe, sep='\t')


# In[6]:

# now just have a few more filtering to make the dataframe smaller so we can process faster
#### first is to select everthing that did hit something
sum_df = sum_df[sum_df.pident_nt > 0]
#### second it needs to pass the filtering steps and has blast hit
sum_df = sum_df[(sum_df.sseqid_nt != False) & (sum_df.passes_filtering == True) & (sum_df.pc_survived == True)] 


# In[7]:

#### thirdly just add one more column of the similarity rate for the blast for change parameters of filtering
sum_df['similarity_rate'] = sum_df.length_nt/sum_df.read_length_pc*100


# In[8]:

#### drop the colums that we do not need
sum_df = sum_df.drop(columns = ['passes_filtering', 'sequence_length_template', 'mean_qscore_template', 
                                'barcode_arrangement', 'barcode_score', 'sseqid_nt', 'kit',
                                'length_nt', 'staxids_nt', 'nident_nt',
                                'variant', 'pc_survived', 'qseqid_nt', 'sacc_nt', 'scomnames_nt'])


# In[38]:

final_df = pd.merge(sum_df, tax_df, how='outer',left_on= 'read_id', right_on='read_id')


# In[39]:

final_df.staxids_nt.fillna(False, inplace=True)
final_df = final_df[final_df.staxids_nt != False]
final_df = final_df.reset_index(drop=True)


# In[45]:

clinical_df = final_df[final_df.barcode_arrangement == 'barcode07']
#clinical_df.head()


# In[ ]:

iterables = [ [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90],
              [50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90] ]

os.mkdir(os.path.join(BASEDIR, 'analysis','blast_pitent_similarity'))
for t in itertools.product(*iterables):
    clinical_df_fn = os.path.join(BASEDIR, 'analysis','blast_pitent_similarity', 'blast_pitent_similarity_%s_%s.txt' % (t[0], t[1]))
    clinical_df[(clinical_df.pident_nt > t[0]) & (clinical_df.similarity_rate > t[1])].to_csv(clinical_df_fn, sep='\t', index=None)


# In[ ]:



