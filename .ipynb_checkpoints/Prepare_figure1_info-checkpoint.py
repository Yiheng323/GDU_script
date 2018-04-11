
# coding: utf-8

#  This is the code to plot out the figures for the WGS data
#  By now you will need to run both YHscript_2 and YHscript_3 to perform the data analysis and construct the final dataframe
#  All the plotting will be performed in the analysis folder under the BASEDIR

# In[52]:

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse


# In[53]:

parser = argparse.ArgumentParser(description='This is a script to get information from figure one. It require to have a sum_df_DATE_FLOWCELLID.tab file in the analysis folder of folder of each run.')
parser.add_argument("BASEDIR", help="base folder, supposed to have all the sub folders processed by WGS script. The same as Indir in YH_script2. remenber DO NOT put the backslash '/' at the end!")
args = parser.parse_args()


# In[84]:

# use argparse to do this
BASEDIR = args.BASEDIR
#BASEDIR = '/home/yiheng/data/Wagga_run2'

# here we define the folder name of the dataframe it created by capturing the folder from the BASDIR
folder_name = os.path.basename(BASEDIR)


# In[85]:

# first check if the analysis folder is there
folder_list = 'analysis  basecalled_data  scripts  tracking  workspace'.split(' ')
for x in range(0,folder_list.count('')):
    folder_list.remove('')
#fix this test
if not set(os.listdir(os.path.abspath(BASEDIR))) == set (folder_list):
    print("Something wrong with basefolder. check it please.")


# In[86]:

# get the dataframe there
dataframe = os.path.join(BASEDIR, 'analysis', 'summary_df_%s.tab' % folder_name)
sum_df = pd.read_csv(dataframe, sep='\t')


# In[87]:

sum_df.columns


# In[88]:

##################################################################
################### CODE for Figure 1#############################
##################################################################


# In[89]:

########## CODE for preparing df for Figure1B

#replace the nan with False for better handling
sum_df.qseqid_rg.fillna(False, inplace=True)
sum_df.qseqid_nt.fillna(False, inplace=True)
pd.set_option('display.max_columns', None)

# define different group of reads for generating graph
ntblasthit_reads = sum_df[(sum_df.qseqid_nt != False) & (sum_df.qseqid_rg == False) & (sum_df.passes_filtering == True)]
rgblasthit_reads = sum_df[(sum_df.qseqid_nt == False) & (sum_df.qseqid_rg != False) & (sum_df.passes_filtering == True)]
noblasthit_reads = sum_df[(sum_df.qseqid_nt == False) & (sum_df.qseqid_rg == False) & (sum_df.passes_filtering == True)]

# creating the column for figure 1B, which is the read distribution proportion by hit on each database for each flowcellID
total_passed_reads = len(sum_df[(sum_df.passes_filtering == True)])

nohit_prop = len(noblasthit_reads)/total_passed_reads
nthit_prop = len(ntblasthit_reads)/total_passed_reads
rghit_prop = len(rgblasthit_reads)/total_passed_reads

# here we define the column name based on the flowcell ID captured from folder_name
column_name = folder_name.split('_')[-1]
blasthit_df = pd.DataFrame([rghit_prop, nthit_prop, nohit_prop])
blasthit_df.columns = [column_name]
blasthit_df.to_csv(r'/home/yiheng/analysis/WGS/%s_hit.tab' % column_name, header=column_name, index=None, sep='\t')


# In[90]:

blasthit_df


# In[91]:

################### CODE for preparing df for Figure1A


# In[92]:

barcode01_reads = sum_df[(sum_df.barcode_arrangement == 'barcode01') & (sum_df.passes_filtering == True)]
barcode02_reads = sum_df[(sum_df.barcode_arrangement == 'barcode02') & (sum_df.passes_filtering == True)]
barcode03_reads = sum_df[(sum_df.barcode_arrangement == 'barcode03') & (sum_df.passes_filtering == True)]
barcode04_reads = sum_df[(sum_df.barcode_arrangement == 'barcode04') & (sum_df.passes_filtering == True)]
barcode05_reads = sum_df[(sum_df.barcode_arrangement == 'barcode05') & (sum_df.passes_filtering == True)]
unclassified_reads = sum_df[(sum_df.barcode_arrangement == 'unclassified') & (sum_df.passes_filtering == True)]

total_barcoded_reads = len(barcode01_reads) + len(barcode02_reads) + len(barcode03_reads) + len(barcode04_reads) + len(barcode05_reads) + len(unclassified_reads)
barcode01_prop = len(barcode01_reads)/total_barcoded_reads
barcode02_prop = len(barcode02_reads)/total_barcoded_reads
barcode03_prop = len(barcode03_reads)/total_barcoded_reads
barcode04_prop = len(barcode04_reads)/total_barcoded_reads
barcode05_prop = len(barcode05_reads)/total_barcoded_reads
unclassified_prop = len(unclassified_reads)/total_barcoded_reads

barcode_df = pd.DataFrame([barcode01_prop, barcode02_prop, barcode03_prop, barcode04_prop, barcode05_prop, unclassified_prop])
barcode_df.columns = [column_name]
barcode_df.to_csv(r'/home/yiheng/analysis/WGS/%s_barcode.tab' % column_name, header=column_name, index=None, sep='\t')


# In[93]:

barcode_df


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



