
# coding: utf-8

# In[14]:

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import subprocess as sub


# In[21]:

parser = argparse.ArgumentParser(description='''
This is a script to extract reads from certain species out from a fasta file with its blast output file. It will need the following files and one program installed:
    1. intall bbmap, make filterbyname.sh excutable
    2. nt taxa file that contain read id, barcode, and taxa linages
    3. input fasta, read ids in the header of each reads

it will separate the reads that hits a desired species/genus/family/order/class and put it into the input folder and name as species.fasta or genus.fasta
By the way, it will generate a tmp.txt file in the input directory when excute this scrip. Do not worry too much since it will be deleted after it is finished.
''')


# In[22]:

parser.add_argument("indir", help="input fasta file directory, contain all the reads that went though blast analysis. MUST have at least one '/', DO NOT add the back slash'/' at the end!")
#parser.add_argument("outdir", help="the directory for output fasta file, no need to name the fasta file.")
parser.add_argument("taxa", help="the taxa file path. It will contain the linages of all the species and their read ids. MUST have at least one '/'.")
parser.add_argument("--species", help="the name of species that you want to extract, if species, put a underscore '_' in between the genus name and species name")
parser.add_argument("--genus", help="the name of genus that you want to extract")
parser.add_argument("--family", help="the name of family that you want to extract")
parser.add_argument("--order", help="the name of order that you want to extract")


args = parser.parse_args()


# In[ ]:

# get the taxa dataframe there
BASEDIR = args.indir
OUTDIR = os.path.dirname(BASEDIR)
dataframe_path = args.taxa
taxa_df = pd.read_csv(dataframe_path, sep='\t')


# In[8]:

species_name = args.species
genus = args.genus
family = args.family
order = args.order


# In[10]:

filter_command = r'filterbyname.sh in=%s out=%s/%s.fasta names=%s/tmp.txt include=t'
if species_name is not None:
    genus_name = species_name.split('_')[-2]
    species_last_name = species_name.split('_')[-1]
    # make the readids list
    species_readids = taxa_df[(taxa_df.species.str.contains(species_last_name)) & (taxa_df.genus.str.contains(genus_name))]
    species_readids['read_id'].to_csv(os.path.join(OUTDIR + '/tmp.txt'), index=None)
    # now filter out reads using filterbyname.sh script
    filter_command_stderr = sub.check_output(filter_command % (BASEDIR, OUTDIR, species_name, OUTDIR), shell=True, stderr=sub.STDOUT)
    
elif genus is not None:
    genus_readids = taxa_df[taxa_df.genus.str.contains(genus)]
    genus_readids['read_id'].to_csv(os.path.join(OUTDIR + '/tmp.txt'), index=None)
    filter_command_stderr = sub.check_output(filter_command % (BASEDIR, OUTDIR, genus, OUTDIR), shell=True, stderr=sub.STDOUT)
    
elif family is not None:
    family_readids = taxa_df[taxa_df.family.str.contains(family)]
    family_readids['read_id'].to_csv(os.path.join(OUTDIR + '/tmp.txt'), index=None)
    filter_command_stderr = sub.check_output(filter_command % (BASEDIR, OUTDIR, family, OUTDIR), shell=True, stderr=sub.STDOUT)
    
elif order is not None:
    order_readids = taxa_df[taxa_df.order.str.contains(order)]
    order_readids['read_id'].to_csv(os.path.join(OUTDIR + '/tmp.txt'), index=None)
    filter_command_stderr = sub.check_output(filter_command % (BASEDIR, OUTDIR, order, OUTDIR), shell=True, stderr=sub.STDOUT)
    
else:
    print('no name detected. check the parameters')



# In[12]:

# last just remove the tmp file
os.remove(os.path.join(OUTDIR + '/tmp.txt'))


# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:




# In[ ]:



