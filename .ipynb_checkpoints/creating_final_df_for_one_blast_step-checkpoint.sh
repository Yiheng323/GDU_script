#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 8
#$ -l h_vmem=5g,virtual_free=5g
#$ -N final_df
set -vx

cd /home/yiheng/script
python creating_final_df_one_blast.py /home/yiheng/test


