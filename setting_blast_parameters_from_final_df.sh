#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 8
#$ -l h_vmem=3g,virtual_free=2.9g
#$ -N set_blast_para
set -vx

cd /home/yiheng/script
python setting_blast_parameters_from_final_df.py

