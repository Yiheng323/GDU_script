#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 8
#$ -l h_vmem=3g,virtual_free=2.9g
#$ -N blast
set -vx
#uses a local nt database downloaded 20170802 to blast all contigs against this database
WORKDIR=~/data/20170617_FAH05432
for fasta in ~/data/20170617_FAH05432/*.fasta
do blastn -query $fasta -db nt -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 8 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > $fasta.20170821_NCBI.nt.ncbiblast_output
done
