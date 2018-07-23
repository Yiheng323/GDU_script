#!/bin/bash
#$ -M yiheng.hu@anu.edu.au
#$ -m a
#$ -cwd
#$ -V
#$ -j y
#$ -pe threads 12
#$ -l h_vmem=3g,virtual_free=2.9g
#$ -N barcode07
set -vx

BASEFOLDER='/home/yiheng/test' # make sure the input file is only a copy
NAME='Hu_FAH05731_albacore202'
DATE='20171025'
BARCODE='barcode07'

cd ${BASEFOLDER}/basecalled_data
### tar -xvf ${NAME}.tar.gz
cd ${NAME}/workspace

### filter CDS

cat pass/${BARCODE}/*.fastq fail/${BARCODE}/*.fastq > ${NAME}.${BARCODE}.unlysed.fastq 
gzip ${NAME}.${BARCODE}.unlysed.fastq
gunzip -c ${NAME}.${BARCODE}.unlysed.fastq.gz | NanoLyse | gzip > ${NAME}.${BARCODE}.fastq.gz
gunzip ${NAME}.${BARCODE}.fastq.gz
rm ${NAME}.${BARCODE}.unlysed.fastq.gz

### move to the workspace folder for further manipulation
cd ${BASEFOLDER}
mkdir -p workspace/${BARCODE}
mv ${BASEFOLDER}/basecalled_data/${NAME}/workspace/${NAME}.${BARCODE}.fastq ${BASEFOLDER}/workspace/${BARCODE}
cd ${BASEFOLDER}/workspace/${BARCODE}

# do porechop to chop out adapter sequence

porechop -i ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.${BARCODE}.fastq -o ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fastq --format fastq --middle_threshold 95

# convert the fastq to fasta

sed '/^@/!d;s//>/;N' ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fastq > ${BASEFOLDER}/workspace/${BARCODE}/${NAME}.chopped.${BARCODE}.fasta

# do blastn for the fasta file

blastn -query ${NAME}.chopped.${BARCODE}.fasta -db rg -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 12 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > ${NAME}_${BARCODE}_chopped.fasta.${DATE}.rgblast_output

### supplimentray information for format6:

### qseqid means Query Seq-id
### sseqid means Subject Seq-id
### evalue means Expect value
### bitscore means Bit score
### length means Alignment length
### pident means Percentage of identical matches
### nident means Number of identical matches
### sgi means Subject GI
### sacc means Subject accession
### staxids means Subject Taxonomy ID(s), separated by a ';'
### sscinames means Subject Scientific Name(s), separated by a ';'
### scomnames means Subject Common Name(s), separated by a ';'
### sskingdoms means Subject Super Kingdom(s), separated by a ';'

# cut the seqid from rgblast_output and write it into a list for separate the fasta hit files

cut -f 1 ${NAME}_${BARCODE}_chopped.fasta.${DATE}.rgblast_output > ${NAME}.rgblast.qseqid.${BARCODE}.txt

# separate the blast hit and nohit reads from the fasta file

filterbyname.sh in=${NAME}.chopped.${BARCODE}.fasta out=${NAME}.chopped.rghityes.${BARCODE}.fasta names=${NAME}.rgblast.qseqid.${BARCODE}.txt include=t
filterbyname.sh in=${NAME}.chopped.${BARCODE}.fasta out=${NAME}.chopped.rghitno.${BARCODE}.fasta names=${NAME}.rgblast.qseqid.${BARCODE}.txt include=f

# do the second blast for the rghitno fasta file

blastn -query ${NAME}.chopped.rghitno.${BARCODE}.fasta -db nt -evalue 0.01 -outfmt '6 qseqid sseqid evalue bitscore length pident nident sgi sacc staxids sscinames scomnames sskingdoms' -show_gis -num_threads 12 | sort -k1,1 -k4,4nr | sort -u -k1,1 --merge > ${NAME}_${BARCODE}_chopped.fasta.${DATE}.ntblast_output

# cut the seqid from the ntblast_output and write it into a list for separate the fasta hit files

cut -f 1 ${NAME}_${BARCODE}_chopped.fasta.${DATE}.ntblast_output > ${NAME}.ntblast.qseqid.${BARCODE}.txt

# separate the blast hit and nohit reads from the fasta file

filterbyname.sh in=${NAME}.chopped.rghitno.${BARCODE}.fasta out=${NAME}.chopped.rghitno.nthityes.${BARCODE}.fasta names=${NAME}.ntblast.qseqid.${BARCODE}.txt include=t
filterbyname.sh in=${NAME}.chopped.rghitno.${BARCODE}.fasta out=${NAME}.chopped.rghitno.nthitno.${BARCODE}.fasta names=${NAME}.ntblast.qseqid.${BARCODE}.txt include=f

