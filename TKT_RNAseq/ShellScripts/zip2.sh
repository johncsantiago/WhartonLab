#!/bin/bash

#SBATCH -J zip2.sh
#SBATCH -o /gpfs/data/kwharton/AnswerALS/data/transcriptomics/1_fastq/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/AnswerALS/data/transcriptomics/1_fastq/error_%j.out

#Submit cmd: sbatch -n 26 -t 16:00:00 --mem=64G zip2.sh

gzip "/gpfs/data/kwharton/AnswerALS/data/transcriptomics/1_fastq/"*".fastq"

