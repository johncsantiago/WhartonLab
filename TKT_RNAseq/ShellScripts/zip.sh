#!/bin/bash

#SBATCH -J zip.sh
#SBATCH -o /gpfs/data/kwharton/azhou/ALS_model_gbb/seq_files/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/azhou/ALS_model_gbb/seq_files/error_%j.out

#Submit cmd: sbatch -n 26 -t 16:00:00 --mem=64G zip.sh

gzip "/gpfs/data/kwharton/azhou/ALS_model_gbb/seq_files/"*".fastq"

