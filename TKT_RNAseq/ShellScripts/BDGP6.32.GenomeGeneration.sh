#!/bin/bash

#SBATCH -J BDGP6.32.GenomeGeneration
#SBATCH -o /gpfs/data/kwharton/jsantiago/BDGP6.32/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/jsantiago/BDGP6.32/error_%j.out

#Submit cmd: sbatch -n 4 -t 16:00:00 --mem=16G BDGP6.32.GenomeGeneration.sh

module load star/2.7.3a

workingDir="/gpfs/data/kwharton/jsantiago/BDGP6.32/"

##gzip -d $workingDir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa.gz"
##gzip -d $workingDir"Drosophila_melanogaster.BDGP6.32.109.gtf.gz"

STAR --runMode genomeGenerate --genomeDir $workingDir"STAR_Index" --genomeFastaFiles $workingDir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa" --sjdbGTFfile $workingDir"Drosophila_melanogaster.BDGP6.32.109.gtf" --sjdbOverhang 149 --genomeSAindexNbases 13