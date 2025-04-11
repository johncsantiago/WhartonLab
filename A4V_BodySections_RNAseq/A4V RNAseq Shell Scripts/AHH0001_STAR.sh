#!/bin/bash

#SBATCH -J AHH0001_STAR.sh
#SBATCH -o /gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/Output_Files/Scratch/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/Output_Files/Scratch/error_%j.out

#Submit cmd: sbatch -n 26 -t 16:00:00 --mem=64G AHH0001_STAR.sh

shelldir="/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Shell_Scripts/" 
sampledir="/gpfs/data/kwharton/azhou/ALS_model_gbb/seq_files/"
stardir="/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/"
genomedir="/gpfs/data/kwharton/jsantiago/BDGP6.32/"

module load star
module load samtools
module load python


sample=$sampledir"AHH0001"


gzip -d $sample".fastq.gz"


cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

python -m venv $sample
source ~/$sample/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"

