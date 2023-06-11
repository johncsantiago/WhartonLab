#!/bin/bash

#SBATCH -J TktOEGR-F3_STAR.sh
#SBATCH -o /gpfs/data/kwharton/jsantiago/30-840029517/STAR_Alignment/Output_Files/Scratch/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/jsantiago/30-840029517/STAR_Alignment/Output_Files/Scratch/error_%j.out

#Submit cmd: sbatch -n 26 -t 16:00:00 --mem=64G TktOEGR-F3_STAR.sh

sample="TktOEGR-F3"

sampledir="/gpfs/data/kwharton/jsantiago/30-840029517/00_fastq/"
stardir="/gpfs/data/kwharton/jsantiago/30-840029517/STAR_Alignment/"
genomedir="/gpfs/data/kwharton/jsantiago/BDGP6.32/"

sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

module load fastqc/0.11.5
fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"

module load star/2.7.3a

STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
        
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"

STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample

module load samtools/1.3.1
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

module load python/3.7.4
module load htseq/0.11.1
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-840029517/Count_Tables/"$sample"_CountTable.txt"

gzip $sample1".fastq"
gzip $sample2".fastq"