#!/bin/bash


#SBATCH -J S9M_STAR.sh

#SBATCH -o /gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/Output_Files/Scratch/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/Output_Files/Scratch/error_%j.out

#Submit cmd: sbatch -n 30 -t 24:00:00 -N 2 --mem=64G S9M_STAR.sh
#/oscar/data/kwharton/jsantiago

shelldir="/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Shell_Scripts/" 
sampledir="/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/00_fastq/"
stardir="/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/STAR_Alignment/"
genomedir="/gpfs/data/kwharton/jsantiago/BDGP6.32/"

module load fastqc
module load star
module load samtools
module load python







sample="S9MA1"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MA2"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MA3"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"















sample="S9MH1"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MH2"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MH3"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"















sample="S9MT1"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MT2"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"




sample="S9MT3"
sample1=$sampledir$sample"_R1_001"
sample2=$sampledir$sample"_R2_001"

fastqc $sample1".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
fastqc $sample2".fastq.gz" -o $stardir"Output_Files/FastQC/" -t 16
gzip -d $sample1".fastq.gz"
gzip -d $sample2".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample1".fastq" $sample2".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sample
samtools index $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam"

cd ~
python -m venv RunningHTseq
source ~/RunningHTseq/bin/activate
htseq-count $stardir"Output_Files/BAM/"$sample"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/30-927564627_A4V_RNAseq/Count_Tables/"$sample"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample1".fastq"
gzip $sample2".fastq"
