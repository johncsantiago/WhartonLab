#!/bin/bash

#SBATCH -J AHH0010-19_STAR.sh
#SBATCH -o /gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/STAR_Alignment/Output_Files/Scratch/outfile_%j.out
#SBATCH -e /gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/STAR_Alignment/Output_Files/Scratch/error_%j.out

#Submit cmd: sbatch -n 26 -t 16:00:00 --mem=64G AHH0010-19_STAR.sh

shelldir="/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Shell_Scripts/" 
sampledir="/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/seq_files/"
stardir="/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/STAR_Alignment/"
genomedir="/gpfs/data/kwharton/jsantiago/BDGP6.32/"


module load star
module load samtools
module load python


##10

sampleID="AHH0010"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##11

sampleID="AHH0011"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##12

sampleID="AHH0012"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##14

sampleID="AHH0014"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##15

sampleID="AHH0015"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##16

sampleID="AHH0016"
sample=$sampledir$sampleID


gzip -d $sample".fastq.gz"


cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"





##17

sampleID="AHH0017"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##18

sampleID="AHH0018"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"







##19

sampleID="AHH0019"
sample=$sampledir$sampleID

gzip -d $sample".fastq.gz"

cd $stardir"first_pass"
STAR --genomeDir $genomedir"STAR_Index" --readFilesIn $sample".fastq" --runThreadN 20
STAR --runMode genomeGenerate --genomeDir $stardir"second_index" --genomeFastaFiles $genomedir"Drosophila_melanogaster.BDGP6.32.dna.toplevel.fa"* --sjdbFileChrStartEnd $stardir"first_pass/SJ.out.tab" --sjdbOverhang 149 --runThreadN 20 --genomeSAindexNbases 13

cd $stardir"second_pass"
STAR --genomeDir $stardir"second_index" --readFilesIn $sample".fastq" --runThreadN 20 --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $stardir"Output_Files/BAM/"$sampleID
samtools index $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam"

python -m venv $stardir"PythonHTseq/"$sampleID
source $stardir"PythonHTseq/"$sampleID/bin/activate
pip install HTSeq
htseq-count $stardir"Output_Files/BAM/"$sampleID"Aligned.sortedByCoord.out.bam" $genomedir"Drosophila_melanogaster.BDGP6.32.109.gtf" -t exon -s no -r pos -i gene_id -f bam > "/gpfs/data/kwharton/jsantiago/GBB_RNAseq_Data/Count_Tables/"$sampleID"_CountTable.txt"
deactivate

cd $shelldir
gzip $sample".fastq"
