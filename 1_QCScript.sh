#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -N 1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user aleu1@uni-koeln.de 
#SBATCH --time=05:00:00
#SBATCH --array=1-4
i=$SLURM_ARRAY_TASK_ID

#Define input and output paths
myDataRaw="/scratch/ari/mockgenome/Raw"
myQCRaw="/scratch/ari/mockgenome/Raw/QC"
myDataTrim="/scratch/ari/mockgenome/Trim"
myQCTrim="/scratch/ari/mockgenome/Trim/QC"

#Define folders where software is installed
FastQCFold="/home/aleu1/software/FastQC"
TrimmFold="/home/aleu1/software/Trimmomatic-0.39"

# here you map against: .fasta
ID=$(ls -1 $myDataRaw |grep "newMixes"|grep "_25_" |grep "_1.fq" | sed -n ''$i'p'  | cut -d_ -f1,2,3,4,5) 

#ID=$(ls -1 $myDataRaw | grep "Mix13_changeAll"| grep "_1.fq" |  sed -n ''$i'p' | cut -d_ -f1,2,3,4)

#ID="Bsub_MockReads_newMixes_change1_5"

#FastQC before
cd $FastQCFold
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_1.fq" -t 8
./fastqc -o $myQCRaw/ $myDataRaw/$ID"_2.fq" -t 8

#trimming the reads
cd $TrimmFold
java -Xmx30G -Xms24G -jar trimmomatic-0.39.jar PE -threads 8 -trimlog $ID.TrimLog $myDataRaw/$ID"_1.fq" $myDataRaw/$ID"_2.fq" $myDataTrim/$ID"_1P".fastq.gz $myDataTrim/$ID"_1U".fastq.gz $myDataTrim/$ID"_2P".fastq.gz $myDataTrim/$ID"_2U".fastq.gz ILLUMINACLIP:$TrimmFold/adapters/TruSeq3-PE-2.fa:2:30:10:4 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#and checking the quality of the trimmed reads
cd $FastQCFold
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_1P.fastq.gz" -t 8
./fastqc -o $myQCTrim/ $myDataTrim/$ID"_2P.fastq.gz" -t 8
exit 0

