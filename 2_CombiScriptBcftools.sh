#!/bin/bash -l
#SBATCH -J Combi
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH -N 1
#SBATCH --account=AG-Maier
#SBATCH --mail-type=END 
#SBATCH --mail-user aleu1@uni-koeln.de  #!
#SBATCH --time=48:00:00
#SBATCH --array=1-4 #!

#Define input and output paths (non-existing directories have to be created before executing this script)
myDataTrim="/scratch/ari/mockgenome/Trim"           # where to find the trimmed sequencing reads 
myDictPath="/projects/ag-advol/dictionaries/BsubNC_000964wt" # where to find the dictionary #!
myDataPath="/scratch/ari/mockgenome/mapped"      # where to save the results #!

#Define folders where software is installed  
bwaFold="/home/aleu1/software/bwa-0.7.17" #!
samFold="/home/aleu1/software/samtools-1.12" #!
bcfFold="/home/aleu1/software/bcftools-1.12" #!
bedFold="/home/aleu1/software/bedtools2/bin" #!

#Define variables 
i=$SLURM_ARRAY_TASK_ID 
dict="Bs166NCe" #!
#ID=$(ls -1 $myDataTrim | grep "_1P.fastq.gz" | sed -n ''$i'p' | cut -d_ -f1,2) #!
ID=$(ls -1 $myDataTrim | grep "_1P.fastq.gz"| grep "newMixes" | grep "_25_" | sed -n ''$i'p' | cut -d_ -f1,2,3,4,5) #!
#ID="Bsub_MockReads_Mut_change"$i
#ID="Bsub_MockReads_Mix4_changeAll_01"

IDout=$ID #!

################################## no futher input necessary #######################
#Mapping the reads on the reference dictionary
cd $bwaFold
./bwa mem -t 8 -B 1 -O 1 $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout".sam"
#hard
#./bwa mem -t 8 -B 100 -O 100 $myDictPath/$dict".fasta" $myDataTrim/$ID"_1P.fastq.gz" $myDataTrim/$ID"_2P.fastq.gz" > $myDataPath/$IDout".sam"

#Sort_Script.sh - sorting the mapped reads for faster following steps
cd $samFold
./samtools sort $myDataPath/$IDout".sam" --threads 8 --reference $myDictPath/$dict.fasta -o $myDataPath/$IDout"_sort.bam"
./samtools index -b $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_sort.bam.bai"



#Mpileup - Information for every position on the dictionary is piled up from the information of every single read
cd $bcfFold
./bcftools mpileup -e 10 -F 0.00001 -h 80 -L 10000 -o 20 -a FORMAT/AD -d 8000 -f $myDictPath/$dict.fasta $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_bcfmp_bcf.vcf"

#Variant calling - Gives out only variants
cd $bcfFold
./bcftools call -vc $myDataPath/$IDout"_bcfmp_bcf.vcf" > $myDataPath/$IDout"_bcfmp_bcfcall.vcf" 

#filter bcf.vcf for further processing with matlab: delete positions with no variants and INDELs
sed -e '/[[:space:]]<\*>/d' -e 's/,<\*>//g' -e 's/..$//g' -e '/INDEL/d' -e '/#/d' $myDataPath/$IDout"_bcfmp_bcf.vcf"  > $myDataPath/$IDout"_bcfmp_bcfilt.vcf"

#Coverage
module unload gnu/4.4.7
module load gnu/7.4.0
cd $bedFold
./genomeCoverageBed -d -ibam $myDataPath/$IDout"_sort.bam" > $myDataPath/$IDout"_coverage.txt"


exit 0

                       
