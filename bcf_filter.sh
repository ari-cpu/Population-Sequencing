path=~/Documents/Master/MATLAB/mock_reads
sample=Bsub_MockReads_Mix_changeAll_bcfmp

infile=$path/$sample"_bcf.vcf"

sed -e '/[[:space:]]<\*>/d' -e 's/,<\*>//g' -e 's/..$//g' -e '/INDEL/d' -e '/#/d' $infile  > $path/$sample"_bcfilt.vcf"
sed 's/[[:space:]]/z/g' $path/$sample"_bcfilt.vcf" | cut -dz -f2,4,5,10  > $path/$sample"_cut_bcfilt.vcf"
sed -e 's/z/[[:space:]]/g' $path/$sample"_cut_bcfilt.vcf" > $path/"_test.vcf" && mv $path"_test.vcf" $path/$sample"_cut_bcfilt.vcf"
sed -e 's/,255*//g' -e 's/0://g' $path/$sample"_cut_bcfilt.vcf" > $path/$sample"_clean_bcfilt.vcf" && rm $path/$sample"_cut_bcfilt.vcf"



 



