#NanoCaller

##bed file should tab separated and bgzip and tabix and shoud be in chr1:start-end

for bam in ${Barcodes} ; do

while read Region ; do
    chr=$(echo $Region | cut -d":" -f1)
    pos=$(echo $Region | cut -d":" -f2)
    start=$(echo $pos | cut -d"-" -f1)
    end=$(echo $pos | cut -d"-" -f2)

mkdir -p VcfFiles/NanoCaller/${bam}/"${bam}$Region/"

docker run  \
    -v `pwd`:/data/ \
    -v  ${RefFolder}:/Ref/ \
    -v /data/NGC_Data/Analysis/NGC_Internal/Amplicon_RUN7/:/bed/ \
    -v `pwd`/Final_VcfFiles/NanoCaller/Raw/:/output/ \
    nanocaller NanoCaller  \
    -bam /data/BamFiles/${bam}_TrimmerLarge_sorted.bam  \
    --ref /Ref/${fasta} \
    --sequencing ont \
    --preset ont \
    --mincov 10 \
    --maxcov 1500 \
    -cpu 4 \
    -sample ${bam} \
    -o /data/VcfFiles/NanoCaller/${bam}/"${bam}$Region/" \
    -chrom ${chr} -start ${start} -end ${end} 

done < ${T_RegionFolder}${MedakaT_Region} ; done


#This option didnt included     --min_allele_freq 0.3 --ins_threshold 0.3 --del_threshold 0.3 \


#Merging vcf files

sudo chmod 777 -R VcfFiles/NanoCaller/

for bam in ${Barcodes} ; do
rm  "${bam}_all_NanoCaller_VCF_files"

while read Region ; do
if [ -f VcfFiles/NanoCaller/${bam}/"${bam}$Region"/variant_calls.final.vcf.gz ]; then
    echo VcfFiles/NanoCaller/${bam}/"${bam}$Region"/variant_calls.final.vcf.gz >> "${bam}_all_NanoCaller_VCF_files"
fi ; done < ${T_RegionFolder}${MedakaT_Region} ; done 

mkdir  -p Final_VcfFiles/NanoCaller/Raw/
sudo chmod 777 -R Final_VcfFiles/NanoCaller/

for bam in ${Barcodes} ; do
bcftools concat -f "${bam}_all_NanoCaller_VCF_files" | bcftools sort > Final_VcfFiles/NanoCaller/Raw/"${bam}_NanoCaller.vcf" ;done

for bam in ${Barcodes} ; do rm "${bam}_all_NanoCaller_VCF_files"  ; done



######QC 

mkdir -p Final_VcfFiles/NanoCaller/QC

for bam in ${Barcodes} ; do 
cat Final_VcfFiles/NanoCaller/Raw/${bam}_NanoCaller.vcf | \
SnpSift filter "( QUAL > ${QUAL} )  &( GEN[0].DP > ${DP} ) & ( GEN[0].FQ > ${AF} )" > temp.vcf

bgzip --force temp.vcf ; tabix --force -p vcf temp.vcf.gz
bcftools view temp.vcf.gz --regions-file  ${T_RegionFolder}/${T_Region} > ${bam}_phased_QC_Target.vcf 

gatk --java-options '-Xmx16g -XX:ParallelGCThreads=4' LeftAlignAndTrimVariants \
-R ${RefFolder}/${fasta} -V ${bam}_phased_QC_Target.vcf -O ${bam}_NanoCaller_QC_Target_Normalised.vcf \
--split-multi-allelics 


bgzip ${bam}_NanoCaller_QC_Target_Normalised.vcf  ; tabix -p vcf ${bam}_NanoCaller_QC_Target_Normalised.vcf.gz 
mv ${bam}_NanoCaller_QC_Target_Normalised.vcf.gz* Final_VcfFiles/NanoCaller/QC/ ; rm temp* ${bam}_phased_QC_Target.vcf ; rm *idx  

done 
