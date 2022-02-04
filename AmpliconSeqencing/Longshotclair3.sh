##Data folder
TrimmerLarge="python3 /data/Software_Resources/ProwlerTrimmer/TrimmerLarge.py"
Mosedepth_plots='python3 /data/Software_Resources/NGC_scripts/plot-dist.py'
MosedepthTargetComparison='python3 /data/Software_Resources/NGC_scripts/MosedepthTargetComparison.py'
TargetRegionBed_toPrecentage='python3 /data/Software_Resources/NGC_scripts/MosedepthTargetRegionBed_toPrecentage.py'

Annovar='/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/'
AnnovarDb="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/"

Data='/data/NGC_Data/Nanopore_Output/2022/Ampliocn_Run8_lib1_20222/Ampliocn_Run8_lib1_20222/20220202_1239_X3_FAR36853_9fa80ceb/fastq_pass/'
RefFolder="/data/Software_Resources/ReferenceFiles/Homo_sapiens/GATK_hg38_Resources/"
fasta='Homo_sapiens_assembly38.fasta'

T_RegionFolder='/data/NGC_Data/Analysis/CDFD_Projects/2022/Dr_Ashwin_Ampliocn_Run8_lib1_20222/'
T_Region="Ampliocn_Run8_Region.bed"
MedakaT_Region="Medaka_Ampliocn_Run8_Region.bed"


DP=10
AF=0.299
QUAL=10
Qscore=10


Barcodes="barcode08"




conda activate NnaoporeQC
#Merge fASTQ ILES
mkdir -p FastQfiles/

for barcod in ${Barcodes} ; do
zcat ${Data}${barcod}/*gz > FastQfiles/${barcod}.fq ; done

#Step1: Prealignment QC 

    #A) QC checking of Raw FastQ data using Fatqc
mkdir -p QC/FastQC/Raw_Fastq_QC/
for barcod in ${Barcodes} ; do
fastqc -t 35 -o QC/FastQC/Raw_Fastq_QC/ FastQfiles/${barcod}.fq ; done 

mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/


#Adapter Removal and Low quality reads and Bases
mkdir  -p CleanedFastQ/PoreChop/

for barcod in ${Barcodes} ; do
porechop -t 15 \
-i FastQfiles/${barcod}.fq \
-o CleanedFastQ/PoreChop/"${barcod}_PoreChop.fq" --discard_middle ; done


mkdir -p CleanedFastQ/TrimmerLarge/
for barcod in ${Barcodes} ; do 
${TrimmerLarge} \
--file=CleanedFastQ/PoreChop/"${barcod}_PoreChop.fq" \
--windowsize=50 \
--trimmode=D \
--qscore=${Qscore} ;done


for barcod in ${Barcodes} ; do 
mv CleanedFastQ/PoreChop/"${barcod}_"Pore*fastq CleanedFastQ/TrimmerLarge/"${barcod}_"Pore_TrimmerLarge.fastq 
mv CleanedFastQ/PoreChop/"${barcod}_"Pore*.csv CleanedFastQ/TrimmerLarge/ ; done


#A) QC checking of QC passed FastQ data using Fatqc
mkdir -p QC/FastQC/TrimmerLarge_Fastq_QC/
for barcod in ${Barcodes} ; do
fastqc -t 35 -o QC/FastQC/TrimmerLarge_Fastq_QC/ CleanedFastQ/TrimmerLarge/"${barcod}_"Pore_TrimmerLarge.fastq ;done 

mkdir -p MultiQC/FastQC/TrimmerLargeFastq_QC
multiqc -o MultiQC/FastQC/TrimmerLargeFastq_QC/ QC/FastQC/TrimmerLarge_Fastq_QC/


##Alignment of reads
mkdir -p BamFiles 

for fastq in ${Barcodes} ; do
minimap2 -t 10 \
    -ax map-ont \
    ${RefFolder}${fasta} CleanedFastQ/TrimmerLarge/"${fastq}_"Pore_TrimmerLarge.fastq | \
    samtools view -bS | samtools sort -o BamFiles/${fastq}_TrimmerLarge_sorted.bam ;done ;wait 


parallel  samtools index ::: BamFiles/*.bam
#Alignmet summary ; pomoxis is used for this; https://github.com/nanoporetech/pomoxis
mkdir -p QC/AlignmentSummary

for bam in ${Barcodes} ; do
stats_from_bam BamFiles/${bam}_TrimmerLarge_sorted.bam \
    > QC/AlignmentSummary/${bam}_TrimmerLarge_sorted.stats & done ;wait 


#C) Mos depth global and target region and gene exon
mkdir -p QC/AlignmentQC/Mosdepth/TargetRegion/

for bam in ${Barcodes} ; do
mosdepth -F 1284 --fast-mode --by ${T_RegionFolder}${T_Region} QC/AlignmentQC/Mosdepth/TargetRegion/${bam}_TargetRegion BamFiles/${bam}_TrimmerLarge_sorted.bam \
--quantize 0:1:4:10:20:50:100:200:500:600 \
--thresholds 1,10,20,30,50,100,200,500,600 

mkdir -p QC/AlignmentQC/flagstat
samtools flagstat BamFiles/${bam}_TrimmerLarge_sorted.bam > QC/AlignmentQC/flagstat/${bam}_flagstat.txt 

${Mosedepth_plots}  QC/AlignmentQC/Mosdepth/TargetRegion/*.dist.txt --output QC/AlignmentQC/Mosdepth/TargetRegion/${bam}_MosdepthCoverage.html 
done
 


mkdir -p MultiQC/AlignmentQC/Mosdepth/TargetRegion/
multiqc -o MultiQC/AlignmentQC/Mosdepth/TargetRegion/ QC/AlignmentQC/Mosdepth/TargetRegion/

mkdir -p MultiQC/AlignmentQC/flagstat/
multiqc -o MultiQC/AlignmentQC/flagstat/ QC/AlignmentQC/flagstat


${MosedepthTargetComparison} \
	-Mfolder QC/AlignmentQC/Mosdepth/TargetRegion/ \
	-Outputfolder  QC/AlignmentQC/Mosdepth/TargetRegion/

${TargetRegionBed_toPrecentage} \
	-Mfolder QC/AlignmentQC/Mosdepth/TargetRegion/



#-----------------------------------------------------longshot-------------------------------------------------------------------
#longshot ; reads with at least 30x coverage

mkdir -p Final_VcfFiles/longshot/Raw

for bam in ${Barcodes} ; do
longshot \
--force_overwrite \
--sample_id ${bam} \
--bam BamFiles/"${bam}_TrimmerLarge_sorted.bam" \
--ref ${RefFolder}${fasta} \
--out Final_VcfFiles/longshot/Raw/"${bam}_longshot.vcf"; done

#--min_alt_frac 0.3 --min_cov 10

#Adding total depth and Allele fraction to info fields
for bam in ${Barcodes} ; do
python /data/Software_Resources/NGC_scripts/longshot_add_DP_AF_Info_tovcf.py \
-VcfFile Final_VcfFiles/longshot/Raw/"${bam}_longshot.vcf" ;done



######QC 

mkdir -p Final_VcfFiles/longshot/QC

for bam in ${Barcodes} ; do 
rm Final_VcfFiles/longshot/Raw/"${bam}"_longshot.vcf

cat Final_VcfFiles/longshot/Raw/${bam}_longshotcorrected.vcf | \
SnpSift filter "( QUAL > ${QUAL} )  &(C_TDP > ${DP} ) & ( C_AF > ${AF} )" > temp.vcf

bgzip --force temp.vcf ; tabix --force -p vcf temp.vcf.gz
bcftools view temp.vcf.gz --regions-file  ${T_RegionFolder}/${T_Region} > ${bam}_phased_QC_Target.vcf 

gatk --java-options '-Xmx16g -XX:ParallelGCThreads=4' LeftAlignAndTrimVariants \
-R ${RefFolder}/${fasta} -V ${bam}_phased_QC_Target.vcf -O ${bam}_longshot_QC_Target_Normalised.vcf \
--split-multi-allelics 

bgzip ${bam}_longshot_QC_Target_Normalised.vcf  ; tabix -p vcf ${bam}_longshot_QC_Target_Normalised.vcf.gz 
mv ${bam}_longshot_QC_Target_Normalised.vcf.gz* Final_VcfFiles/longshot/QC/ ; rm temp* ${bam}_phased_QC_Target.vcf ; rm *idx  

done 



###ANNOTATION OF VCF FILE USING aNNOVAR

for bam in ${Barcodes} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/longshot/QC/${bam}_longshot_QC_Target_Normalised.vcf.gz \
		-outfile ${bam}.avinput \
		-includeinfo -withzyg ;done

for SAMPLENAME in ${Barcodes} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME}_longshot_QC \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 


##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes} ; do 
	python ${Annovar}AnnovarFiltering.py  --AnnovarFile ${SAMPLENAME}_longshot_QC.hg38_multianno.txt
done 

#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes} ; do 
	find="Otherinfo4"
	replace=$(zgrep "#CHROM"  Final_VcfFiles/longshot/QC/${SAMPLENAME}_longshot_QC_Target_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_QC.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_QC.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(zgrep "#CHROM" Final_VcfFiles/longshot/QC/${SAMPLENAME}_longshot_QC_Target_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_QC.hg38_multianno_DesiredColumns_Filtered.tsv 
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_QC.hg38_multianno_DesiredColumns_Filtered.tsv

done

for SAMPLENAME in ${Barcodes} ; do 
python /data/Software_Resources/NGC_scripts/longshot_multianno_DesiredColumns_Filtered.py \
	 --AnnovarOutputFile ${SAMPLENAME}_longshot_QC.hg38_multiannoDesiredColumns.tsv  

python /data/Software_Resources/NGC_scripts/longshot_multianno_DesiredColumns_Filtered.py \
	         --AnnovarOutputFile ${SAMPLENAME}_longshot_QC.hg38_multianno_DesiredColumns_Filtered.tsv
rm ${SAMPLENAME}.avinput ${SAMPLENAME}_longshot_QC.hg38_multianno.txt

done




#############Raw data
for bam in ${Barcodes} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/longshot/Raw/${bam}_longshotcorrected.vcf \
		-outfile ${bam}.avinput \
		-includeinfo -withzyg ;done

for SAMPLENAME in ${Barcodes} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME}_longshot_Raw \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 



##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes} ; do 
	python ${Annovar}AnnovarFiltering.py  --AnnovarFile ${SAMPLENAME}_longshot_Raw.hg38_multianno.txt &
done 


#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes} ; do 
	find="Otherinfo4"
	replace=$(grep "#CHROM"  Final_VcfFiles/longshot/Raw/${bam}_longshotcorrected.vcf)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_Raw.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_Raw.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(grep "#CHROM" Final_VcfFiles/longshot/Raw/${bam}_longshotcorrected.vcf)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_Raw.hg38_multianno_DesiredColumns_Filtered.tsv  
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_Raw.hg38_multianno_DesiredColumns_Filtered.tsv 
done



for SAMPLENAME in ${Barcodes} ; do 
python /data/Software_Resources/NGC_scripts/longshot_multianno_DesiredColumns_Filtered.py \
	 --AnnovarOutputFile ${SAMPLENAME}_longshot_Raw.hg38_multiannoDesiredColumns.tsv 

python /data/Software_Resources/NGC_scripts/longshot_multianno_DesiredColumns_Filtered.py \
	         --AnnovarOutputFile ${SAMPLENAME}_longshot_Raw.hg38_multianno_DesiredColumns_Filtered.tsv

rm ${SAMPLENAME}.avinput ${SAMPLENAME}_longshot_Raw.hg38_multianno.txt

done



#######-----------------------------------------------------Clair3---------------------------------------------------------------------------


#Clair3

cp ${RefFolder}/${fasta} `pwd`/BamFiles/ 
cp ${RefFolder}/${fasta}.fai `pwd`/BamFiles/ 
cp ${T_RegionFolder}${T_Region} `pwd`/BamFiles/TargetRegion.bed 


for bam in ${Barcodes} ; do
mkdir -p VcfFiles/Clair3/${bam}/

INPUT_DIR=`pwd`/BamFiles/      
OUTPUT_DIR="`pwd`/VcfFiles/Clair3/${bam}/"     
THREADS="5"           
MODEL_NAME="r941_prom_sup_g5014"         

docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn=${INPUT_DIR}/${bam}_TrimmerLarge_sorted.bam --ref_fn=${INPUT_DIR}/${fasta} \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/opt/models/${MODEL_NAME}" \
  --enable_phasing \
  --remove_intermediate_dir \
  --bed_fn=${INPUT_DIR}/TargetRegion.bed \
  --sample_name=${bam} \
  --output=${OUTPUT_DIR}  ; done          

#These options didnt work   --snp_min_af=0.3 --indel_min_af=0.3 --qual=10 \


for bam in ${Barcodes} ; do
mkdir -p Final_VcfFiles/Clair3/Raw/
cp "`pwd`/VcfFiles/Clair3/${bam}/"phased_merge_output.vcf.gz  Final_VcfFiles/Clair3/Raw/${bam}_phased_merge_output.vcf.gz 
cp `pwd`/VcfFiles/Clair3/${bam}/run_clair3.log  Final_VcfFiles/Clair3/Raw/${bam}_run_clair3.log ; done 



mkdir -p Final_VcfFiles/Clair3/QC

for bam in ${Barcodes} ; do 
zcat Final_VcfFiles/Clair3/Raw/${bam}_phased_merge_output.vcf.gz | \
SnpSift filter "( QUAL > ${QUAL} )  & ( GEN[0].DP > ${DP} ) & ( GEN[0].AF > ${AF} )" > temp.vcf

bgzip --force temp.vcf ; tabix --force -p vcf temp.vcf.gz
bcftools view temp.vcf.gz --regions-file   ${T_RegionFolder}${T_Region} > ${bam}_phased_QC_Target.vcf ; 

gatk --java-options '-Xmx16g -XX:ParallelGCThreads=4' LeftAlignAndTrimVariants \
-R ${RefFolder}/${fasta} -V ${bam}_phased_QC_Target.vcf -O ${bam}_clair3_QC_Target_Normalised.vcf \
--split-multi-allelics 

bcftools sort ${bam}_clair3_QC_Target_Normalised.vcf > output.vcf ; mv output.vcf  ${bam}_clair3_QC_Target_Normalised.vcf

bgzip ${bam}_clair3_QC_Target_Normalised.vcf  ; tabix -p vcf ${bam}_clair3_QC_Target_Normalised.vcf.gz 
mv ${bam}_clair3_QC_Target_Normalised.vcf.gz* Final_VcfFiles/Clair3/QC/ ; rm temp* ${bam}_phased_QC_Target.vcf ; rm *idx  

done 

## Annotation


###ANNOTATION OF VCF FILE USING aNNOVAR
for bam in ${Barcodes} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/Clair3/QC/${bam}_clair3_QC_Target_Normalised.vcf.gz \
		-outfile ${bam}.avinput \
		-includeinfo -withzyg 
done


for SAMPLENAME in ${Barcodes} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME} \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 


##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes} ; do 
	python ${Annovar}AnnovarFiltering.py  --AnnovarFile ${SAMPLENAME}.hg38_multianno.txt &
done 

#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes} ; do 
	find="Otherinfo4"
	replace=$(zgrep "#CHROM"  Final_VcfFiles/Clair3/QC/${Barcodes}_clair3_QC_Target_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(zgrep "#CHROM" Final_VcfFiles/Clair3/QC/${Barcodes}_clair3_QC_Target_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv 
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv 
done

for SAMPLENAME in ${Barcodes} ; do 
python /data/Software_Resources/NGC_scripts/Clair3_barcode08.hg38_multianno_DesiredColumns_Filtered.py \
	 --AnnovarOutputFile ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv 

python /data/Software_Resources/NGC_scripts/Clair3_barcode08.hg38_multianno_DesiredColumns_Filtered.py \
	         --AnnovarOutputFile ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv

mv ${SAMPLENAME}.hg38_multiannoDesiredColumns.tsv  ${SAMPLENAME}_Clair3_QC.hg38_multiannoDesiredColumns.tsv 
mv ${SAMPLENAME}.hg38_multianno_DesiredColumns_Filtered.tsv  ${SAMPLENAME}_Clair3_QC.hg38_multianno_DesiredColumns_Filtered.tsv
rm ${SAMPLENAME}.avinput ${SAMPLENAME}.hg38_multianno.txt

done


#rAW VCF Annotation
for bam in ${Barcodes} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/Clair3/Raw/${bam}_phased_merge_output.vcf.gz \
		-outfile ${bam}.avinput \
		-includeinfo -withzyg 
done


for SAMPLENAME in ${Barcodes} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME}_Raw_Clair3_ \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 


##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes} ; do 
	python ${Annovar}AnnovarFiltering.py  --AnnovarFile ${SAMPLENAME}_Raw_Clair3_.hg38_multianno.txt
done 


#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes} ; do 
	find="Otherinfo4"
	replace=$(zgrep "#CHROM"  Final_VcfFiles/Clair3/Raw/${SAMPLENAME}_phased_merge_output.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_Raw_Clair3_.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_Raw_Clair3_.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(zgrep "#CHROM" Final_VcfFiles/Clair3/Raw/${SAMPLENAME}_phased_merge_output.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_Raw_Clair3_.hg38_multianno_DesiredColumns_Filtered.tsv 
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_Raw_Clair3_.hg38_multianno_DesiredColumns_Filtered.tsv 
done


for SAMPLENAME in ${Barcodes} ; do 
python /data/Software_Resources/NGC_scripts/Clair3_RawData_hg38_multianno_DesiredColumns_Filtered.py \
	 --AnnovarOutputFile ${SAMPLENAME}_Raw_Clair3_.hg38_multiannoDesiredColumns.tsv 

python /data/Software_Resources/NGC_scripts/Clair3_RawData_hg38_multianno_DesiredColumns_Filtered.py \
	         --AnnovarOutputFile  ${SAMPLENAME}_Raw_Clair3_.hg38_multianno_DesiredColumns_Filtered.tsv

mv ${SAMPLENAME}_Raw_Clair3_.hg38_multiannoDesiredColumns.tsv  ${SAMPLENAME}_Raw_Clair3_hg38_multiannoDesiredColumns.tsv
mv ${SAMPLENAME}_Raw_Clair3_.hg38_multianno_DesiredColumns_Filtered.tsv  ${SAMPLENAME}_Raw_Clair3_hg38_multianno_DesiredColumns_Filtered.tsv

rm ${SAMPLENAME}_Raw_Clair3_.hg38_multianno.txt ${SAMPLENAME}.avinput
done


