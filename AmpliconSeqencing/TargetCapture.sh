
Data='/data/NGC_Data/Nanopore_Output/2022/AMPLICON_RUN10_7LIB_160322/AMPLICON_RUN10_7LIB_160322/20220316_1304_X2_FAR39699_a67b0d22/fastq_pass/'
RefFolder="/data/Software_Resources/ReferenceFiles/Homo_sapiens/GATK_hg38_Resources/"
fasta='Homo_sapiens_assembly38.fasta'

T_RegionFolder='/data/NGC_Data/Analysis/CDFD_Projects/2022/Dr_Ashwin/Amplicon_Run10/'
T_Region="Ampliocn_Run8_Region.bed"
MedakaT_Region="Medaka_Ampliocn_Run8_Region.bed"
OutputDir=$PWD

DP=10
AF=0.299
QUAL=10
Qscore=10
THREADS=10

Barcodes="barcode60 barcode61 barcode62 barcode63 barcode64 barcode65 barcode66"
Barcodes=("barcode60" "barcode61" "barcode62" "barcode63" "barcode64" "barcode65" "barcode66")



##Data folder
TrimmerLarge="python3 /data/Software_Resources/ProwlerTrimmer/TrimmerLarge.py"
Mosedepth_plots='python3 /data/Software_Resources/NGC_scripts/plot-dist.py'
MosedepthTargetComparison='python3 /data/Software_Resources/NGC_scripts/MosedepthTargetComparison.py'
TargetRegionBed_toPrecentage='python3 /data/Software_Resources/NGC_scripts/MosedepthTargetRegionBed_toPrecentage.py'
longshot_add_DP_AF="python /data/Software_Resources/NGC_scripts/longshot_add_DP_AF_Info_tovcf.py"
longshot_multianno="python /data/Software_Resources/NGC_scripts/longshot_multianno_DesiredColumns_Filtered.py"
AnnovarFiltering="python ${Annovar}AnnovarFiltering.py"
Clair3_IFO_FIELDS_ADDITION="python /data/Software_Resources/NGC_scripts/Clair3_IFO_FIELDS_ADDITION.py"
Longshot_clair_AnnovarOutputMerg="python /data/Software_Resources/NGC_scripts/Longshot_clair_AnnovarOutputMerging.py"
Amplicon_sequencing_demultiplex="/data/Software_Resources/NGC_scripts/Amplicon_sequencing_demultiplex_March21.py"
Amplicon_Variant_SampleSummary="python /data/Software_Resources/NGC_scripts/Amplicon_sequencing_demultiplex_Variant_SampleSummary_March21.py"

VT="/data/Software_Resources/vt/vt"
Annovar='/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/'
AnnovarDb="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/"





conda activate NnaoporeQC
#Merge fASTQ ILES
mkdir -p FastQfiles/

for barcod in ${Barcodes[@]} ; do
zcat ${Data}${barcod}/*gz > FastQfiles/${barcod}.fq ; done

#Step1: Prealignment QC 

for Samplename in ${Barcodes[@]} ; do
mkdir -p ${OutputDir}/QC/LongQC/Raw/
docker run  \
-v ${OutputDir}:/input \
-v ${OutputDir}/QC/LongQC/Raw/:/output \
longqc sampleqc -x ont-rapid -m 2 -p 5 -o /output/${Samplename}/ /input/FastQfiles/${Samplename}.fq 

done

#A) QC checking of Raw FastQ data using Fatqc
mkdir -p QC/FastQC/Raw_Fastq_QC/
for barcod in ${Barcodes[@]} ; do
fastqc -t 35 -o QC/FastQC/Raw_Fastq_QC/ FastQfiles/${barcod}.fq ; done 

mkdir -p MultiQC/FastQC/RawFastq_QC
multiqc -o MultiQC/FastQC/RawFastq_QC/ QC/FastQC/Raw_Fastq_QC/


#Adapter Removal and Low quality reads and Bases
mkdir  -p CleanedFastQ/PoreChop/

for barcod in ${Barcodes[@]} ; do
porechop -t 15 \
-i FastQfiles/${barcod}.fq \
-o CleanedFastQ/PoreChop/"${barcod}_PoreChop.fq" --discard_middle ; done


mkdir -p CleanedFastQ/TrimmerLarge/
for barcod in ${Barcodes[@]} ; do 
${TrimmerLarge} \
--file=CleanedFastQ/PoreChop/"${barcod}_PoreChop.fq" \
--windowsize=50 \
--trimmode=D \
--qscore=${Qscore} ;done


for barcod in ${Barcodes[@]} ; do 
mv CleanedFastQ/PoreChop/"${barcod}_"Pore*fastq CleanedFastQ/TrimmerLarge/"${barcod}_"Pore_TrimmerLarge.fastq 
mv CleanedFastQ/PoreChop/"${barcod}_"Pore*.csv CleanedFastQ/TrimmerLarge/ ; done


#Post Trimming QC
for Samplename in ${Barcodes[@]} ; do
mkdir -p ${OutputDir}/QC/LongQC/TrimmerLarge/
docker run  \
-v ${OutputDir}:/input \
-v ${OutputDir}/QC/LongQC/TrimmerLarge/:/output \
longqc sampleqc -x ont-rapid -m 2 -p 5 -o /output/${Samplename}/ /input/CleanedFastQ/TrimmerLarge/${Samplename}_Pore_TrimmerLarge.fastq
done



#A) QC checking of QC passed FastQ data using Fatqc
mkdir -p QC/FastQC/TrimmerLarge_Fastq_QC/
for barcod in ${Barcodes[@]} ; do
fastqc -t 35 -o QC/FastQC/TrimmerLarge_Fastq_QC/ CleanedFastQ/TrimmerLarge/${Samplename}_Pore_TrimmerLarge.fastq ;done 

mkdir -p MultiQC/FastQC/TrimmerLargeFastq_QC
multiqc -o MultiQC/FastQC/TrimmerLargeFastq_QC/ QC/FastQC/TrimmerLarge_Fastq_QC/



##Alignment of reads
mkdir -p BamFiles 

for fastq in ${Barcodes[@]} ; do
minimap2 -t 10 \
    -ax map-ont \
    ${RefFolder}${fasta} CleanedFastQ/TrimmerLarge/"${fastq}_"Pore_TrimmerLarge.fastq | \
    samtools view -bS | samtools sort -o BamFiles/${fastq}_TrimmerLarge_sorted.bam ;done ;wait 


parallel  samtools index ::: BamFiles/*.bam
#Alignmet summary ; pomoxis is used for this; https://github.com/nanoporetech/pomoxis
mkdir -p QC/AlignmentSummary

for bam in ${Barcodes[@]} ; do
stats_from_bam BamFiles/${bam}_TrimmerLarge_sorted.bam \
    > QC/AlignmentSummary/${bam}_TrimmerLarge_sorted.stats & done ;wait 



#C) Mos depth global and target region and gene exon
mkdir -p QC/AlignmentQC/Mosdepth/TargetRegion/

for bam in ${Barcodes[@]} ; do
mosdepth -F 1284 --fast-mode --by ${T_RegionFolder}${T_Region} \
QC/AlignmentQC/Mosdepth/TargetRegion/${bam}_TargetRegion BamFiles/${bam}_TrimmerLarge_sorted.bam \
--thresholds 1,10,20,30,50,100,200,500,600 

mkdir -p QC/AlignmentQC/flagstat
samtools flagstat BamFiles/${bam}_TrimmerLarge_sorted.bam > QC/AlignmentQC/flagstat/${bam}_flagstat.txt 

${Mosedepth_plots}  QC/AlignmentQC/Mosdepth/TargetRegion/*.dist.txt \
--output QC/AlignmentQC/Mosdepth/TargetRegion/${bam}_MosdepthCoverage.html 
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

for bam in ${Barcodes[@]} ; do
longshot \
--force_overwrite \
--sample_id ${bam} \
--bam BamFiles/"${bam}_TrimmerLarge_sorted.bam" \
--ref ${RefFolder}${fasta} \
--out Final_VcfFiles/longshot/Raw/"${bam}_longshot.vcf"; done

#--min_alt_frac 0.3 --min_cov 10

#Adding total depth and Allele fraction to info fields
for bam in ${Barcodes[@]} ; do
${longshot_add_DP_AF} \
-VcfFile Final_VcfFiles/longshot/Raw/"${bam}_longshot.vcf" ;done



for Samplename in ${Barcodes[@]} ; do
gatk VariantFiltration \
	    -V Final_VcfFiles/longshot/Raw/"${Samplename}_longshotcorrected.vcf" \
	    -filter "QUAL < 30.0" --filter-name "QUAL30" \
	    -filter "C_AF < 0.3" --filter-name "AF_LT_0.3" \
	    -filter "C_TDP < 20.0" --filter-name "DP_LT_20" \
	   -O ${Samplename}_longshot_filtered.vcf.gz 

tabix --force -p vcf ${Samplename}_longshot_filtered.vcf.gz

gatk --java-options '-Xmx16g -XX:ParallelGCThreads=4' LeftAlignAndTrimVariants \
-R ${RefFolder}/${fasta} -V ${Samplename}_longshot_filtered.vcf.gz -O ${Samplename}_temp.vcf \
--split-multi-allelics
bgzip ${Samplename}_temp.vcf ; tabix -p vcf ${Samplename}_temp.vcf.gz
#${VT} decompose ${Samplename}_longshot_filtered.vcf.gz > temp.vcf
#${VT} normalize -r ${RefFolder}/${fasta} temp.vcf > ${Samplename}_temp.vcf

#bcftools view ${Samplename}_longshot_filtered_Normalised.vcf.gz --regions-file  ${T_RegionFolder}/${T_Region} > ${Samplename}_temp.vcf 
bcftools sort -Oz  ${Samplename}_temp.vcf.gz -o ${Samplename}_longshot_filtered_Normalised.vcf.gz
tabix -p vcf ${Samplename}_longshot_filtered_Normalised.vcf.gz 
done 


rm -rf Final_VcfFiles/longshot/Raw/*
mv *longshot_filtered_Normalised.vcf.gz* Final_VcfFiles/longshot/Raw/ ; rm *vcf*


###ANNOTATION OF VCF FILE USING aNNOVAR
for Samplename in ${Barcodes[@]} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/longshot/Raw/${Samplename}_longshot_filtered_Normalised.vcf.gz \
		-outfile ${Samplename}.avinput \
		-includeinfo -withzyg ;done


for SAMPLENAME in ${Barcodes[@]} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME}_longshot_filtered_Normalised \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 



##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes[@]} ; do 
	${AnnovarFiltering}  --AnnovarFile ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multianno.txt
done 


#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes[@]} ; do 
    echo ${SAMPLENAME}
	find="Otherinfo4"
	replace=$(zgrep "#CHROM"  Final_VcfFiles/longshot/Raw/${SAMPLENAME}_longshot_filtered_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(zgrep "#CHROM" Final_VcfFiles/longshot/Raw/${SAMPLENAME}_longshot_filtered_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv
done


for SAMPLENAME in ${Barcodes[@]} ; do 
 ${longshot_multianno} \
 -AnnovarOutputFile ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multiannoDesiredColumns.tsv  

 ${longshot_multianno} \
 -AnnovarOutputFile ${SAMPLENAME}_longshot_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv
#rm ${SAMPLENAME}.avinput ${SAMPLENAME}_longshot_QC.hg38_multianno.txt
done


mkdir -p  Annovar_Output/LongShoot/Filtered
mkdir -p  Annovar_Output/LongShoot/NotFiltered

mv *_longshot_filtered_Normalised.hg38_multiannoDesiredColumns.tsv Annovar_Output/LongShoot/NotFiltered
mv *_longshot_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv Annovar_Output/LongShoot/Filtered
rm barcode*



#######-----------------------------------------------------Clair3---------------------------------------------------------------------------


#Clair3

cp ${RefFolder}/${fasta} `pwd`/BamFiles/ 
cp ${RefFolder}/${fasta}.fai `pwd`/BamFiles/ 
cp ${T_RegionFolder}${T_Region} `pwd`/BamFiles/TargetRegion.bed 


for bam in ${Barcodes[@]} ; do
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
  --sample_name=${bam} \
   --gvcf  \
  --output=${OUTPUT_DIR}  ; done        #if want to generate  gvcf use --gvcf    


##This steps is for gvcf generation
#for Samplename in ${Barcodes[@]} ; do
#echo "4.B) GenotypeGVCFs ; create single sample vcf file ; https://gatk.broadinstitute.org/hc/en-us/articles/360036899732-GenotypeGVCFs"
#gatk --java-options "-Xmx16g -XX:ParallelGCThreads=12"  GenotypeGVCFs \
# -R ${RefFolder}/${fasta} \
# -V VcfFiles/Clair3/${Samplename}/merge_output.gvcf.gz \
# -O ${Samplename}_raw.vcf.gz \
# -A StrandBiasBySample \
# -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation & 
# done


for Samplename in ${Barcodes[@]} ; do
${Clair3_IFO_FIELDS_ADDITION}  -VcfFile VcfFiles/Clair3/${Samplename}/merge_output.vcf.gz
done


for Samplename in ${Barcodes[@]} ; do

#If you used GATK used below and comment next
#gatk VariantFiltration \
#	    -V "${Samplename}_raw.vcf.gz" \
#	    -filter "QUAL < 30.0" --filter-name "QUAL30" \
#	    -filter "AF < 0.3" --filter-name "AF_LT_0.3" \
#	    -filter "DP < 20.0" --filter-name "DP_LT_20" \
#	   -O ${Samplename}_Clair3_filtered.vcf.gz 

gatk VariantFiltration \
	    -V VcfFiles/Clair3/${Samplename}/merge_output.vccorrected.vcf \
	    -filter "QUAL < 20.0" --filter-name "QUAL20" \
	    -filter "Reported_AF < 0.3" --filter-name "AF_LT_0.3" \
	    -filter "DP < 20.0" --filter-name "DP_LT_20" \
	   -O ${Samplename}_Clair3_filtered.vcf.gz 

tabix --force -p vcf ${Samplename}_Clair3_filtered.vcf.gz

#gatk --java-options '-Xmx16g -XX:ParallelGCThreads=4' LeftAlignAndTrimVariants \
#-R ${RefFolder}/${fasta} -V ${Samplename}_Clair3_filtered.vcf.gz -O ${Samplename}_Clair3_filtered_Normalised.vcf \
#--split-multi-allelics


${VT} decompose ${Samplename}_Clair3_filtered.vcf.gz > temp.vcf
${VT} normalize -r ${RefFolder}/${fasta} temp.vcf > ${Samplename}_temp.vcf
bgzip --force ${Samplename}_temp.vcf ; tabix --force -p vcf ${Samplename}_temp.vcf.gz

#bcftools view ${Samplename}_Clair3_filtered_Normalised.vcf.gz --regions-file  ${T_RegionFolder}/${T_Region} > ${Samplename}_temp.vcf 
#bcftools sort -Oz  ${Samplename}_temp.vcf -o ${Samplename}_Clair3_filtered_Normalised_Target.vcf.gz

bcftools sort -Oz ${Samplename}_temp.vcf.gz -o ${Samplename}_Clair3_filtered_Normalised.vcf.gz
tabix --force -p vcf ${Samplename}_Clair3_filtered_Normalised_Target.vcf.gz 

done 


mkdir -p Final_VcfFiles/Clair3/Raw/*
mv *Clair3_filtered_Normalised.vcf.gz* Final_VcfFiles/Clair3/Raw/ ; rm *vcf*




###ANNOTATION OF VCF FILE USING aNNOVAR

for Samplename in ${Barcodes[@]} ; do 
	perl ${Annovar}convert2annovar.pl \
		-format vcf4 Final_VcfFiles/Clair3/Raw/${Samplename}_Clair3_filtered_Normalised.vcf.gz \
		-outfile ${Samplename}.avinput \
		-includeinfo -withzyg ;done


for SAMPLENAME in ${Barcodes[@]} ; do 
perl ${Annovar}table_annovar.pl \
${SAMPLENAME}.avinput ${AnnovarDb} \
-buildver hg38 \
-out ${SAMPLENAME}_Clair3_filtered_Normalised \
-remove -otherinfo \
-protocol refGene,\
gnomad211_exome,1000g2015aug_all,esp6500siv2_all,gnomad30_genome,kaviar_20150923,\
avsnp150,dbnsfp42a,dbscsnv11,intervar_20180118,revel,clinvar_20210501,genomicSuperDups,dgvMerged,gwasCatalog \
-operation gx,f,f,f,f,f,f,f,f,f,f,f,r,r,r \
-nastring . -polish -xref ${AnnovarDb}gene_fullxref.txt 
done 



##Filter annovar output file based on the defined parameters
for SAMPLENAME in ${Barcodes[@]} ; do 
	python ${Annovar}AnnovarFiltering.py  --AnnovarFile ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multianno.txt
done 


#Edit column headers and Incorporate additional columns
for SAMPLENAME in ${Barcodes[@]} ; do 
	find="Otherinfo4"
	replace=$(zgrep "#CHROM"  Final_VcfFiles/Clair3/Raw/${SAMPLENAME}_Clair3_filtered_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multiannoDesiredColumns.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multiannoDesiredColumns.tsv 

	find="Otherinfo4"
	replace=$(zgrep "#CHROM" Final_VcfFiles/Clair3/Raw/${SAMPLENAME}_Clair3_filtered_Normalised.vcf.gz)
	sed -i "s+${find}+${replace}+g" ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv
	sed -i 's/\S*Otherinfo\S*//g' ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv

done


for SAMPLENAME in ${Barcodes[@]} ; do 
 python /data/Software_Resources/NGC_scripts/Clair3_RawData_hg38_multianno_DesiredColumns_Filtered.py \
 -AnnovarOutputFile ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multiannoDesiredColumns.tsv  

 python /data/Software_Resources/NGC_scripts/Clair3_RawData_hg38_multianno_DesiredColumns_Filtered.py \
 -AnnovarOutputFile ${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv
#rm ${SAMPLENAME}.avinput ${SAMPLENAME}_longshot_QC.hg38_multianno.txt
done


mkdir -p  Annovar_Output/Clair3/Filtered
mkdir -p  Annovar_Output/Clair3/NotFiltered

mv *_Clair3_filtered_Normalised.hg38_multiannoDesiredColumns.tsv Annovar_Output/Clair3/NotFiltered
mv *_Clair3_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv Annovar_Output/Clair3/Filtered
rm barcode*



for SAMPLENAME in ${Barcodes[@]} ; do 
${Longshot_clair_AnnovarOutputMerg} \
-ClairAnnovarFile Annovar_Output/Clair3/NotFiltered/${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multiannoDesiredColumns.tsv \
-LongshotAnnovarFile Annovar_Output/LongShoot/NotFiltered/${SAMPLENAME}_longshot_filtered_Normalised.hg38_multiannoDesiredColumns.tsv 
done


for SAMPLENAME in ${Barcodes[@]} ; do 
${Longshot_clair_AnnovarOutputMerg} \
-ClairAnnovarFile Annovar_Output/Clair3/Filtered/${SAMPLENAME}_Clair3_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv \
-LongshotAnnovarFile Annovar_Output/LongShoot/Filtered/${SAMPLENAME}_longshot_filtered_Normalised.hg38_multianno_DesiredColumns_Filtered.tsv 
done

mkdir -p Annovar_Output/Clair3_Longshot_Merged/Filtered
mkdir -p Annovar_Output/Clair3_Longshot_Merged/Not_Filtered
mv *hg38_multianno_DesiredColumns_Filtered_clair_long.csv Annovar_Output/Clair3_Longshot_Merged/Filtered
mv *hg38_multiannoDesiredColumns_clair_long.csv Annovar_Output/Clair3_Longshot_Merged/Not_Filtered



python -W ignore ${Amplicon_sequencing_demultiplex} \
-SampleDetails SampleDetails.tsv \
-Library AmpliconLibraryDistribution.tsv \
-Coverage QC/AlignmentQC/Mosdepth/TargetRegion/20X.csv \
-Folder Annovar_Output/Clair3_Longshot_Merged/Not_Filtered/ \
--EndofFile hg38_multiannoDesiredColumns_clair_long.csv

${Amplicon_Variant_SampleSummary}

