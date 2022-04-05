#!/home/cdfdhyd/miniconda3/envs/Exome/bin/python

import pandas as pd
import glob
import argparse



parser=argparse.ArgumentParser(description="It is for Incorporating additional information and Initial filtering of Annovar out put file")
parser.add_argument('-AnnovarFile','--AnnovarFile', help="Name of Annovar output File", required=True)
args=parser.parse_args()

AnnoFile=args.AnnovarFile



IndianFile="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/hg38_Generic_All-Indian_Variants.txt"
InhouseFile="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/hg38_Generic_In-house_Variants.txt"
gmeFile="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/hg38_gme.txt"
abraomFile="/data/Software_Resources/ReferenceFiles/Homo_sapiens/annovar/humandb_hg38/hg38_abraom.txt"



Indiandb=pd.read_csv(IndianFile,sep="\t")
Indiandb['#Chr'] = 'chr' + Indiandb['#Chr'].astype(str)
Indiandb[['#Chr', 'Start', 'End']] = Indiandb[['#Chr', 'Start', 'End']].astype('str')


Inhousedb=pd.read_csv(InhouseFile,sep="\t")
Inhousedb[['#Chr', 'Start', 'End']] = Inhousedb[['#Chr', 'Start', 'End']].astype('str')

Annodb=pd.read_csv(AnnoFile,sep="\t")
Annodb.rename(columns={"Chr": "#Chr"}, inplace=True)
Annodb[['#Chr', 'Start', 'End']] = Annodb[['#Chr', 'Start', 'End']].astype('str')


gmedb=pd.read_csv(gmeFile,sep="\t")
gmedb['#Chr'] = 'chr' + gmedb['#Chr'].astype(str)
gmedb[['#Chr', 'Start', 'End']] = gmedb[['#Chr', 'Start', 'End']].astype('str')

abraomdb=pd.read_csv(abraomFile,sep="\t")
abraomdb['#Chr'] = 'chr' + abraomdb['#Chr'].astype(str)
abraomdb[['#Chr', 'Start', 'End']] = abraomdb[['#Chr', 'Start', 'End']].astype('str')


#Merge Files
Annodb_India=pd.merge(Annodb,Indiandb,on=["#Chr","Start","End", "Ref","Alt"],how='left')
Annodb_India_Inhouse_1=pd.merge(Annodb_India,Inhousedb,on=["#Chr","Start","End", "Ref","Alt"],how='left')

Annodb_India_Inhouse_1=pd.merge(Annodb_India_Inhouse_1,gmedb,on=["#Chr","Start","End", "Ref","Alt"],how='left')
Annodb_India_Inhouse=pd.merge(Annodb_India_Inhouse_1,abraomdb,on=["#Chr","Start","End", "Ref","Alt"],how='left')





######################-------------------MAF---------Based Analysis-----------------------------------------------------------------------------------
#GNOMad Exome
GnomAdExome=['AF', 'AF_popmax', 'AF_male', 'AF_female', 'AF_raw', 'AF_afr', 'AF_sas', 'AF_amr', 'AF_eas', 'AF_nfe', 'AF_fin', 
            'AF_asj', 'AF_oth', 'non_topmed_AF_popmax', 'non_neuro_AF_popmax','non_cancer_AF_popmax', 'controls_AF_popmax']

Annodb_India_Inhouse[GnomAdExome]=Annodb_India_Inhouse[GnomAdExome].apply(pd.to_numeric, errors='coerce').fillna(0)
for x in GnomAdExome:
  Annodb_India_Inhouse.rename(columns={x:"GnomAdExome_"+x}, inplace=True)


#GNOMAD_GENOME
GnomAdGenome=['AF.1', 'AF_raw.1', 'AF_male.1', 'AF_female.1', 'AF_afr.1', 'AF_amr.1', 'AF_asj.1', 'AF_eas.1', 'AF_fin.1', 'AF_nfe.1', 
            'AF_oth.1', 'AF_sas.1']

Annodb_India_Inhouse[GnomAdGenome+['AF_ami']]=Annodb_India_Inhouse[GnomAdGenome+['AF_ami']].apply(pd.to_numeric, errors='coerce').fillna(0)

Annodb_India_Inhouse.rename(columns={'AF_ami':"GnomAdGenome_AF_ami"},inplace=True)
for x in GnomAdGenome:
  Annodb_India_Inhouse.rename(columns={x:"GnomAdGenome_"+x[:-2]}, inplace=True)


Annodb_India_Inhouse['GnomAdGenomeAF_popmax']=Annodb_India_Inhouse[['GnomAdGenome_AF_afr', 'GnomAdGenome_AF_ami', 'GnomAdGenome_AF_amr', 
                                                                    'GnomAdGenome_AF_asj', 'GnomAdGenome_AF_eas', 'GnomAdGenome_AF_fin', 
                                                                    'GnomAdGenome_AF_nfe', 'GnomAdGenome_AF_oth', 'GnomAdGenome_AF_sas']].max(axis=1)


GME=['GME_NWA','GME_NEA','GME_AP','GME_Israel','GME_SD','GME_TP','GME_CA']
Annodb_India_Inhouse[GME+['GME_AF']]=Annodb_India_Inhouse[GME+['GME_AF']].apply(pd.to_numeric, errors='coerce').fillna(0)
Annodb_India_Inhouse['GME_AF_popmax']=Annodb_India_Inhouse[GME].max(axis=1)

## Identifying maximum reported MAF frquency
MaxMAF_Columns=['GnomAdExome_AF', 'GnomAdExome_AF_afr','GnomAdExome_AF_sas', 'GnomAdExome_AF_amr', 
'GnomAdExome_AF_eas', 'GnomAdExome_AF_nfe', 'GnomAdExome_AF_fin', 'GnomAdExome_AF_asj', 'GnomAdExome_AF_oth',
'1000g2015aug_all', 'esp6500siv2_all', 
'GnomAdGenome_AF','GnomAdGenome_AF_afr', 'GnomAdGenome_AF_ami', 'GnomAdGenome_AF_amr', 'GnomAdGenome_AF_asj','GnomAdGenome_AF_eas', 
'GnomAdGenome_AF_fin', 'GnomAdGenome_AF_nfe', 'GnomAdGenome_AF_oth','GnomAdGenome_AF_sas', 
'Kaviar_AF', 'abraom_freq', 
'GME_NWA', 'GME_NEA', 'GME_AP', 'GME_Israel', 'GME_SD', 'GME_TP', 'GME_CA']

Annodb_India_Inhouse[MaxMAF_Columns]=Annodb_India_Inhouse[MaxMAF_Columns].apply(pd.to_numeric, errors='coerce').fillna(0)

Annodb_India_Inhouse['MAX_MAF']=Annodb_India_Inhouse[MaxMAF_Columns].max(axis=1)
Annodb_India_Inhouse['MAX_MAF_ReportedDatabase']=Annodb_India_Inhouse[MaxMAF_Columns].idxmax(axis="columns")

Annodb_India_Inhouse['Zygosity']=Annodb_India_Inhouse['Otherinfo1']
Annodb_India_Inhouse['QUAL']=Annodb_India_Inhouse['Otherinfo2']
Annodb_India_Inhouse['Depth']=Annodb_India_Inhouse['Otherinfo3']
Annodb_India_Inhouse['FILTER']=Annodb_India_Inhouse['Otherinfo10']

Annodb_India_Inhouse.drop('Otherinfo1',axis=1,inplace=True)
Annodb_India_Inhouse.drop('Otherinfo2',axis=1,inplace=True)
Annodb_India_Inhouse.drop('Otherinfo3',axis=1,inplace=True)

###################------------------------------------------Select Desired Columns---------------------------------------------------------------------------
MAF_COLUMNS=['GnomAdExome_AF','GnomAdExome_AF_popmax','GnomAdExome_controls_AF_popmax','GnomAdGenome_AF','GnomAdGenomeAF_popmax','1000g2015aug_all', 
'esp6500siv2_all','Kaviar_AF','abraom_freq', 'abraom_filter','GME_AF_popmax',
'GIV_Indian_AAF', 'GIV_Indian_HET', 'GIV_Indian_HOM-VAR', 'GIV_CDFD_AAF', 'GIV_CDFD_HET', 'GIV_CDFD_HOM-VAR'] 



##Select annovar file with desired columns
OtherInfoColumns=[x for x in Annodb_India_Inhouse.columns if "Otherinfo" in x]

PrimaryColumns=['#Chr', 'Start','End','Ref','Alt','Gene.refGene', 
                'Func.refGene','ExonicFunc.refGene','GeneDetail.refGene','AAChange.refGene', 
                'MAX_MAF','MAX_MAF_ReportedDatabase','GIV_Indian_HOM-VAR','GIV_CDFD_HOM-VAR','Zygosity','Depth','QUAL','FILTER',
                'CADD_phred','MIM_disease.refGene', 'MIM_id.refGene', 'InterVar_automated',
                'Orphanet_disorder.refGene', 'Orphanet_association_type.refGene', 'Trait_association(GWAS).refGene', 
                'HPO_id.refGene', 'HPO_name.refGene', 'CLNDN', 'CLNREVSTAT', 'CLNSIG', 'Tissue_specificity(Uniprot).refGene', 
                'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 
                'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'MetaSVM_pred',
                 'MetaLR_pred', 'MetaRNN_pred', 'M-CAP_pred', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred', 
                 'ClinPred_pred', 'LIST-S2_pred', 'Aloft_pred', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 
                 'Interpro_domain', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE',  
                 'Gene damage prediction (all Mendelian disease-causing genes).refGene', 
                 'Gene damage prediction (Mendelian AD disease-causing genes).refGene', 
                 'Gene damage prediction (Mendelian AR disease-causing genes).refGene', 'Expression(egenetics).refGene', 
                 'Expression(GNF/Atlas).refGene', 'LoFtool_score.refGene', 'Essential_gene.refGene', 'Essential_gene_CRISPR.refGene', 
                 'Essential_gene_CRISPR2.refGene', 'Essential_gene_gene-trap.refGene', 'Gene_indispensability_score.refGene', 
                 'Gene_indispensability_pred.refGene', 'MGI_mouse_gene.refGene', 'MGI_mouse_phenotype.refGene', 'ZFIN_zebrafish_gene.refGene', 
                 'ZFIN_zebrafish_structure.refGene', 'ZFIN_zebrafish_phenotype_quality.refGene', 'ZFIN_zebrafish_phenotype_tag.refGene']


Annodb_India_Inhouse_withDesColumns=Annodb_India_Inhouse[PrimaryColumns+MAF_COLUMNS+OtherInfoColumns]

Annodb_India_Inhouse_withDesColumns = Annodb_India_Inhouse_withDesColumns.loc[:,~Annodb_India_Inhouse_withDesColumns.columns.duplicated()]
Annodb_India_Inhouse_withDesColumns[['GIV_Indian_HOM-VAR','GIV_CDFD_HOM-VAR']]=Annodb_India_Inhouse_withDesColumns[['GIV_Indian_HOM-VAR',
    'GIV_CDFD_HOM-VAR']].fillna("NA")

Annodb_India_Inhouse_withDesColumns.to_csv(AnnoFile[:-4]+"DesiredColumns.tsv",index=False,sep="\t")

##Filter variants


FiltDf=Annodb_India_Inhouse_withDesColumns

#print(FiltDf["Func.refGene"].unique())

FiltDf=FiltDf[FiltDf["Func.refGene"].isin(['exonic','splicing','exonic;splicing'])]
#print(FiltDf["Func.refGene"].unique())

NonCodVariation=['synonymous SNV','nonframeshift insertion', 'nonframeshift deletion','unknown']
FiltDf=FiltDf[~FiltDf['ExonicFunc.refGene'].isin(NonCodVariation)]


FiltDf=FiltDf[FiltDf['MAX_MAF']<0.01]

FiltDf.to_csv(AnnoFile[:-4]+"_DesiredColumns_Filtered.tsv",index=False,sep="\t")


