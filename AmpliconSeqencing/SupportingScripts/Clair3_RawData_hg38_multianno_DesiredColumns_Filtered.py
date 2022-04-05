

import os
import subprocess
import pandas as pd
import argparse




parser=argparse.ArgumentParser(description="It is for Correcting the genotype")
parser.add_argument('-AnnovarOutputFile','--AnnovarOutputFile', help="Name of Annovar output file", required=True)

args=parser.parse_args()
AnnovarFile=args.AnnovarOutputFile


#Read the annovar files
AnnovarFile=pd.read_csv(AnnovarFile,sep="\t")

try:
    AnnovarFile=AnnovarFile.drop([x for x in AnnovarFile.columns if "Unnamed" in x],axis=1)
except:
    pass

try:
    AnnovarFile.rename(columns={'O#CHROM':'#CHROM'},inplace=True)
except:
    pass

#Identify the sample name
SampleName=" ".join(pd.DataFrame(AnnovarFile.iloc[:,-1]).columns)

Depth=pd.DataFrame(AnnovarFile[SampleName].str.split(":", expand=True)[[2]])
Depth.rename(columns={2:"TotalDepth"},inplace=True)

#Depth[["ReferenceAlleleDepth","AltAlleleDepth"]]=AnnovarFile[SampleName].str.split(":", 
#                                                                        expand=True)[[1]][1].str.split(",",expand=True)
                                                                        
cols=Depth.columns

#Depth[cols] = Depth[cols].apply(pd.to_numeric, errors='coerce', axis=1)
Depth['VariantAlleleFrequency']=AnnovarFile[SampleName].str.split(":",expand=True)[[4]] 
Depth["VariantAlleleFrequency"] = pd.to_numeric(Depth["VariantAlleleFrequency"])


#Defining the genotype
Depth.loc[Depth['VariantAlleleFrequency']>0.7, 
          'RefinedGenotype'] = "HomozygousVariant"

Depth.loc[(Depth['VariantAlleleFrequency']<=0.7) & (Depth['VariantAlleleFrequency']>=0.30),
          'RefinedGenotype'] = "HeterozygousVariant"

Depth.loc[Depth['VariantAlleleFrequency']<0.3, 
          'RefinedGenotype'] = "WT/SomaticVariant"

AnnovarFile[['VariantAlleleFrequency','RefinedGenotype']]=Depth[['VariantAlleleFrequency', 'RefinedGenotype']]

ColumnOrder=['#Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'GeneDetail.refGene',
 'AAChange.refGene', 'MAX_MAF', 'MAX_MAF_ReportedDatabase', 'GIV_Indian_HOM-VAR', 'GIV_CDFD_HOM-VAR', 
 'RefinedGenotype','Zygosity','Depth', 'VariantAlleleFrequency', 
 'QUAL', 'FILTER', 'CADD_phred', 'MIM_disease.refGene', 'MIM_id.refGene', 'InterVar_automated', 'Orphanet_disorder.refGene',
 'Orphanet_association_type.refGene', 'Trait_association(GWAS).refGene', 'HPO_id.refGene', 'HPO_name.refGene', 'CLNDN',
 'CLNREVSTAT', 'CLNSIG', 'Tissue_specificity(Uniprot).refGene', 'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred',
 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred',
 'MetaSVM_pred', 'MetaLR_pred','MetaRNN_pred', 'M-CAP_pred', 'PrimateAI_pred', 'DEOGEN2_pred', 'BayesDel_addAF_pred',
 'ClinPred_pred', 'LIST-S2_pred', 'Aloft_pred', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'Interpro_domain',
 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'Gene damage prediction (all Mendelian disease-causing genes).refGene',
 'Gene damage prediction (Mendelian AD disease-causing genes).refGene', 'Gene damage prediction (Mendelian AR disease-causing genes).refGene',
 'Expression(egenetics).refGene', 'Expression(GNF/Atlas).refGene', 'LoFtool_score.refGene', 'Essential_gene.refGene',
 'Essential_gene_CRISPR.refGene', 'Essential_gene_CRISPR2.refGene', 'Essential_gene_gene-trap.refGene',
 'Gene_indispensability_score.refGene', 'Gene_indispensability_pred.refGene', 'MGI_mouse_gene.refGene', 
 'MGI_mouse_phenotype.refGene', 'ZFIN_zebrafish_gene.refGene', 'ZFIN_zebrafish_structure.refGene',
 'ZFIN_zebrafish_phenotype_quality.refGene', 'ZFIN_zebrafish_phenotype_tag.refGene', 'GnomAdExome_AF', 'GnomAdExome_AF_popmax',
 'GnomAdExome_controls_AF_popmax', 'GnomAdGenome_AF', 'GnomAdGenomeAF_popmax', '1000g2015aug_all', 'esp6500siv2_all',
 'Kaviar_AF', 'abraom_freq', 'abraom_filter', 'GME_AF_popmax', 'GIV_Indian_AAF', 'GIV_Indian_HET', 'GIV_CDFD_AAF',
 'GIV_CDFD_HET', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL.1', 'FILTER.1', 'INFO', 'FORMAT', SampleName]

AnnovarFile=AnnovarFile[ColumnOrder]
AnnovarFile=AnnovarFile.replace(to_replace =".",value ="NA")


AnnovarFile.to_csv(args.AnnovarOutputFile,sep="\t",index=None)

