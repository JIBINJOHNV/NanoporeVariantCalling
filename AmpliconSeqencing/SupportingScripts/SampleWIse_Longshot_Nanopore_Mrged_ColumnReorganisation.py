
import pandas as pd
import glob
import argparse
import os


EndofFile="multiannoDesiredColumns.csv"
files = glob.glob("SampleWise_Annovar/*multiannoDesiredColumns.csv")

for file in files:
    df=pd.read_csv(file)
    sampleName=df.columns[-1]
    
    DesiredColumns=['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'Gene.refGene', 'Func.refGene', 'ExonicFunc.refGene', 'GeneDetail.refGene', 
'AAChange.refGene', 'MAX_MAF', 'MAX_MAF_ReportedDatabase', 'GIV_Indian_HOM-VAR', 'GIV_CDFD_HOM-VAR', 'clair|long_RefinedGenotype', 
'clair|long_Zygosity', 'clair|long_Depth', 'clair_long_VariantAlleleFrequency', 'clair|long_QUAL', 'clair|long_FILTER', 'CADD_phred', 
'MIM_disease.refGene', 'MIM_id.refGene', 'InterVar_automated', 'Orphanet_disorder.refGene', 'Orphanet_association_type.refGene', 
'Trait_association(GWAS).refGene', 'HPO_id.refGene', 'HPO_name.refGene', 'CLNDN', 'CLNREVSTAT', 'CLNSIG', 'Tissue_specificity(Uniprot).refGene',
 'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred', 'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred', 'MutationAssessor_pred', 
 'FATHMM_pred', 'PROVEAN_pred', 'MetaSVM_pred', 'MetaLR_pred', 'MetaRNN_pred', 'M-CAP_pred', 'PrimateAI_pred', 'DEOGEN2_pred', 
 'BayesDel_addAF_pred', 'ClinPred_pred', 'LIST-S2_pred', 'Aloft_pred', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred', 'Interpro_domain', 
 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE', 'Gene damage prediction (all Mendelian disease-causing genes).refGene', 
 'Gene damage prediction (Mendelian AD disease-causing genes).refGene', 'Gene damage prediction (Mendelian AR disease-causing genes).refGene', 
 'Expression(egenetics).refGene', 'Expression(GNF/Atlas).refGene', 'LoFtool_score.refGene', 'Essential_gene.refGene', 
 'Essential_gene_CRISPR.refGene', 'Essential_gene_CRISPR2.refGene', 'Essential_gene_gene-trap.refGene', 'Gene_indispensability_score.refGene', 
 'Gene_indispensability_pred.refGene', 'MGI_mouse_gene.refGene', 'MGI_mouse_phenotype.refGene', 'ZFIN_zebrafish_gene.refGene', 
 'ZFIN_zebrafish_structure.refGene', 'ZFIN_zebrafish_phenotype_quality.refGene', 'ZFIN_zebrafish_phenotype_tag.refGene', 'GnomAdExome_AF', 
 'GnomAdExome_AF_popmax', 'GnomAdExome_controls_AF_popmax', 'GnomAdGenome_AF', 'GnomAdGenomeAF_popmax', '1000g2015aug_all', 'esp6500siv2_all', 
 'Kaviar_AF', 'abraom_freq', 'abraom_filter', 'GME_AF_popmax', 'GIV_Indian_AAF', 'GIV_Indian_HET', 'GIV_CDFD_AAF', 'GIV_CDFD_HET', '#CHROM', 
 'POS', 'ID', 'REF', 'ALT', 'clair|long_INFO', 'clair|long_FORMAT', 'clair|long_'+sampleName[:-6], 'Chromosome', 'Start', 'End', 'f_NAME', 'Library', 
 'Sample ID', 'Sample Name', 'Gene', 'No_of_Fragments', 'Disease', '#chrom', 'start', 'end', 'TotalRegion']+[sampleName]
    df=df[DesiredColumns]
    df.to_csv(file,index=None)