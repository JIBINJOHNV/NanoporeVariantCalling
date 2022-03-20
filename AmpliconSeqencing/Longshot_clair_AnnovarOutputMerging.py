

import pandas as pd
import numpy as np
import argparse

parser=argparse.ArgumentParser(description="It is for generating Manhatton and qq plot from genesis software; eg Annovar_Output/Clair3/Filtered/barcode60_Clair3_filtered_Normalised_Target.hg38_multianno_DesiredColumns_Filtered.tsv\n")
parser.add_argument('-ClairAnnovarFile','--ClairAnnovarFile', help="Files generated from the clair vcf annovar annotation ; eg: Annovar_Output/LongShoot/Filtered/barcode60_longshot_filtered_Normalised_Target.hg38_multianno_DesiredColumns_Filtered.tsv", required=True)
parser.add_argument('-LongshotAnnovarFile','--LongshotAnnovarFile', help="Files generated from the Longshot vcf annovar annotation ; eg barcode60_Normalised_Target.hg38_multianno_DesiredColumns_Filtered", required=True)

parser.add_argument('-OutputFilename','--OutputFilename', help="Output file name", required=False)
args=parser.parse_args()


clar=args.ClairAnnovarFile
long=args.LongshotAnnovarFile
outputname="_".join(clar.split("/")[-1].split("_")[2:])[:-4]


#clar="Annovar_Output/Clair3/Filtered/barcode60_Clair3_filtered_Normalised_Target.hg38_multianno_DesiredColumns_Filtered.tsv"
#long="Annovar_Output/LongShoot/Filtered/barcode60_longshot_filtered_Normalised_Target.hg38_multianno_DesiredColumns_Filtered.tsv"


clairdf=pd.read_csv(clar,sep="\t")
longdf=pd.read_csv(long,sep="\t")

CSampleName=" ".join(pd.DataFrame(clairdf.iloc[:,-1]).columns)
LSampleName=" ".join(pd.DataFrame(longdf.iloc[:,-1]).columns)


claircolumns=['RefinedGenotype', 'Zygosity', 'Depth','VariantAlleleFrequency', 'QUAL', 'FILTER','#CHROM','POS','ID', 'REF', 'ALT', 'INFO','FORMAT', CSampleName]
longcolumns=['RefinedGenotype', 'Zygosity', 'Depth','VariantAlleleFrequency', 'QUAL', 'FILTER','#CHROM', 'POS', 'ID', 'REF', 'ALT', 'INFO','FORMAT', LSampleName]




CommonColumnsdf = pd.concat([clairdf.drop(claircolumns+['QUAL.1','FILTER.1'],axis=1),longdf.drop(longcolumns+['QUAL.1','FILTER.1'],axis=1)]).drop_duplicates()
CommonColumnsdf.sort_values(by=['#Chr', 'Start', 'End', 'Ref', 'Alt'], inplace = True)
CommonColumnsdf.drop_duplicates(subset=['#Chr', 'Start', 'End', 'Ref', 'Alt'], keep='first',inplace = True)
CommonColumnsdf[['#Chr', 'Start', 'End', 'Ref', 'Alt']]=CommonColumnsdf[['#Chr', 'Start', 'End', 'Ref', 'Alt']].astype(str)


clairdfSpe=clairdf[claircolumns]
longdfSpe=longdf[longcolumns]


clairdfSpe.columns = ['clair_' +'RefinedGenotype', 'clair_' +'Zygosity', 'clair_' +'Depth','clair_' +'VariantAlleleFrequency', 'clair_' +'QUAL', 'clair_' +'FILTER',
                        '#CHROM','POS','ID', 'REF', 'ALT', 'clair_' +'INFO', 'clair_' +'FORMAT', 'clair_' +CSampleName]


longdfSpe.columns = ['long_' +'RefinedGenotype', 'long_' +'Zygosity', 'long_' +'Depth','long_' +'VariantAlleleFrequency', 'long_' +'QUAL', 'long_' +'FILTER',
                        '#CHROM','POS','ID', 'REF', 'ALT', 'long_' +'INFO', 'long_' +'FORMAT', 'long_' +LSampleName]

clairdfSpe[clairdfSpe.columns]=clairdfSpe[clairdfSpe.columns].astype(str)
longdfSpe[longdfSpe.columns]=longdfSpe[longdfSpe.columns].astype(str)


Spec_Clairlong=pd.merge(longdfSpe,clairdfSpe,on=['#CHROM','POS','ID', 'REF', 'ALT'],how="outer")
Spec_Clairlong.fillna("NA", inplace = True) 



Spec_Clairlong['clair|long_RefinedGenotype'] = Spec_Clairlong[['clair_RefinedGenotype','long_RefinedGenotype']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_Zygosity'] = Spec_Clairlong[['clair_Zygosity','long_Zygosity']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_Depth'] = Spec_Clairlong[['clair_Depth','long_Depth']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair_long_VariantAlleleFrequency'] = Spec_Clairlong[['clair_VariantAlleleFrequency','long_VariantAlleleFrequency']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_QUAL'] = Spec_Clairlong[['clair_QUAL','long_QUAL']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_FILTER'] = Spec_Clairlong[['clair_FILTER','long_FILTER']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_INFO'] = Spec_Clairlong[['clair_INFO','long_INFO']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_FORMAT'] = Spec_Clairlong[['clair_FORMAT','long_FORMAT']].apply(lambda x: '|'.join(x), axis = 1)
Spec_Clairlong['clair|long_'+CSampleName] = Spec_Clairlong[['clair_'+CSampleName,'long_'+LSampleName]].apply(lambda x: '|'.join(x), axis = 1)


Spec_Clairlong=Spec_Clairlong[['#CHROM','POS', 'ID', 'REF', 'ALT','clair|long_'+CSampleName,'clair|long_RefinedGenotype', 'clair|long_Zygosity', 'clair|long_Depth', 'clair|long_QUAL', 'clair|long_FILTER', 
                                'clair|long_INFO', 'clair|long_FORMAT','clair_long_VariantAlleleFrequency']]


Merged=pd.merge(CommonColumnsdf,Spec_Clairlong,left_on=['#Chr', 'Start', 'Ref', 'Alt'],right_on=['#CHROM','POS', 'REF', 'ALT'])


Merged=Merged[['#Chr', 'Start', 'End', 'Ref', 'Alt', 'Gene.refGene', 'Func.refGene',
       'ExonicFunc.refGene', 'GeneDetail.refGene', 'AAChange.refGene',
       'MAX_MAF', 'MAX_MAF_ReportedDatabase', 'GIV_Indian_HOM-VAR',
       'GIV_CDFD_HOM-VAR','clair|long_RefinedGenotype', 'clair|long_Zygosity',
       'clair|long_Depth', 'clair_long_VariantAlleleFrequency','clair|long_QUAL', 'clair|long_FILTER',
       'CADD_phred','MIM_disease.refGene', 'MIM_id.refGene', 'InterVar_automated',
       'Orphanet_disorder.refGene', 'Orphanet_association_type.refGene',
       'Trait_association(GWAS).refGene', 'HPO_id.refGene', 'HPO_name.refGene',
       'CLNDN', 'CLNREVSTAT', 'CLNSIG', 'Tissue_specificity(Uniprot).refGene',
       'SIFT_pred', 'SIFT4G_pred', 'Polyphen2_HDIV_pred',
       'Polyphen2_HVAR_pred', 'LRT_pred', 'MutationTaster_pred',
       'MutationAssessor_pred', 'FATHMM_pred', 'PROVEAN_pred', 'MetaSVM_pred',
       'MetaLR_pred', 'MetaRNN_pred', 'M-CAP_pred', 'PrimateAI_pred',
       'DEOGEN2_pred', 'BayesDel_addAF_pred', 'ClinPred_pred', 'LIST-S2_pred',
       'Aloft_pred', 'fathmm-MKL_coding_pred', 'fathmm-XF_coding_pred',
       'Interpro_domain', 'dbscSNV_ADA_SCORE', 'dbscSNV_RF_SCORE',
       'Gene damage prediction (all Mendelian disease-causing genes).refGene',
       'Gene damage prediction (Mendelian AD disease-causing genes).refGene',
       'Gene damage prediction (Mendelian AR disease-causing genes).refGene',
       'Expression(egenetics).refGene', 'Expression(GNF/Atlas).refGene',
       'LoFtool_score.refGene', 'Essential_gene.refGene',
       'Essential_gene_CRISPR.refGene', 'Essential_gene_CRISPR2.refGene',
       'Essential_gene_gene-trap.refGene',
       'Gene_indispensability_score.refGene',
       'Gene_indispensability_pred.refGene', 'MGI_mouse_gene.refGene',
       'MGI_mouse_phenotype.refGene', 'ZFIN_zebrafish_gene.refGene',
       'ZFIN_zebrafish_structure.refGene',
       'ZFIN_zebrafish_phenotype_quality.refGene',
       'ZFIN_zebrafish_phenotype_tag.refGene', 'GnomAdExome_AF',
       'GnomAdExome_AF_popmax', 'GnomAdExome_controls_AF_popmax',
       'GnomAdGenome_AF', 'GnomAdGenomeAF_popmax', '1000g2015aug_all',
       'esp6500siv2_all', 'Kaviar_AF', 'abraom_freq', 'abraom_filter',
       'GME_AF_popmax', 'GIV_Indian_AAF', 'GIV_Indian_HET', 'GIV_CDFD_AAF',
       'GIV_CDFD_HET', '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'clair|long_INFO', 'clair|long_FORMAT', 'clair|long_'+CSampleName]]

Merged.to_csv(outputname+'clair_long.csv',index=None)

