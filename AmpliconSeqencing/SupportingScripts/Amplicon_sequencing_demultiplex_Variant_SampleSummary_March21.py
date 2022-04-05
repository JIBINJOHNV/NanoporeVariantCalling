

import pandas as pd
import glob
import argparse
import os
import numpy as np
pd.options.display.float_format = "{:,.2f}".format






PrioritisedDf=pd.DataFrame()
files = glob.glob("SampleWise_Annovar/Present/*")
for file in files:
    df=pd.read_csv(file)
    df.rename(columns={" ".join([x for x in df.columns if "barcode" in x]):"VCF_Genotype"},inplace=True)
    PrioritisedDf=PrioritisedDf.append(df)

N_sampleWithVariant=str(len(set(PrioritisedDf['Sample ID'])))
N_sample=str(PrioritisedDf.shape[0])

PrioritisedDf=PrioritisedDf.drop(['Chromosome','Start.1','End.1'],axis=1)
TotalSample=str(len(PrioritisedDf['Sample ID'].unique())) #Findng total number of samples without prioitised mutation
PrioritisedDf.rename(columns={"f_NAME":"FragmentsinWHichMutationObserved"},
                     inplace=True)
PrioritisedDf.drop(["Library"],axis=1,inplace=True)


PrioritisedDf=PrioritisedDf[list(PrioritisedDf.iloc[:,:15].columns)+['Sequenced_Library','Sample Name','Sample ID','Disease','Gene','FragmentsinWHichMutationObserved','No_of_Fragments','Sequenced_Fragemnts_Name_count','Sequenced_Fragemnts_Name',
    'Amplicon_with_lt<95%_20Xcoverage', 'Amplicon_with_Mt>=95%_20Xcoverage']+['#chrom', 'start', 'end','TotalRegion']+list(PrioritisedDf.iloc[:,15:-15].columns)]
PrioritisedDf.replace({'{}':np.nan},inplace=True)
PrioritisedDf=PrioritisedDf.reset_index().drop('index',axis=1)
PrioritisedDf.to_csv("Total_"+TotalSample+"_Samples_With_PrioritisedVariannts_Annovar.csv")






#________________________________
NoPrioritisedVariants_df=pd.DataFrame()

files2 = glob.glob("SampleWise_Annovar/NotPresent/*")

for file in files2:
    df=pd.read_csv(file)
    #df.rename(columns={" ".join([x for x in df.columns if "barcode" in x]):"VCF_Genotype"},inplace=True)
    NoPrioritisedVariants_df=NoPrioritisedVariants_df.append(df)

N_sampleWithVariant=str(len(set(NoPrioritisedVariants_df['Sample ID'])))
N_sample=str(NoPrioritisedVariants_df.shape[0])

TotalSample=str(len(NoPrioritisedVariants_df['Sample ID'].unique())) #Findng total number of samples without prioitised mutation

#Creating New column to say whether require attention or not
NoPrioritisedVariants_df.loc[NoPrioritisedVariants_df['NofFragments']!=NoPrioritisedVariants_df['Sequenced_Fragemnts_Name_count'], 
                             'AttentionRequire_OrNot'] = "Yes"
NoPrioritisedVariants_df.loc[NoPrioritisedVariants_df['Amplicon_with_lt<95%_20Xcoverage']!='{}', 
                             'AttentionRequire_OrNot'] = "Yes"
NoPrioritisedVariants_df['AttentionRequire_OrNot'].replace({np.nan:'No'},inplace=True)

#Select Desired COlumns
NoPrioritisedVariants_df=NoPrioritisedVariants_df[['Sample ID', 'Sample Name', 'Disease','Gene','AttentionRequire_OrNot',
                                                   'Sequenced_Library','NofFragments','Sequenced_Fragemnts_Name_count',
                                                   'Sequenced_Fragemnts_Name','Amplicon_with_lt<95%_20Xcoverage',
                                                   'Amplicon_with_Mt>=95%_20Xcoverage']]

NoPrioritisedVariants_df['AttentionRequire_OrNot'].replace({np.nan:'No'},inplace=True)
NoPrioritisedVariants_df=NoPrioritisedVariants_df.reset_index().drop('index',axis=1)
NoPrioritisedVariants_df.replace({'{}':np.nan},inplace=True)
attention=str(len(NoPrioritisedVariants_df[NoPrioritisedVariants_df['AttentionRequire_OrNot']=='Yes']))

NoPrioritisedVariants_df.to_csv("Total_"+TotalSample+"_Samples_With_No_PrioritisedVariannts_"+attention+"_Samples_RequireaAttention.csv")
