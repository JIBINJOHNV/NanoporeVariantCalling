

import pandas as pd
import glob
import argparse
pd.options.display.float_format = "{:,.2f}".format


# 

parser=argparse.ArgumentParser(description="It is for summarising Mosedepth TargetRegion.thresholds outpt")
parser.add_argument('-Mfolder','--MosedepthOutPutFolder', help="Mosedepth Output Folder", required=True)
args=parser.parse_args()

Folder=args.MosedepthOutPutFolder

# In[80]:

files = glob.glob(Folder+"*thresholds.bed.gz")


for i in files :
    
    File=pd.read_csv(i,sep="\t")
    
    File["TotalBases"]=File["end"]-File["start"]
    
    DropCol=['#chrom', 'start', 'end', 'region','TotalBases']
    Columns=File.columns
    Columns.drop(DropCol)
    
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].apply(pd.to_numeric)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].div(File.TotalBases, axis=0)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].mul(100)
    File[Columns.drop(DropCol)]=File[Columns.drop(DropCol)].astype(int)
    
    File.to_csv(i[:-7]+"_percentage_bed.tsv",sep="\t",float_format='%.2f',index=False)
