#!/usr/bin/env python
# coding: utf-8

# In[79]:


import pandas as pd
import glob
import argparse
pd.options.display.float_format = "{:,.2f}".format


# In[ ]:


parser=argparse.ArgumentParser(description="It is for summarising Mosedepth TargetRegion.thresholds outpt")
parser.add_argument('-Mfolder','--MosedepthOutPutFolder', help="Mosedepth Output Folder", required=True)
parser.add_argument('-Outputfolder','--Outputfolder', help="OutputFolderTo save the results; folder should end with / ; then file name can provid", required=True)
args=parser.parse_args()

Folder=args.MosedepthOutPutFolder
Outputfolder=args.Outputfolder

# In[80]:


files = glob.glob(Folder+"*thresholds.bed.gz")


# In[149]:


def greetings(coverage):
    FName=coverage
    FName=pd.read_csv(files[0],sep="\t").iloc[:,0:4]
    FName['TotalRegion']=FName['end']-FName['start']
    
    for file in files:
        df=pd.read_csv(file,sep="\t")[['#chrom','start','end','region',coverage]]
        NewColnames=str(file.rsplit('/', 1)[1]).rsplit('_', 1)[0]+":>="+coverage
        df = df.rename(columns = {coverage:NewColnames})
        FName=pd.merge(FName,df,on=['#chrom','start','end','region'])
        FName[NewColnames]=FName[NewColnames]/FName['TotalRegion']*100
        FName.to_csv(Outputfolder+coverage+".csv",sep="\t",float_format='%.2f',index=False)

    return FName


# In[150]:

#try:
coverage10x=greetings('10X')
#except:
 #   pass
try:
    overage20x=greetings('20X')
except:
    pass
try:
    coverage50x=greetings('50X')
except:
    pass
try:
    coverage100x=greetings('100X')
except:
    pass



try:
        coverage100x=greetings('200X')
except:
        pass

try:
        coverage100x=greetings('500X')
except:
        pass


try:
        coverage100x=greetings('800X')
except:
        pass

try:
            coverage100x=greetings('1000X')
except:
            pass


# In[ ]:




