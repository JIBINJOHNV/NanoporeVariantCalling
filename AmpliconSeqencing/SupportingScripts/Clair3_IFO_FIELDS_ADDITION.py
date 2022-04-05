


import argparse
from cyvcf2 import VCF, Writer

parser=argparse.ArgumentParser(description="It is for generating Manhatton and qq plot from genesis software")
parser.add_argument('-VcfFile','--VcfFile', help="Vcf fie with path", required=True)
args=parser.parse_args()

VcfFile=args.VcfFile

#VcfFile="VcfFiles/Clair3/barcode60/merge_output.vcf.gz"
#VcfFile="merge_output.vcf.gz"
vcf = VCF(VcfFile)


vcf.add_info_to_header({'ID':'Reported_AF', 'Description':'Alternate Allelic fraction','Type':'Float', 'Number':'1'})
vcf.add_info_to_header({'ID':'DP', 'Description':'Totall Depth Reference and varant alle together','Type':'Float', 'Number':'1'})
vcf.add_info_to_header({'ID':'Wt_Ad', 'Description':'Totall Depth of Reference allele','Type':'Float', 'Number':'1'})
vcf.add_info_to_header({'ID':'V_Ad', 'Description':'Totall Depth of varant alle ','Type':'Float', 'Number':'1'})
vcf.add_info_to_header({'ID':'CalculatedAf', 'Description':'calculated by V_Ad//(Wt_Ad+V_Ad) ','Type':'Float', 'Number':'1'})

fname=VcfFile
fname=fname[:-4]+"corrected.vcf"
w = Writer(fname, vcf)

for variant in vcf:
    DP=int(str(variant.format('DP')[:,0]).split('[')[1].split("]")[0])
    Reported_AF=float(str(variant.format('AF')[:,0]).split('[')[1].split("]")[0])
    Wt_Ad=int(variant.format('AD')[0][0])
    V_Ad=int(variant.format('AD')[0][1])
    CalculatedAf=V_Ad/(Wt_Ad+V_Ad)
    variant.INFO['DP']=DP
    variant.INFO['Reported_AF']=Reported_AF
    variant.INFO['Wt_Ad']=Wt_Ad
    variant.INFO['V_Ad']=V_Ad
    variant.INFO['CalculatedAf']=CalculatedAf
    w.write_record(variant)



w.close(); vcf.close()

