
import argparse
from cyvcf2 import VCF, Writer

parser=argparse.ArgumentParser(description="It is for generating Manhatton and qq plot from genesis software")
parser.add_argument('-VcfFile','--VcfFile', help="Vcf fie with path", required=True)
args=parser.parse_args()

VcfFile=args.VcfFile

vcf = VCF(VcfFile)


vcf.add_info_to_header({'ID':'C_AF', 'Description':'Alternate Allelic fraction','Type':'Float', 'Number':'1'})
vcf.add_info_to_header({'ID':'C_TDP', 'Description':'Totall Depth Reference and varant alle together','Type':'Float', 'Number':'1'})

fname=VcfFile
fname=fname[:-4]+"corrected.vcf"
w = Writer(fname, vcf)

for variant in vcf:
        RefDP=variant.INFO.get('AC')[0]
        VDP=variant.INFO.get('AC')[1]
        variant.INFO['C_AF']=str(VDP/(RefDP+VDP))
        variant.INFO['C_TDP']=str(RefDP+VDP)
        w.write_record(variant)

w.close(); vcf.close()

