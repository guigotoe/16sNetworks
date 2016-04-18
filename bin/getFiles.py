#!usr/bin/python

######################################################
## By Guillermo Torres PhD. - IKMB                 ##
####################################################

## Libraries importing

import sys,os,re,codecs

## Parsing Options ##


## Global variables ##sys.argv[1]#'
inf = '/home/torres/ikmb_storage/Metagenome/rawdata/Metagenome/FF_Aging-Foxo-Microbiomes-LLI_1,/home/torres/ikmb_storage/Metagenome/rawdata/Metagenome/last9/FF_Aging-Foxo-Microbiomes-LLI_1' #sys.argv[1] #paths array, separated by comma '/home/torres/Documents/Projects/Metagenome/resutls/all/16s.file' #sys.argv[1]
inf = inf.split(',')

## Main ##

def main():
    gzfiles = [];fqfiles = {};togzip=[]#;contigs = {}
    for i in inf:
        f1 = os.listdir(i)
        for j in f1:
            try:
                f2 = os.listdir("%s/%s"%(i,j))#;contigs[j] = [j]
                fqfiles[j] = [j,'','']
                for k in f2:
                    if re.search("\.gz$",k):
                        gzfiles.append("%s/%s/%s"%(i,j,k))
                        if re.search("R1_001\.fastq\.gz$",k):fq= re.sub("\.gz$","",k);fqfiles[j][1] = "%s/%s/%s"%(i,j,fq)#;rp = re.sub("R1_001\.fastq\.gz$",'1_P.fq',k)
                        elif re.search("R2_001\.fastq\.gz$" ,k):fq = re.sub("\.gz$","",k);fqfiles[j][2] = "%s/%s/%s"%(i,j,fq)#;rp = re.sub("R2_001\.fastq\.gz$" ,'2_U.fq',k)
                    #if re.search("\.fastq$",k):
                    #    togzip.append("%s/%s/%s"%(i,j,k))
                    #    if re.search("R1_001\.fastq$",k):rp = re.sub("R1_001\.fastq$",'1_P.fq',k);fqfiles[j][1] = "%s/%s/%s"%(i,j,rp)
                    #    elif re.search("R2_001\.fastq$",k):rp = re.sub("R2_001\.fastq$" ,'2_U.fq',k);fqfiles[j][2] = "%s/%s/%s"%(i,j,rp)
                    #    fqfiles[j].append("%s/%s/%s"%(i,j,k))#;contigs[j].append("%s/%s/%s"%(i,j,rp))
            except OSError:
                pass
    #print fqfiles

    gzout = open('gzfile.txt','w');fqout = open('fqfiles.txt','w')#;togzout = open('togzip.txt','w')#;fcout = open('contigs.txt','w')
    gzout.write("\n".join(str(elem) for elem in gzfiles))
    #togzout.write("\n".join(str(elem) for elem in togzip))
    for key,val in fqfiles.items():fqout.write("\t".join(val[0:3])+"\n")

if __name__ == '__main__':main()