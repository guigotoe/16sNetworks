#!usr/bin/python

######################################################
## By Guillermo Torres PhD. - IKMB                 ##
####################################################

## Libraries importing

import sys,os,re,codecs

## Parsing Options ##


## Global variables ##
inf = sys.argv[1] #'/home/torres/Documents/Projects/Metagenome/resutls/all/16s.file' #sys.argv[1]
outf = inf+'_qc.txt'

## Main ##

def main():
    out = open(outf,'w')

    for line in codecs.open(inf,'r','utf-8'):
        l = line.strip("\n").split("\t")
        r1 = l[1]
        r2 = l[2]
        rp1 = re.sub("1\.fastq","1_P.fq" ,l[1])
        rp2 = re.sub("2\.fastq" ,'2_P.fq' ,l[2])
        ru1 = re.sub("1\.fastq","1_U.fq" ,l[1])
        ru2 = re.sub("2\.fastq" ,'2_U.fq' ,l[2])
        path = '/'.join(r1.split('/')[:-1])

        trimm = "java -jar ~/Bin/trimmomatic-0.33.jar PE -threads 4 %s %s %s %s %s %s ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:5:20 MINLEN:100 2> %s/trimmomatic.log " \
        % (r1,r2,rp1,ru1,rp2,ru2,path)
        print trimm
        os.system(trimm)
        out.write("%s\n"%'\t'.join([l[0],rp1,rp2]))
        # HEADCROP:12
    print 'Done! new file list named as: %s' %outf
## Functions ##


if __name__ == '__main__':main()