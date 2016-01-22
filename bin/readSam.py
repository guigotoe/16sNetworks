#!/usr/bin/python
####################################################
# By Guillermo Torres PhD.c - IKMB                 #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #
####################################################
# Last update: September 2015

import re, sys, os, math,pysam
from optparse import OptionParser

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to read SAM/BAM files and out some statistics.
\t%prog -i sam/bam_file -o outname [options] arg \n\t%prog {-h --help}\n
\t It requires pysam libray
As part of Toolbox scripts by Guillermo G. Torres (ggtorrese@unal.edu.co) -
IKMB - Christian-Albrechts-Universitat zu Kiel (CAU) """
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--infile",type="string",action="store",dest="inf",default='samread/test.sam',
                      help="Input folder path - absolute path is recomended")
parser.add_option("-o","--out_file",type="string",action="store",dest="out",default='out',
                      help="Out files name - default name = cat.fna")
parser.add_option("-D","--diamond", action="store_true",dest="diamond",default=True,
                      help = "Diamond SAM file - default False")
parser.add_option("-v", action="store_true",dest="verbose",default=True,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.inf == None):
    parser = OptionParser(usage="\n\t%prog -i sam/bam_file -o outname [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
path = os.path.abspath(o.inf)
#out = open(o.out,'w')

def main():
    if o.diamond == True:
        sed = 'sed -i \'1s/^/@SQ	SN:	LN:\\n/\' '+path
        unsed = 'sed -i \'1d\' '+path
        os.system(sed)
        sam = pysam.AlignmentFile(o.inf,'r',check_sq=False)
        os.system(unsed)
        for read in sam:
            print read
            print read.query_alignment_length
            print read.query_length

    #fails = []
    #sys.stderr.write("Concatenating...\n")
    #for i in os.listdir(path):
    #    keys,files = {},[]
    #    try:
    #        for j in os.listdir(path+'/'+i):
    #            keys[j] = re.findall("\d+\.fna",j)[0].strip('.fna')
    #            files.append(re.findall("\d+\.fna",j)[0].strip('.fna'))
    #        for k,v in keys.iteritems():
    #            if v == max(files): cat_files.append(path+'/'+i+'/'+k)
    #    except IndexError:
    #        if len(os.listdir(path+'/'+i)) == 1:
    #            cat_files.append(path+'/'+i+'/'+os.listdir(path+'/'+i)[0])
    #        else:
    #            fails.append(path+'/'+i)
    #out = open(o.out,'w')
    #fl = len(cat_files)
    #lc = 0
    #for f in cat_files:
    #    os.system("cat %s >> %s"%(f,o.out))
    #    lc+=1
    #    if o.verbose==True:
    #        sys.stderr.write('\r'+'' *0)
    #        sys.stderr.write(str(int(lc*100/fl))+'%')
    #        sys.stdout.flush()
    #if len(fails) != 0:sys.stderr.write("\nNot concatenated!:\n%s"%fails)
    #sys.stderr.write("\n%s files Successfully concatenated!\n"%lc)

if __name__ == "__main__": main()

