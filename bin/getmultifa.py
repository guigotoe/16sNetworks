#!/usr/bin/python
####################################################
# By Guillermo Torres PhD.c - IKMB                 #
# Institue of Clinical Molecular Biology (IKMB)    #
# Christian-Albrechts-Universitat zu Kiel (CAU)    #
####################################################
# Last update: September 2015

import re, sys, os, gzip
from optparse import OptionParser

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to concatenate fna files from input_folder folder.
\t%prog -f allprokaryote_folder -o outname [options] arg \n\t%prog {-h --help}\n
As part of Toolbox scripts by Guillermo G. Torres (ggtorrese@unal.edu.co) -
IKMB - Christian-Albrechts-Universitat zu Kiel (CAU) """
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-f","--folder",type="string",action="store",dest="folder",default='allprok_test/dfg',
                      help="Input folder path - absolute path is recomended")
parser.add_option("-o","--out_file",type="string",action="store",dest="out",default='cat.fna',
                      help="Out files name - default name = cat.fna")
parser.add_option("-z","--gz",action="store_true",dest="gz",default=False,
                      help="input file in .gz format")
parser.add_option("-t","--targz",action="store_true",dest="tgz",default=False,
                      help="input file in .tar.gz format")
parser.add_option("-v", action="store_true",dest="verbose",default=True,
                      help = "Verbose")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.folder == None):
    parser = OptionParser(usage="\n\t%prog -f input_folder -o outname [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")
   
## temp raw file ##
path = os.path.abspath(o.folder)
out = open(o.out,'w')

def main():
    files = []
    for i in os.listdir(path):
        if o.gz == True:
            if re.findall("\.gz$",i): files.append(i)
            ext = ".gz"
        elif o.tgz == True:
            if re.findall("\.tar\.gz$",i): files.append(i)
            ext = ".tar.gz"
        else:
            if re.findall("\.(fna)|(fa)|(fasta)$",i): files.append(i)
            ext = ".fasta"
    if len(files) == 0: print ("No files with extension '%s' found"%(ext))
    else:
        sys.stderr.write("Concatenating...\n")
        fl = len(files)
        lc = 0
        for j in files:
            if o.gz == True:
                with gzip.open("%s/%s"%(path,j),mode="rb") as f:
                    for l in f: out.write(l)
            else:
                for l in open("%s/%s"%(path,j),'r'):
                    print l
                    out.write(l)
            lc+=1
            if o.verbose==True:
                sys.stderr.write('\r'+'' *0)
                sys.stderr.write(str(int(lc*100/fl))+'%')
                sys.stdout.flush()

if __name__ == "__main__": main()

