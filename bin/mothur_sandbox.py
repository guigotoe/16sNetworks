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
\t%prog -f sam/bam_file -o outname [options] arg \n\t%prog {-h --help}\n
\t It requires pysam libray
As part of Toolbox scripts by Guillermo G. Torres (ggtorrese@unal.edu.co) -
IKMB - Christian-Albrechts-Universitat zu Kiel (CAU) """
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-f","--fasta",type="string",action="store",dest="fasta",default='./mothur/16s.file_qc.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta',
                      help="Input folder path - absolute path is recomended")
parser.add_option("-n","--name",type="string",action="store",dest="name",default='./mothur/16s.file_qc.trim.contigs.good.unique.good.filter.names',
                      help="Input folder path - absolute path is recomended")
parser.add_option("-g","--group",type="string",action="store",dest="group",default='./mothur/16s.file_qc.contigs.groups',
                      help="Input folder path - absolute path is recomended")
parser.add_option("-o","--out_file",type="string",action="store",dest="out",default='new',
                      help="Out files name - default name = cat.fna")
parser.add_option("-D","--diamond", action="store_true",dest="diamond",default=True,
                      help = "Diamond SAM file - default False")
parser.add_option("-v", action="store_true",dest="verbose",default=True,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.fasta == None):
    parser = OptionParser(usage="\n\t%prog -f sam/bam_file -o outname [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
#path = os.path.abspath(o.inf)
#out = open(o.out,'w')

def main():
    out_n=open("%s.names"%o.out,'w')
    out_g=open("%s.groups"%o.out,'w')
    ids = []
    names,groups={},{}
    notf =[]
    notf_g =[]
    newname,newgroup= {},{}
    for j in open(o.name,'r'):
        nline = j.strip("\n").split("\t")
        names[nline[0]]=nline[1]
    for i in open(o.fasta,'r'):
        if re.findall('>',i):
            id = i.strip(">\n")
            ids.append(id)
            try:
                if names[id]:
                    newname[id] = names[id]
            except KeyError:
                notf.append(id)

    for q in open(o.group,'r'):
        gline = q.strip("\n").split("\t")
        groups[gline[0]]=gline[1]
    for k in newname.keys():
        out_n.write("%s\t%s\n"%(k,newname[k]))
        seqs = newname[k].split(',')
        for i in range(len(seqs)):
            try:
                if groups[seqs[i]]:
                    newgroup[seqs[i]]=groups[seqs[i]]
            except KeyError:
                notf_g.append(seqs[i])
    print len(newgroup.keys())
    for r in newgroup.keys():
        out_g.write("%s\t%s\n"%(r,newgroup[r]))


    print 'Done!'

if __name__ == "__main__": main()

