#!/usr/bin/python
# Guillermo Torres MSc,
# ggtorrese@unal.edu.co, guigotoe@gmail.com
# Bioinformatic Research Group
# Biotechnology Institute of National University of Colombia
# last update: November 2013
# **HISS integrates: Comparison + Section processes
# This is written as part of IN SILCO MICROARRAY pipeline, but could be
# splitted to serve different porpuoses.
#/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Requirements:-Biopython >= v1.52; -New Blast installed (BLAST+)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

import re, sys, os, math 
from tempfile import mkstemp
from collections import Counter
import subprocess as sub
from optparse import OptionParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIStandalone
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbiblastnCommandline


## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to get hit blast. The Results are presented in tabular file.\n REQUIREMENTS:\n** Biopython >= v 1.52\n** BLAST \n
\t%prog --prb PROBES_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
Part of in silico Hybridization system by Guillermo G. Torres (ggtorrese@unal.edu.co) - Bioinformatics
Group - IBUN - National University of Colombia"""
parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-b","--blast_file",type="string",action="store",dest="blast",#default='blast.out',
                      help="Input blast out file in default format")
parser.add_option("-o","--out_file",type="string",action="store",dest="out",default='hits.txt',
                      help="Out tabular files name - default name = hits.txt")
parser.add_option("-i","--gene_list",type="string",action="store",dest="genes",#default='subcluster_18_log2_medianCentered_fpkm.matrix',#L1_1XDgR',#nirsdb',
                      help="gene list")
parser.add_option("-v", action="store_true",dest="verbose",default=True,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.blast == None and o.genes == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")
   
## temp raw probes file ##
path = os.getcwd()
out = open(o.out,'w')

def main():

    handle = open(o.genes, "rU")
    fl = int(open(o.genes,"rU").read().count("\n")) # length of file
    lc=0 # line counter
    if o.verbose==True:sys.stderr.write("""
Your parameters:\ninput gene list = %s\nblast file = %s\n
Your out file name = %s\nReading gene list...\n"""%(o.genes,o.blast,o.out))
    gene = {}
    try:
        flag = 0
        for i in handle:
            if flag != 0:
                gene[i.split("\t")[0]] = i.split("\t")[0]
            flag +=1
            ## Progress counter for verbose mode ##
            lc+=1
            if o.verbose==True:
                sys.stderr.write('\r'+'' *0)
                sys.stderr.write(str(int(lc*100/fl))+'%')
                sys.stdout.flush()
    except TypeError:sys.stderr.write('Error!... in gene list file')

    getHits(gene)
    sys.stderr.write("\nSearching done!\n")
   
def getHits(gene):
    ''' BLAST parser using Biopython
    Input: name of blast out file in standard ouput format
    Outputs: 2 files 
    ''' 
    inf = open(o.blast,'rU')
    parser = NCBIStandalone.BlastParser()
    error_parser = NCBIStandalone.BlastErrorParser(inf)
    iterator = NCBIStandalone.Iterator(inf, error_parser)
    err_iterator = NCBIStandalone.Iterator(inf, error_parser)
    #next_record =
    
    ## *** Parsing *** ##
    lg = len(gene)
    if o.verbose==True:
        sys.stderr.write("\nGetting hits...\n")
    for record in iterator:
        query = record.query.split(" ")[0]

        if query in gene:
            out.write("%s\n"%gene[query])
            if record.alignments is []:
                out.write("%s\tNA\tNA\tNA\n"%gene[query])
            else:
                flag = 0
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        #-->## ** Selection Process **##
                        if float(hsp.expect) < 0.0001 and flag < 3:
                            out.write("%s\t%s\t%s\tHigh\n"%(gene[query], alignment.title.split(">")[1],float(hsp.expect)))
                            flag += 1
                        elif float(hsp.expect) < 1.0 and flag < 3:
                            out.write("%s\t%s\t%s\tLow\n"%(gene[query], alignment.title.split(">")[1],float(hsp.expect)))
                            flag += 1
                        elif float(hsp.expect) < 5.0 and flag < 3:
                            out.write("%s\t%s\t%s\tScare\n"%(gene[query], alignment.title.split(">")[1],float(hsp.expect)))
                            flag += 1
                        elif float(hsp.expect) > 1.0 and flag < 1:
                            out.write("%s\tNA\tNA\tNA\n"%gene[query])
                            flag += 1
            del gene[query]
            if o.verbose==True:
                sys.stderr.write('\r'+'' *0)
                sys.stderr.write(str(int((lg-len(gene))*100/lg))+'%')
                sys.stdout.flush()
        else: pass

    if (lg-len(gene)) != len(gene):
        sys.stderr.write ("\nGenes not found:\n%s"%gene.keys())
    
if __name__ == "__main__": main()