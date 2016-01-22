#!/usr/bin/python
# Guillermo Torres MSc,
# g.torres@ikmb.uni-kiel.de
# Genetics of human longevity group
# Institute of Clinical Molecular Biology
# last update: August 2015
# **netMetrics.py integrates: mothur Networkx for Network analysis
# This is written as part of Metagenome analysis pipeline, but could be
# splitted to serve different porpuoses.
#/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Requirements:
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
__author__ = 'Guillermo G. Torres (PhD)'

import re, sys, os, math
import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt
import cPickle as pickle
from optparse import OptionParser
import networkx as nx

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to Bla Bla.\n REQUIREMENTS:\n** Bla\n** bla \n
\t%prog --prb PROBES_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
By Guillermo G. Torres (ggtorrese@unal.edu.co)
- IKMB - Christian-Albrechts-Universitat zu Kiel"""

parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--input",type="string",action="store",dest="tensor",default='tensor',
                      help="Input file")
parser.add_option("-z","--lib_size",type="string",action="store",dest="libS",default='/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome/Madrid/data/MiSeq_SOP/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.count.summary',
                      help="mothur .an.unique_list.count.summary file")
parser.add_option("-S", action="store_true",dest="mothur",default=True,
                      help = "True if input is '.cons.tax.summary' from mothur - default False")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.tensor == None):# and o.reads_db == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
path = os.getcwd()


def main():
    tensor = pickle.load(open(o.tensor, 'r'))
    for k in tensor.keys():
        metrics(tensor[k])
        break

    print "All done!"
def metrics(matrix):
    G = nx.from_numpy_matrix(matrix.values)
    print nx.all_pairs_dijkstra_path(G)

if __name__ == "__main__": main()