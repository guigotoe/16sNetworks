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
import cPickle as pickle
from optparse import OptionParser
import networkx as nx

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

path = "/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome/networks"#os.getcwd()
path2 = os.getcwd()

def main():
    df = pickle.load(open("%s/SO_males"%path2, 'r'))
    print df.ix[0:5,0:5]



def main2():
    inf = "%s/net57.5.txt"%path

    file = {}
    c = 0
    for i in open(inf):
        i = i.strip("\n").split("\t")
        if c == 0:file[(i[0],i[1])] = i
        else:
            if (i[1],i[0]) in file.keys():pass
            else:file[(i[0],i[1])] = i
        c+=1
        #if c == 20: break
    outf = open("%s/net57_5.txt"%path,'w')
    print len(file)
    for k in file.keys():
        outf.write ("%s\t%s\t%s\n"%(file[k][0],file[k][1],file[k][2]))


if __name__ == "__main__": main()