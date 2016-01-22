#!/usr/bin/python
# Guillermo Torres MSc,
# g.torres@ikmb.uni-kiel.de
# Genetics of human longevity group
# Institute of Clinical Molecular Biology
# last update: August 2015
# **16sTensor.py integrates: mothur 16s output + Network analysis
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
from cStringIO import StringIO

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to Bla Bla.\n REQUIREMENTS:\n** Bla\n** bla \n
\t%prog --prb PROBES_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
By Guillermo G. Torres (ggtorrese@unal.edu.co)
- IKMB - Christian-Albrechts-Universitat zu Kiel"""

parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--input",type="string",action="store",dest="inS",default='/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome/Madrid/data/MiSeq_SOP/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.pick.tx.1.cons.tax.summary',
                      help="Input file")
parser.add_option("-z","--lib_size",type="string",action="store",dest="libS",default='/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome/Madrid/data/MiSeq_SOP/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.an.unique_list.count.summary',
                      help="mothur .an.unique_list.count.summary file")
parser.add_option("-S", action="store_true",dest="mothur",default=True,
                      help = "True if input is '.cons.tax.summary' from mothur - default False")
parser.add_option("-d","--design",type="string",action="store",dest="design",default='/Users/guillermotorres/Documents/Proyectos/Doctorado/Metagenome/Madrid/data/MiSeq_SOP/mose.age.design.txt',
                      help="extra info of samples 'age'")
parser.add_option("-t","--tax_lev",type="string",action="store",dest="tl",default='6',
                      help="taxonomic_level")
parser.add_option("-o","--out",type="string",action="store",dest="out",default='tensor',
                      help="out files name - default: infile.out")
parser.add_option("-g","--graph",type="string",action="store",dest="img",#default='',#L1_1XDgR',#nirsdb',
                      help="Reference data base to evaluate the specificity of gene regions from which will be designed the probes")
parser.add_option("-n","--processors_num",type="int",action="store",dest="proc",metavar=" Number of processor used by BLAST (int)",default=1,
                      help="Number of processor used by BLAST - default value = 1")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.inS == None):# and o.reads_db == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
path = os.getcwd()


def main():
    if o.mothur==True:
        assert o.libS is not None, "\n Library size file missing! (.an.unique_list.count.summary file)"
        assert o.design is not None, "\n Design file missing!"
        tensor = org_corrMatrix(retrieve_matrix(o.inS,o.libS,o.tl,o.design))    # dictionary tensor, store all matrixes
        pickle.dump(tensor,open(o.out,'w'), pickle.HIGHEST_PROTOCOL)
        for k in tensor.keys():
                matrix_image(tensor[k])
                break
        if o.img == True:
            for k in tensor.keys():
                matrix_image(tensor[k])
    else: org_corrMatrix(o.inS)
    print 'All done!'

def matrix_image(corrMat):
    fig, ax = plt.subplots()
    im = ax.imshow(corrMat,interpolation='None',cmap=plt.get_cmap('gray_r'),vmin=0, vmax=1)
    fig.colorbar(im)
    #ax.set_title('network')
    ax.plot()
    plt.show()

def retrieve_matrix(inf,libs,tl,design):
    '''Input: '.cons.tax.summary' file from mothur output, This function retrieve the input matrix
        at some taxonomic level (tl); by default at level of 6 -> gender
    '''
    matrix, lib_size = [],[]
    lib_factor = {}
    for line in open(inf,'r'):
        line = line.strip("\n").split("\t")
        if line[0] == "taxlevel":colnames=line
        elif line[0] == tl:
            line[2] = line[2].strip("\"").split(" ")[0]
            matrix.append(line)
    df = pd.DataFrame(matrix, columns=colnames)

    for i in open(libs,'r'):
        lib_size.append(float((i.strip("\n").split("\t")[1])))
        lib_factor[i.strip("\n").split("\t")[0]] = [1,float((i.strip("\n").split("\t")[1]))]
    for k in lib_factor.keys():
        lib_factor[k][0] = round(lib_factor[k][1]/min(lib_size),4)
        df[k] = df[k].apply(lambda x: round(float(x)/float(lib_factor[k][0]),4))

    df = df.ix[:,[2]+range(5,len(colnames))]

    df.loc['age'] = ['NA' for n in range(len(df.columns))]
    for line in open(design,'r'): df.loc['age',line.strip("\n").split("\t")[0]] = line.strip("\n").split("\t")[1]

    df = df.transpose()
    df.columns = map(list, df.loc[["taxon"]].values)[0]
    df = df.sort('NA')
    df = df.drop(df.index[[-1]])
    df.pop('NA')
    df.pop('unclassified')
    return df.transpose()

def org_corrMatrix(df):
    '''Input: Matrix i,j {organisms,individuals}'''
    tensor = {}
    windowsize = 4
    ti = 0
    while ti+windowsize != len(df.columns)+1:     # run an overlapping window
        tensor[ti] = corrByOrg(org_window(df.ix[:,ti:ti+windowsize]))
        ti+=1
    return tensor

def org_window(win):
    w_norm = win.copy()
    for i in range(len(win.columns)):
        ind =  w_norm.ix[:,i]
        w_norm.ix[:,i] = ind.apply(lambda x: round((x - ind.min())/(ind.max() - ind.min()),4) )
    return w_norm

def corrByOrg(df):
    '''
    Correlation by organism distance
    corrOij = |Oi - Oj|
    '''

    permut = [x for x in it.combinations(df.index.values.tolist(),2)]
    corrMat = pd.DataFrame(index=df.index.values.tolist(), columns=df.index.values.tolist()).fillna(0)
    exp = [x for x in it.combinations([1,2,3],2)]
    for i in permut:
        org_v0 = []
        org_v1 = []
        df.loc[i[0]].apply(lambda x: org_v0.append(x))
        df.loc[i[1]].apply(lambda x: org_v1.append(x))
        org_diff = [abs(org_v0[x] - org_v1[x]) for x in range(len(org_v0))]
        if max(org_diff) == 0 : org_assoc_weight = float(1)
        else: org_assoc_weight = (max(org_diff) - np.mean(org_diff))/max(org_diff)
        corrMat.ix[i[0],i[1]] = org_assoc_weight
    corrMat = corrMat.add(corrMat.transpose())
    return corrMat


if __name__ == "__main__": main()