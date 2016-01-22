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

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to Bla Bla.\n REQUIREMENTS:\n** Bla\n** bla \n
\t%prog --prb PROBES_fileNAME --refDB path/DB_NAME [options] arg \n\t%prog {-h --help}\n
By Guillermo G. Torres (ggtorrese@unal.edu.co)
- IKMB - Christian-Albrechts-Universitat zu Kiel"""

parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--input",type="string",action="store",dest="tensor",default='OTUp_SO_males',
                      help="Input file. Pandas dataframe (numpy based matrix)")
parser.add_option("-a","--atrib",type="string",action="store",dest="atrib",default='SO_males',
                      help="SOsin matrix")
parser.add_option("-o","--out",type="string",action="store",dest="out",default='OTUp_males',
                      help="SOsin matrix")
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

    files = o.tensor.split(',')
    so_files = o.atrib.split(',')
    out = o.out.split(',')

    j=0
    for i in files:
        tensor = pickle.load(open(i, 'r'))
        so = pickle.load(open(so_files[j], 'r'))
        for k in tensor.keys():
            tensor[k] = tensor[k][tensor[k]>0.1].fillna(0)  ## threshold to 0.1 for links removement
            #tensor[k].to_csv("%s_%s_t01.txt"%(out[j],k),sep='\t')

        nets = netmetrics(tensor)
        #nets.wdegree().to_csv("%s_t01_wdegree.txt"%out[j],sep="\t")
        #nets.aveclustercoef().to_csv("%s_avclustcoef.txt"%out[j],sep="\t")
        #nets.eigvectc().to_csv("%s_eigvectc.txt"%out[j],sep="\t")
        #nets.betwcent().to_csv("%s_betwcent.txt"%out[j],sep="\t")
        #nets.aveshortest().to_csv("%s_avshortest.txt"%out[j],sep="\t")
        #nets.degassortcoef().to_csv("%s_assortcoef.txt"%out[j],sep="\t")
        #nets.closeness().to_csv("%s_closeness.txt"%out[j],sep="\t")
        nets.netsout([33.3,57.5])
        j+=1

    sys.stderr.write('All Done!')

class netmetrics():
    """ Calculates some network metrics.
    """
    def __init__(self,tensor):
        assert type(tensor) is dict
        self.nets = tensor
        self.nodes = list(tensor[tensor.keys()[0]].columns.values)
        self.indiv = tensor.keys()
    def wdegree(self):
        '''
        It returns weighted degree dictionary of the nodes for each sample
        '''
        deg=pd.DataFrame(np.nan,index=self.indiv,columns=self.nodes)
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            wdegree = G.degree(weight='weight')
            wd = [wdegree[node] for node in wdegree]
            deg.loc[k] = wd
        return deg
    def aveclustercoef(self):
        '''
        It returns average cluster coeficient from each network
        '''
        acc=pd.DataFrame(np.nan,index=self.indiv,columns=['ave_cluster_coef'])
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            acc.loc[k] = nx.average_clustering(G,weight='weight')
        return acc
    def eigvectc(self):
        '''
        It returns eigenvector centrality from each network
        '''
        evc=pd.DataFrame(np.nan,index=self.indiv,columns=self.nodes)
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            centrality = nx.eigenvector_centrality_numpy(G,weight='weight')
            nevc = [centrality[node] for node in centrality]
            evc.loc[k] = nevc
        return evc
    def betwcent(self):
        '''
        It returns betweenness centrality from each network.
        *For weighted graphs the edge weights must be greater than zero.
        Zero edge weights can produce an infinite number of equal length paths between pairs of nodes.
        '''
        bc=pd.DataFrame(np.nan,index=self.indiv,columns=self.nodes)
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            betw = nx.betweenness_centrality(G,weight='weight')
            bt = [betw[node] for node in betw]
            bc.loc[k] = bt
        return bc
    def closeness(self):
        '''
        It returns closeness centrality from each network.
        '''
        cc=pd.DataFrame(np.nan,index=self.indiv,columns=self.nodes)
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            closeness = nx.closeness_centrality(G,distance='weight')
            c = [closeness[node] for node in closeness]
            cc.loc[k] = c
        return cc
    def aveshortest(self):
        '''
        It returns average shortest path from each network
        '''
        asp=pd.DataFrame(np.nan,index=self.indiv,columns=['ave_shortest_path'])
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            asp.loc[k] = nx.average_shortest_path_length(G,weight='weight')
        return asp
    def degassortcoef(self):
        '''
        It returns degree assortativity coefficient path from each network
        '''
        dac=pd.DataFrame(np.nan,index=self.indiv,columns=['ave_cluster_coef'])
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            dac.loc[k] = nx.degree_pearson_correlation_coefficient(G,weight='weight')
            print nx.degree_pearson_correlation_coefficient(G,weight='weight')
        return dac
    def netsout(self,ids):
        '''
        It returns degree assortativity coefficient path from each network
        '''
        if ids == 'all':print 'alls'
        else:
            for id in ids:
                print self.nets[id].ix[0:5,0:5]
                df = self.nets[id][self.nets[id]==0].fillna()

                #G = nx.from_numpy_matrix(self.nets[id].values)
                #G.edges(data=True)
                #print G.edges
                sys.exit('stop')
                #out=pd.DataFrame(np.nan,index=self.indiv,columns=['ave_cluster_coef',])
                G = nx.from_numpy_matrix(self.nets[i].values)
                #G.edges(data=True)
                #dac.loc[k] = nx.degree_pearson_correlation_coefficient(G,weight='weight')
                #print nx.degree_pearson_correlation_coefficient(G,weight='weight')
                #out.to_csv("%s_closeness.txt"%out[j],sep="\t")



def metrics(matrix):
    print so.ix[:5,:5]
    G = nx.from_numpy_matrix(matrix.values)
    G.edges(data=True)
    #print G.nodes()
    #print G.edges(data=True)
    #for i in G.nodes():
    #    print G.degree(i,weight='weight')
    print G.degree(weight='weight')
    #nx.draw(G)
    #plt.show()

if __name__ == "__main__": main()