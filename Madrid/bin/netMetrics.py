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
import itertools as it
import numpy as np
import pandas as pd
import cPickle as pickle
from optparse import OptionParser
import networkx as nx
import matplotlib.pyplot as plt

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
parser.add_option("-i","--input",type="string",action="store",dest="tensor",default='OTUp_SO_females',
                      help="Input file. Pandas dataframe (numpy based matrix)")
parser.add_option("-a","--atrib",type="string",action="store",dest="atrib",default='SO_females',
                      help="SOsin matrix")
parser.add_option("-o","--out",type="string",action="store",dest="out",default='OTUp_females',
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
    th = 0.05
    j=0
    for i in files:
        tensor = pickle.load(open(i, 'r'))
        #so = pickle.load(open(so_files[j], 'r'))
        ft = {} # filtered tensor with all otus
        st = {} # filtered tensor with filtered otus
        q=0
        for k in tensor.keys():
            #tensor[k] = tensor[k][tensor[k]>0.1].fillna(0)  ## threshold to 0.1 for links removement
            df=tensor[k]
            temp = []
            permut = [x for x in it.combinations(df.index.values.tolist(),2)]   # All comparisons tuples (Oi,Oj)
            m = pd.DataFrame(index=df.index.values.tolist(), columns=df.index.values.tolist()).fillna(0)  # Empty new matrix of network
            for i in permut:temp.append([i[0],i[1],df.ix[i[0],i[1]]])
            temp = pd.DataFrame(np.array(temp)).sort(2,ascending=False).reset_index().drop('index', axis=1)
            #for i in range(int(0.05*len(permut))):m.ix[temp.ix[i,][0],temp.ix[i,][1]] = float(temp.ix[i,][2])
            rows=[]
            for i in range(int(th*len(permut))):rows.append([temp.ix[i,][0],temp.ix[i,][1],float(temp.ix[i,][2])])
            st[k] = rows
            rows = pd.DataFrame(np.array(rows)).sort(2,ascending=False).reset_index().drop('index', axis=1)
            rows.to_csv("%s_%s_best5p_net.txt"%(out[j],k),sep='\t',index=False)
            #print m.ix[0:5,0:5]
            #m = m.add(m.transpose())
            #m.to_csv("%s_%s_best5p.txt"%(out[j],k),sep='\t')
            #ft[k] = m
        nets = netmetrics(st)
        nets.wdegree().to_csv("%s_best5p_wdegree_st.txt"%out[j],sep="\t")
        #nets.aveclustercoef().to_csv("%s_avclustcoef.txt"%out[j],sep="\t")
        #nets.eigvectc().to_csv("%s_eigvectc.txt"%out[j],sep="\t")
        #nets.betwcent().to_csv("%s_betwcent.txt"%out[j],sep="\t")
        #nets.aveshortest().to_csv("%s_avshortest.txt"%out[j],sep="\t")
        #nets.degassortcoef().to_csv("%s_assortcoef.txt"%out[j],sep="\t")
        #nets.closeness().to_csv("%s_closeness.txt"%out[j],sep="\t")
        #nets.graph()
        j+=1

    sys.stderr.write('All Done!')

class netmetrics():
    """ Calculates some network metrics.
    """
    def __init__(self,tensor):
        assert type(tensor) is dict
        self.nets = tensor
        if type(tensor[tensor.keys()[0]]) != list:
            self.pd = True
            self.nodes = list(tensor[tensor.keys()[0]].columns.values)
            self.indiv = tensor.keys()
        else:
            self.pd = False
            self.indiv = tensor.keys()

    def wdegree(self):
        '''
        It returns weighted degree dictionary of the nodes for each sample
        '''
        if self.pd :
            deg=pd.DataFrame(np.nan,index=self.indiv,columns=self.nodes)
            for k in self.indiv:
                G = nx.from_numpy_matrix(self.nets[k].values)
                G.edges(data=True)
                wdegree = G.degree(weight='weight')
                wd = [wdegree[node] for node in wdegree]
                deg.loc[k] = wd
        else:
            deg = pd.DataFrame()
            id = {}
            for k in self.nets.keys():
                G= nx.DiGraph()
                G.add_weighted_edges_from(self.nets[k])
                wdegree = G.degree(weight='weight')
                #wd = node[wdegree[node] for node in wdegree]
                wd={}
                for node in wdegree:wd[node]=wdegree[node]
                deg = deg.append(wd,ignore_index=True)
                #print deg.index[-1]
                id[deg.index[-1]]=k
            deg.rename(index=id,inplace=True)
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
        It returns average shortest path from each network. needs full conected network
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
    def graph(self):
        #print so.ix[:5,:5]
        for k in self.indiv:
            G = nx.from_numpy_matrix(self.nets[k].values)
            G.edges(data=True)
            #print G.nodes()
            #print G.edges(data=True)
            #for i in G.nodes():
            #    print G.degree(i,weight='weight')
            #print G.degree(weight='weight')
            nx.draw(G)
            plt.show()

if __name__ == "__main__": main()