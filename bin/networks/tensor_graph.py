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
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
import itertools as it
import matplotlib.pyplot as plt
import cPickle as pickle
from optparse import OptionParser
from cStringIO import StringIO

pd.set_option('display.height', 1000)
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

## PARSING... ##
 #* Options *#
parser = OptionParser()
usage = """\nThis scripts is able to Bla Bla.\n REQUIREMENTS:\n** Bla\n** bla \n
\t%prog -i infile -z libs [options] arg \n\t%prog {-h --help}\n
By Guillermo G. Torres (ggtorrese@unal.edu.co)
- IKMB - Christian-Albrechts-Universitat zu Kiel"""

parser = OptionParser(usage=usage,version="%prog 1.0")
parser.add_option("-i","--tensors",type="string",action="store",dest="tens",default='tensor_all,tensor_females,tensor_males',
                      help="mothur .cons.tax.summary file")
parser.add_option("-v", action="store_true",dest="verbose",default=False,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.tens == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
path = os.getcwd()

def main():
    tens = o.tens.split(',')
    for t in tens:
        print t
        tensor = pickle.load(open(t, 'r'))
        matrix_image(tensor,t)
    print 'All done!'

def matrix_image(tensor,name):
    fl=1
    #f, ((ax1,ax2,ax3,ax4,ax5),(ax6,ax7,ax8,ax9,ax10)) = plt.subplot(2,2,sharex='col',shrey='row')
    fig = plt.figure()
    print len(tensor.keys())
    counter = 0
    for k in OrderedDict(sorted(tensor.items())):
        if fl == 1:
            ax = fig.add_subplot(3,4,fl)
            im = ax.imshow(tensor[k],interpolation='None',cmap=plt.get_cmap('RdYlBu'),vmin=0, vmax=1)
            plt.title(k)
            ax.tick_params(labelsize=7)
        elif fl > 1:
            ax = fig.add_subplot(3,4,fl,sharex=ax, sharey=ax)
            im = ax.imshow(tensor[k],interpolation='None',cmap=plt.get_cmap('RdYlBu'),vmin=0, vmax=1)
            #fig.colorbar(im)
            plt.title(k)
            ax.tick_params(labelsize=7)
        if fl == 12 :
            plt.subplots_adjust(left=0.05, right=0.90, top=0.95, bottom=0.05,wspace=0.2,hspace=0.3)
            plt.colorbar(im,cax=plt.axes([0.92, 0.1, 0.03, 0.8]))
            #plt.show()
            fl = 0
            counter +=1
            pdf = PdfPages("/home/torres/Documents/Projects/Metagenome/resutls/all/Networks/%s_%s.pdf"%(name,counter))
            pdf.savefig()  # saves the current figure into a pdf page
            pdf.close()
            plt.close()
            fig = plt.figure()
        fl +=1
    if fl != 0:
        plt.subplots_adjust(left=0.05, right=0.90, top=0.95, bottom=0.05,wspace=0.2,hspace=0.3)
        plt.colorbar(im,cax=plt.axes([0.92, 0.1, 0.03, 0.8]))
        #plt.show()
        pdf = PdfPages("/home/torres/Documents/Projects/Metagenome/resutls/all/Networks/%s_%s.pdf"%(name,counter+1))
        pdf.savefig()  # saves the current figure into a pdf page
        pdf.close()
        plt.close()


if __name__ == "__main__": main()



        #im = ax.imshow(tensor[k],interpolation='None',cmap=plt.get_cmap('RdYlBu'),vmin=0, vmax=1)

        #fig = plt.subplot(2,5,fl+1)
        #ax = fig.subplot()
        #fig, ax = plt.subplot(10,10,fl+1)
        #fig, ax = plt.subplots()
        #im = ax.imshow(tensor[k],interpolation='None',cmap=plt.get_cmap('RdYlBu'),vmin=0, vmax=1)
        #im = imshow(tensor[k],interpolation='None',cmap=plt.get_cmap('RdYlBu'),vmin=0, vmax=1)
        #fig.colorbar(im)
        #ax.set_title('network')
        #ax.plot()
        #plt.title(k)
        #plt.plot(im)
        #plt.show()
