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
import scipy.stats  as stats

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
parser.add_option("-i","--taxonomy_sum",type="string",action="store",dest="inS",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/16S.an.cons.taxonomy.summary',
                      help="mothur .cons.tax.summary file")
parser.add_option("-u","--otu_taxonomy",type="string",action="store",dest="inT",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/16S.an.0.03.cons.taxonomy',
                      help="mothur .cons.tax file")
parser.add_option("-n","--norm_otu_reads",type="string",action="store",dest="inR",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/16S.an.0.03.subsample.shared',
                      help="mothur .an.shared file")
parser.add_option("-a","--otu_reads",type="string",action="store",dest="inRa",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/16S.an.shared',
                      help="mothur .an.shared file")
parser.add_option("-z","--lib_size",type="string",action="store",dest="libS",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/16S.an.count.summary',
                      help="mothur .count.summary file")
parser.add_option("-S", action="store_true",dest="mothur",default=True,
                      help = "True if input is '.cons.tax.summary' from mothur - default False")
parser.add_option("-d","--design",type="string",action="store",dest="design",default='/Users/guillermotorres/Dropbox/Doctorado/Microbiome/Networks/design.txt',
                      help="extra info of samples 'age'")
parser.add_option("-l","--tax_lev",type="string",action="store",dest="tl",default='6',
                      help="taxonomic_level")
parser.add_option("-o","--out",type="string",action="store",dest="out",default='tensor',
                      help="out files name - default: infile.out")
parser.add_option("-g","--graph",action="store_true",dest="img",default=True,#L1_1XDgR',#nirsdb',
                      help="Graph tensor")
parser.add_option("-c","--scaling",action="store_true",dest="scaling",default=False,#L1_1XDgR',#nirsdb',
                      help="Scaling the data")
parser.add_option("-w","--window_size",type="int",action="store",dest="win",default=10,
                      help="Number of individuals in the window that computes mean idividual - default 0.15")
parser.add_option("-t","--taxonomy_threshold",action="store",type="float",dest="taxt",default=0.15,
                      help="""If is float, then taxon should be present in as many individuals as
the percentage (float) of the total batch. If is int, then taxon should be present in as many individuals as the int value (int >= 10).
""")
parser.add_option("-v", action="store_true",dest="verbose",default=True,
                      help = "Verbose - default False")
(o,args) = parser.parse_args()
if len(args)!=0 or (o.inS == None) and (o.inT == None) and (o.inR == None):
    parser = OptionParser(usage="\n\t%prog -p path/IN_ProbefileNAME -r path/DB_NAME [options] arg \n\t%prog {-h --help}")
    parser.error("incorrect number of arguments\n")

## temp raw probes file ##
path = os.getcwd()
ids = {} # ['5090939SPC', 'female', '81']
totu, rotu = {},{}

def main():
    if o.mothur==True:

        sys.stderr.write("Reading files: \n")
        assert o.libS is not None, "\n Library file is missing! (.count.summary file)"
        assert o.design is not None, "\n Design file is missing!"

        sys.stderr.write("Preparing files: \n")
        matrices = retrieve_matrix(o.inS,o.libS,o.tl,o.design)  # [allindiv,males,females]
        sufix = ['all','males','females']

        sys.stderr.write("Tensor calculation: \n")

        #for i in range(len(matrices)):
        #    pickle.dump(matrices[i],open("SO_%s"%(sufix[i]),'w'), pickle.HIGHEST_PROTOCOL) # Exporting tensor dictionary into binary file to use it afterwards
        #    tensor = org_corrMatrix(matrices[i],i+1)
        #    pickle.dump(tensor,open("%s_%s"%(o.out,sufix[i]),'w'), pickle.HIGHEST_PROTOCOL) # Exporting tensor dictionary into binary file to use it afterwards

        sys.stderr.write("OTU profile calculation: \n")

        for j in ["SO_males","SO_females"]:
            if os.path.exists(j) == False: print "no esta" #df = otu_reads(df,taxt)
            else: df = pickle.load(open(j, 'r'))
            k = dfs_log(df)
            tensor = OTU_profile(k)
            pickle.dump(tensor,open("OTUp_%s"%(j),'w'), pickle.HIGHEST_PROTOCOL) # Exporting tensor dictionary into binary file to use it afterwards
        #for j in ["SO_males","SO_females","SO_all"]:
        #    if os.path.exists(j) == False: print "no esta" #df = otu_reads(df,taxt)
        #    else: df = pickle.load(open(j, 'r'))
        #    k = dfs_log(df)
        #    tensor = OTU_profile_2(k)
        #    pickle.dump(tensor,open("OTUp_2_%s"%(j),'w'), pickle.HIGHEST_PROTOCOL) # Exporting tensor dictionary into binary file to use it afterwards

    else: org_corrMatrix(o.inS)

    sys.stderr.write("\nTensors were successfuly generated")

def OTU_profile(df):
    '''
    It generates OTU abndance profile
    :return:
    '''

    permut = [x for x in it.combinations(df.index.values.tolist(),2)]   # All comparisons tuples (Oi,Oj)
    tensor = {}
    if o.win < 9: o.win = 10
    windowsize = o.win
    fl= 1+len(df.columns)-o.win
    s = 0
    print len(permut)
    print fl

    while s+windowsize != len(df.columns)+1:     # run an overlapping window
        corrMat = pd.DataFrame(index=df.index.values.tolist(), columns=df.index.values.tolist()).fillna(0)  # Empty comparison matrix
        for i in permut:
            frames = df.ix[(i[0],i[1]),s:s+windowsize]
            if sum(frames.ix[0,]) == 0 : cosd = 0.0
            elif sum(frames.ix[1,]) == 0 : cosd = 0.0
            else: cosd = np.dot(frames.ix[0,:],frames.ix[1,:])/(np.linalg.norm(frames.ix[0,:])*np.linalg.norm(frames.ix[1,:])) # cosin distance
            corrMat.ix[i[0],i[1]] = cosd
            #f = frames.transpose()
            #corrMat.ix[i[0],i[1]] = frames.corr(method='pearson').ix[0,1]  # pearson distance but it needed frames.transpose()
        x = []
        for q in frames.columns:x.append(int(q.split('_')[0]))  # mean age of the window
        M = np.mean(x)
        corrMat = corrMat.add(corrMat.transpose())
        tensor[M] = corrMat
        #if s == 10:return tensor;break ## test just 10 windows
        s+=1
        if o.verbose==True:
            sys.stderr.write('\r'+'' *0)
            sys.stderr.write(str(int(s*100/fl))+"%")
            sys.stdout.flush()

    return tensor
def OTU_profile_2(df):
    '''
    It generates OTU abndance profile
    :return:
    '''

    permut = [x for x in it.combinations(df.index.values.tolist(),2)]   # All comparisons tuples (Oi,Oj)
    tensor = {}
    p = 0
    corrInfo = {}
    for i in permut:
        if o.win < 9: o.win = 10
        windowsize = o.win
        s = 0
        corr_array = []
        window = []
        while s+windowsize != len(df.columns)+1:     # run a sliding window
            frames = df.ix[(i[0],i[1]),s:s+windowsize]
            dcos = np.dot(frames.ix[0,:],frames.ix[1,:])/(np.linalg.norm(frames.ix[0,:])*np.linalg.norm(frames.ix[1,:]))
            x = []
            for q in frames.columns:x.append(int(q.split('_')[0]))  # mean age of the window
            M = np.mean(x)
            corr_array.append(dcos)
            window.append(M)
            s += 1
        corrInfo[(i[0],i[1])] = corr_array
        #if p == 10 : break
        p+=1
        if o.verbose==True:
            sys.stderr.write('\r'+'' *0)
            sys.stderr.write(str(int(p*100/len(permut)))+'%')
            sys.stdout.flush()

    for j in range(len(window)):
        corrMat = pd.DataFrame(index=df.index.values.tolist(), columns=df.index.values.tolist()).fillna(0)  # Empty comparison matrix
        for p in permut:
            corrMat.ix[p[0],p[1]] = corrInfo[p][j]
        corrMat = corrMat.add(corrMat.transpose())
        tensor[window[j]] = corrMat
    return tensor

def dfs_log(df):
    df_log = df.copy()
    for i in range(len(df.columns)):
        #x =  map(list,w_norm.iloc[[i]].values)[0]
        # individual normalization + log2 transformation -> in order to mask very high values.
        y = df_log.ix[:,i]
        ylog =  y.apply(lambda x: np.log2(x) if x != 0 else None).fillna(0)
        df_log.ix[:,i] = ylog
    return df_log

def retrieve_matrix(inf,libs,tl,design):
    '''Input: '.cons.tax.summary' file from mothur output and design file. This function retrieve the input matrix
        at some taxonomic level (tl); by default at level of 6 -> gender
    '''

    ## getting the OTU count matrix from mothur - without replicate names ##
    matrix, lib_size = [],[]
    lib_factor,names = {},{}
    colnames = []

    for line in open(inf,'r'):
        line = line.strip("\n").split("\t")
        if line[0] == "taxlevel":
            for i in line: colnames.append(re.sub(' ','',i))
        elif line[0] == str(int(tl)-1): prefix = line[2].strip("\"").split(" ")[0]  # name prefix whether it is replicate
        elif line[0] == tl:
            line[2] = line[2].strip("\"").split(" ")[0]
            line[2] = '_'.join([prefix,line[2]])
            #try:                                                # the following is to handle replicate names
            #    if names[line[2]]:
            #        line[2] = '_'.join([prefix,line[2]])
            #except KeyError:
            #    names[line[2]] = line[2]
            matrix.append(line)
    df = pd.DataFrame(matrix, columns=colnames)
    df.index = list(df.loc[:,'taxon'])
    df.drop(['taxlevel','rankID','taxon','total','daughterlevels',''],axis=1,inplace=True)
    df = df.convert_objects(convert_numeric=True)

    ### Quality control - just leave taxon present in as many as o.taxt #
    # This is why we do not want too many 0 in the matrix, as well as
    # this OTUs (organisms) would be a chymeric, artifacts or just really low abundance - casualities.
    ##

    if o.taxt < 1: o.taxt = int((len(df.columns.values))*o.taxt)
    elif o.taxt >= 1: o.taxt = int((len(df.columns.values))*0.15)
    elif o.taxt > 9: int(o.taxt)
    taxt = o.taxt
    o.taxt = int(len(df.columns.values))-o.taxt

    #df = df.loc[(df.sum(axis=1,numeric_only=True) >= taxt)]
    df = df.loc[(df==0).astype(int).sum(axis=1) < o.taxt]
    df.drop('unclassified_unclassified',axis=0,inplace=True)

    ###* If we like to work with Reads. *###
    ###* changing gender otus table to otus reads by sample *###

    if os.path.exists("otu_reads_by_sample_%s.df"%taxt) == False: df = otu_reads(df,taxt)
    else: df = pickle.load(open("otu_reads_by_sample_%s.df"%taxt, 'r'))

    if o.scaling == True: df = scaling(df)

    sys.stderr.write("-> SO matrix shape: %s \n-> Subject threshold: %s\n"%(df.shape,taxt))

     ## Columns labeling according design #
    global ids
    for line in open(o.design,'r'):
        line = line.strip("\n").split("\t")
        ids[line[0]] = line     # ['5090939SPC', 'female', '81']

    ### Ordering according to age and geting 3 matrices for males, females and all indv. ###

    #df.loc['age'] = ['NA' for n in range(len(list(df.columns.values)))]
    col_names = []
    for j in df.iloc[:]:
        name = re.findall("\d+[BSF]..",j)[0]
        id = ids[name][0]
        age = ids[name][2]
        gender = ids[name][1]
        if gender == 'male': gender = 'M'
        elif gender == 'female': gender = 'F'
        else:gender = 'NA'
        col_names.append("%s_%s_%s"%(age,gender,id))

    df.columns = col_names
    df = df[sorted(df.columns)]

    samples = formating(df) # [all,males,females] -> list of colnames

    df_m = df[samples[1]]
    df_f = df[samples[2]]

    return [df,df_m,df_f]

def otu_reads(ref,taxt):

    global totu; global rotu

    if o.scaling == True: R = o.inRa
    else: R = o.inR

    if  os.path.exists('sample_reads.df') == False: crossRef(o.inT,R);
    elif os.path.exists('tax_otu.dic') == False: crossRef(o.inT,R);
    else:totu = pickle.load(open('tax_otu.dic', 'r'));rotu = pickle.load(open('sample_reads.df', 'r'))

    matrix = []
    idx = []

    for i in ref.index:
        for j in totu[i]:
            id = "%s_%s"%(i,j)
            if j in rotu.index:
                matrix.append(list(rotu.ix[j,:]))
                idx.append(id)
            else:pass

    df = pd.DataFrame(matrix, columns=list(rotu.columns.values),index=idx)
    df = df.convert_objects(convert_numeric=True)
    df = df.loc[(df==0).astype(int).sum(axis=1) < o.taxt]
    #df = df.loc[(df.sum(axis=1,numeric_only=True) >= taxt)]
    pickle.dump(df,open("otu_reads_by_sample_%s.df"%taxt,'w'), pickle.HIGHEST_PROTOCOL)            # Exporting dictionary into binary file to use it afterwards
    return df

def scaling(df):

    ## Normalizing according the library size "scaling" ##

    for i in open(libs,'r'):
        lib_size.append(float((i.strip("\n").split("\t")[1])))
        lib_factor[i.strip("\n").split("\t")[0]] = [1,float((i.strip("\n").split("\t")[1]))]
    for k in lib_factor.keys():
        lib_factor[k][0] = round(lib_factor[k][1]/min(lib_size),4)
        df[k] = df[k].apply(lambda x: round(float(x)/float(lib_factor[k][0]),4))

    return df.ix[:,[2]+range(5,len(colnames))]

def formating(df):

    ### changing names - males and females groups###

    names,males,females = [],[],[]
    for i in range(len(df.columns.values)):
        name = df.columns.values[i]
        id = re.findall("\d+[BSF]..",name)[0]
        try:
            if ids[id][0] == id : names.append(name)
            if ids[id][1] == 'male': males.append(name)
            elif ids[id][1] == 'female': females.append(name)
        except KeyError:
            print name
            pass

    # the following is to check which design samples are missing #
    #k = ids.keys()
    #f = 0
    #for i in names:
    #   k = [x for x in k if x != i]
    #   f+=1
    #print len(k)
    #for i in k: print i
    ###

    return [names,males,females]

def org_corrMatrix(df,numM):
    ''' Input: Matrix i,j {organisms,individuals}
        This function calculates the distance matrix between organisms using
        o.win idividuals to calculate the mean individual-> by default 10.
    '''

    sys.stderr.write ("\n-> Matrix %s calculating...\n"%numM)

    tensor = {}
    if o.win < 9: o.win = 10
    windowsize = o.win

    fl= 1+len(df.columns)-o.win
    ti = 0
    while ti+windowsize != len(df.columns)+1:     # run an overlapping window
        pretensor = distByOrg(org_window(df.ix[:,ti:ti+windowsize]))
        try:
            if tensor[pretensor[1]]:
                name = newname(tensor,pretensor[1])
                tensor[pretensor[1]] = name
        except KeyError:
            tensor[pretensor[1]] = pretensor[0]

        ti+=1
        if o.verbose==True:
            sys.stderr.write('\r'+'' *0)
            sys.stderr.write(str(int(ti*100/fl))+'%')
            sys.stdout.flush()
    return tensor

def newname(tensor,name):
    name = name+0.01
    if tensor[name]:newname(tensor,name)
    else: return name


def org_window(win):
    w_norm = win.copy()
    for i in range(len(win.columns)):
        #x =  map(list,w_norm.iloc[[i]].values)[0]
        # individual normalization + log2 transformation -> in order to mask very high values.
        y = w_norm.ix[:,i]
        ylog =  y.apply(lambda x: np.log2(x) if x != 0 else None).fillna(0)
        w_norm.ix[:,i] = ylog.apply(lambda x: round((x - ylog.min())/(ylog.max() - ylog.min()),4) )
    return w_norm

def distByOrg(df):
    ''' Input: matrix of {organisms,individuals(o.win size)}
        Calculates the distance between organism using a mean individual
        Dist_Oij = |Oi - Oj|
    '''
    #exp = [x for x in it.combinations([1,2,3],2)]
    permut = [x for x in it.combinations(df.index.values.tolist(),2)]   # All comparisons tuples (Oi,Oj)
    corrMat = pd.DataFrame(index=df.index.values.tolist(), columns=df.index.values.tolist()).fillna(0)  # Empty comparison matrix
    x = []  # age list of the window participants
    #print df.ix[:5,:10],"\n***\n"

    for i in df.columns:x.append(int(i.split('_')[0]))  # mean age of the window
    M = np.mean(x)                                      # mean age of the window
    lc=0
    for i in permut:
        org_v0 = []
        org_v1 = []
        df.loc[i[0]].apply(lambda x: org_v0.append(float(x)))
        df.loc[i[1]].apply(lambda x: org_v1.append(float(x)))
        #yij = np.sqrt()
        #org_diff = [org_v0[x] - org_v1[x] for x in range(len(org_v0))]      # Distance calculating
        org_diffabs = [abs(org_v0[x] - org_v1[x]) for x in range(len(org_v0))] # absolut Distance calculating

        #org_assoc_weight = 1 - np.std(org_diff,ddof=1)  # 1 -> == ; 0 -> !=

        ## Test to select the best metric which describes the dynamics of the data
        #** selected DIFF STD  and Scaled 0-1 with: max(diffabs)  **#

        if max(org_diffabs) == 0 : org_assoc_weight = 0;
        else:
            org_assoc_weight = 1-abs(np.mean(org_diffabs))/max(org_diffabs)
        #
        #if max(org_diffabs) == 0 : org_assoc_weightabs = 0;orgcomb = 0
        #else:
        #    org_assoc_weightabs = 1-np.mean(org_diffabs)/max(org_diffabs)
        #    orgcomb = 1-abs(np.mean(org_diff))/max(org_diffabs)
        #    mb = 1-abs(np.std(org_diff,ddof=1))/max(org_diffabs)
        #    mcnab = abs(np.std(org_diff,ddof=1)/(abs(max(org_diff))-abs(min(org_diff))))
        #
        #std_scal_0 = max(org_diffabs)-min(org_diffabs)
        #std_scal_1 = max(org_diff)-min(org_diff)
        #std_scal_2 = math.pow(max(org_diff) -min(org_diff),2)
        #std_scal_x = np.mean(org_diff)
        #if std_scal_x == 0: std_x = 0
        #else: std_x = np.std(org_diffabs,ddof=1)/std_scal_x
        #if std_scal_0 == 0: std_0 = 0
        #else: std_0 = np.std(org_diffabs,ddof=1)/std_scal_0
        #if std_scal_1 == 0: std_1 = 0
        #else: std_1 = np.std(org_diff,ddof=1)/std_scal_1
        #if std_scal_2 == 0: std_2 = 0
        #else:std_2 = math.sqrt(math.pow(np.std(org_diff,ddof=1),2)/std_scal_2)
        #
        ##print org_v0
        ##print org_v1
        #print 'diffs: ', org_diffabs, ', ',org_diff#,"\n"
        ##print 'mean: ',np.mean(org_diffabs),', ',np.mean(org_diff)
        #print 'std: ',np.std(org_diffabs,ddof=1),', ',np.std(org_diff,ddof=1),np.mean(max(org_diff))
        ##print 'max: ',max(org_diffabs),', ',max(org_diff), 'min: ',min(org_diffabs),', ',min(org_diff)
        ##print 'ratio: ',np.mean(org_diffabs)/max(org_diffabs), ', ',np.median(org_diff)/max(org_diff)
        ##print 'scal ', std_scal_0,std_scal_1,std_scal_2
        #print std_0, std_1,std_2,std_x,"\n"
        ##print 'weight: ',org_assoc_weightabs,', ',org_assoc_weight, ', ', orgcomb, ' std: ',ma,', ',mb,', ',mcab, mcnab ,'(',mcnabo,')'
        ##print 'CV: ',cvab,cvn,' Quantile: ',qvab,qvn,np.percentile(org_diffabs,75)
        ##print org_assoc_weightabs,', ',org_assoc_weight,', ',mcab,', ',mcnab,"\n"
        #lc+=1
        #if lc==100:sys.exit('stop')


        corrMat.ix[i[0],i[1]] = org_assoc_weight
    corrMat = corrMat.add(corrMat.transpose())
    return [corrMat,M]


def crossRef(t,r):
    fl = 0
    global totu; global rotu

    T,R = {},[]
    for i in open(t,'r'):
        if fl > 0:
            i = i.strip("\n \;").split("\t")
            t = i[2].split(';')
            gen = re.sub("\(\d+\)",'',t[-1]).strip("\"").split(" ")[0]
            fly = re.sub("\(\d+\)",'',t[-2]).strip("\"").split(" ")[0]
            name = '_'.join([fly,gen])
            try:
                if T[name]:
                    T[name].append(i[0])
            except KeyError: T[name] = [i[0]]
        fl +=1
    fl = 0
    for j in open(r,'r'):
        j = j.split("\t")
        if fl == 0:colnames=j[1:2]+j[3:]
        if fl>0: R.append(j[1:2]+j[3:])
        fl +=1
    dfr = pd.DataFrame(R,columns=colnames)
    dfr.index = dfr.ix[:,0]
    dfr.drop('Group',axis=1,inplace=True)
    dfr = dfr.transpose()
    pickle.dump(T,open("tax_otu.dic",'w'), pickle.HIGHEST_PROTOCOL)            # Exporting ictionary into binary file to use it afterwards
    pickle.dump(dfr,open("sample_reads.df",'w'), pickle.HIGHEST_PROTOCOL)            # Exporting dataframe into binary file to use it afterwards

    totu = T
    rotu = dfr

if __name__ == "__main__": main()