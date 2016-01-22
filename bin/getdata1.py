#!usr/bin/python

######################################################
## By Guillermo Torres PhD. - IKMB                 ##
####################################################

## Libraries importing

import pandas as pd
import numpy as np
import re,codecs

## Parsing Options ##


## Global variables ##
log_out = open('groups.txt', 'w')
log_out.close()     # Re-starting the log file...


## Main ##

def main():
    f1 = codecs.open('../Focus16s.txt','r','utf-8')
    f2 = codecs.open('../BSPSPC16s.txt','r','utf-8')
    f3 = codecs.open('../Metagenomes.txt','r','utf-8')

    info={}
    c = 0

    for line in f1:
        if c != 0:
            line=line.strip("\n\r").split("\t")
            if re.findall("FOC$",line[0]): sample = 'FOC'
            elif re.findall("SPC$",line[0]): sample = 'SPC'
            elif re.findall("BSP$",line[0]): sample = 'BSP'
            info[line[0]]=[sample,line[1],line[3],line[4],0,line[-1],line[-2]] # [sample(Focus/BSP/SPC),rec_type,popgen,16S,microbiome,gender,age]
        c+=1
    c = 0
    for line in f2:
        if c != 0:
            line=line.strip("\n\r").split("\t")
            if re.findall("FOC$",line[0]): sample = 'FOC'
            elif re.findall("SPC$",line[0]): sample = 'SPC'
            elif re.findall("BSP$",line[0]): sample = 'BSP'
            info[line[0]]=[sample,'NA',1,1,0,line[-1],line[-2]] # [sample(Focus/BSP/SPC),rec_type,popgen,16S,microbiome,gender,age]
        c+=1
    c = 0
    for line in f3:
        if c != 0:
            line=line.strip("\n\r").split("\t");
            try: info[line[0]][4] = 1
            except KeyError: print "not in the list %s"%line[0]
        c+=1

    info = pd.DataFrame.from_dict(info,orient="index")
    info = info.rename(columns = {0:'Sample',1:'Rec_Type',2:'Popgen',3:'16S',4:'Metagenome',5:'Gender',6:'Age_years'})
    info = info.replace(to_replace='NA', value=np.nan)
    info['status'] = pd.cut(info["Age_years"].fillna(0).astype(int), [-1,10,40,59,79,120],labels=['NaN','G1[40-]','G2[40-60]','G3[60-80]','G4[80+]'])   # with numpy instead the array-> np.arange(20, 100+20, 20)
    info['Rec_Type'] = info['Rec_Type'].fillna('NaN')
    info['Popgen'] = info['Popgen'].fillna('NaN')
    info['16S'] = info['16S'].fillna('NaN')
    #y = pd.DataFrame(columns=list(info.columns))
    for x in info.groupby(by='Rec_Type'):
        if x[0] == 'clinic':info_cli = x[1]
        if x[0] == 'NaN':info_nan = x[1]
        if x[0] == 'registry':info_reg = x[1]
    info_reg = info_reg.append(info_nan,ignore_index=False)
    for x in info_reg.groupby(by='status'):
        if x[0] == 'NaN':Nan = x[1]
        elif x[0] == 'G1[40-]':G1 = x[1]
        elif x[0] == 'G2[40-60]':G2 = x[1]
        elif x[0] == 'G3[60-80]':G3 = x[1]
        elif x[0] == 'G4[80+]':G4 = x[1]
    G1a,G1b,G2a,G2b,G3a,G3b,G4a,G4b = [pd.DataFrame(columns=list(info.columns))]*8
    for x in G1.groupby(by='Metagenome'):
        if x[0] == 1: G1a = x[1]
        if x[0] == 0: G1b = x[1]
    for x in G2.groupby(by='Metagenome'):
        if x[0] == 1: G2a = x[1]
        if x[0] == 0: G2b = x[1]
    for x in G3.groupby(by='Metagenome'):
        if x[0] == 1: G3a = x[1]
        if x[0] == 0: G3b = x[1]

    log_out = open('groups.txt', 'a')
    log_out.write("Group G1[40-]!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G1a),G1a['Gender'].value_counts()['male'],G1a['Gender'].value_counts()['female']))
    G1a.sort(columns=['Sample','Age_years']).to_csv(log_out,sep="\t",index_label='ID')
    log_out.write("\nGroup G2[40-60]!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G2a),G2a['Gender'].value_counts()['male'],G2a['Gender'].value_counts()['female']))
    G2a.sort(columns=['Sample','Age_years']).to_csv(log_out,sep="\t",index_label='ID')
    log_out.write("\nGroup G3[60-80]!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G3a),G3a['Gender'].value_counts()['male'],G3a['Gender'].value_counts()['female']))
    G3a.sort(columns=['Sample','Age_years']).to_csv(log_out,sep="\t",index_label='ID')
    log_out.write("\nIf We need more metagenomes for this group G3[60-80], there are enough with just 16S ...\
\nGroup G3[60-80]! Without Metagenome!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G3b),G3b['Gender'].value_counts()['male'],G3b['Gender'].value_counts()['female']))
    #G3b.sort(columns=['Sample','Age_years']).to_csv(log_out,sep="\t",index_label='ID')
    log_out.write("\nGroup G4[80+]! Without Metagenome!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G4),G4['Gender'].value_counts()['male'],G4['Gender'].value_counts()['female']))
    G4.sort(columns=['Sample','Age_years']).to_csv(log_out,sep="\t",index_label='ID')



    print 'Done!'

## Functions ##


if __name__ == '__main__':main()