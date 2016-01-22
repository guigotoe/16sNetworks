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
pd.set_option('display.max_columns', 100)
pd.set_option('display.width', 1000)
outf = 'data2.txt'
out = open(outf, 'w')
out.close()
#log_out = open('individuals.txt', 'w')
#log_out.close()     # Re-starting the log file...


## Main ##

def main():
    f1 = codecs.open('../docs/focus_lukas.csv','r','utf-8')
    f2 = codecs.open('../docs/bspspc_lukas.csv','r','utf-8')

    info={}
    c = 0

##  Conditions  ##
# # 0
# all_asthma 1
# BMI 2         => < 25 and bmi > 20
# $$ cancer 3
# chd 4
# chr_diarrh 5
# CRP 6
# diabetes 7
# dvt_pe_rec 8  => Has a deep venous thrombosis or pulmonary embolism been diagnosed during the last 2 months?
# FOC_rec_type 9
# glucose 10    => < 115 and glucose > 55
# $$ heart_att 11
# ibd_new 12
# ibs 13
# liverdisease 14
# $ neuropathy 15
# orgtrans 16
# $ phlebitis 17
# quality_check 18
# $ respi_dis 19
# sample_age 20 => > 60 and age_years < 80
# $ varices 21
# $ veninsuff 22

## codification:
# 0	unknown
# 1	yes
# 2	no
# 7777777	unreadable/inconsistent
# 9999999	missing
# NA   NA


    #for line in f1:
    #    line = line.strip("\n").split("\t")
        #if c == 0:
        #    for index, g in enumerate(line): print "#",g,index
    #    if c != 0:
            #if 60<=line[20]<80:#and  20<= line[2]<25 and 55<=line[10]<115
    #        if line[20]==80:
    #                print line[20]
    #    c+=1

    focus = pd.DataFrame.from_csv('../docs/focus_lukas.csv',header=0, sep='\t').fillna(value='NA')

    bsp = pd.DataFrame.from_csv('../docs/bspspc_lukas.csv',header=0, sep='\t').fillna(value='NA')

    ## *** FOCUS Group 3: Between 60 - 79 *** ##

    focusG3 = focus[(focus['sample_age'].isin(range(60, 80))) & (focus['BMI']>= 20)&(focus['BMI']< 25) & (focus['glucose']>=55)&(focus['glucose']<115)\
    & (focus['all_asthma']==2) & (focus['cancer']==2) & (focus['chd']==2) & (focus['chr_diarrh']==2) & (focus['diabetes']==2) \
    & (focus['heart_att']==2) & (focus['ibd_new']==2) & (focus['ibs']==2) & (focus['liverdisease']==2) & (focus['neuropathy']==2)\
    & (focus['orgtrans']==2) & (focus['phlebitis']==2) & (focus['neuropathy']==2) & (focus['dvt_pe_rec'].isin([2,'NA'])) & (focus['respi_dis']=='2') \
    & (focus['varices']==2) & (focus['veninsuff']==2)]

    out = open(outf, 'a')
    out.write("codification:\nNA = NA\t 0 = unknown\t 1 = yes\t 2 =	no\t 7777777 = unreadable/inconsistent\t 9999999 = missing\n\
Threshold of constraints: 115 > glucose > 55 \t25 > bmi > 20\n")

    #out.write("\nGroup G2[40-60]!\tTotal No. %s;\tMale No. %s;\tFemale No. %s:\n"%(len(G2a),G2a['Gender'].value_counts()['male'],G2a['Gender'].value_counts()['female']))

    out.write("FOCUS Very healthy: %s individuals\n"%len(focusG3))
    focusG3.to_csv(out,sep="\t")
    ## *** FOCUS Group 4: from 80 *** ##

    focusG5t = focus[(focus['sample_age'].isin(range(80, 150))) & (focus['ibd_new']==2) & (focus['ibs']==2) & (focus['chr_diarrh']==2) & (focus['diabetes']==2) ]#& (focus['BMI']>= 20)&(focus['BMI']< 25) & (focus['glucose']>=55)&(focus['glucose']<115)]

    out.write("\n\nFOCUS Elderly > 80 without without chronic_diarrhea, IBD, IBS and Diabetes: %s individuals\n"%len(focusG5t))
    focusG5t.to_csv(out,sep="\t")

     ## *** BSP/SPC Group 3: Between 60 - 79 *** ##

    bspG3 = bsp[(bsp['t37_sample_age_F1'].isin(range(60, 80))) & (bsp['t276_BMI_BL']>= 20)&(bsp['t276_BMI_BL']< 25) & (bsp['t10_glucose_mg/dl_F1']>=55)&(bsp['t10_glucose_mg/dl_F1']<115)\
    & (bsp['t635_allergic_asthma_F1']==2) & (bsp['t903_diabetes_F1']==2) & (bsp['t789_IBD_F1']==2) & (bsp['t790_chronic_diarrhea_F1']==2) \
    & (bsp['t792_IBS_F1']==2) & (bsp['t854_cancer_F1']==2) & (bsp['t907_heart_attack_F1']==2)   & (bsp['t906_liver_disease_F1']==2) & (bsp['t886_neuropathy_F1']==2)\
    & (bsp['t813_organtransplantation_F1']==2) & (bsp['t802_phlebitis_F1']==2) & (bsp['t913_respiratory_disease_F1']==2) \
    & (bsp['t911_varices_F1']==2) & (bsp['t803_venous_insufficiency_F1']==2)]

    out.write("\n\nBSP/SPC Very healthy: %s individuals\n"%len(bspG3))
    bspG3.to_csv(out,sep="\t")
    ## *** BSP/SPC Group 4: from 80 *** ##

    #bspG5t = bsp[(bsp['t37_sample_age_F1'].isin(range(80, 150))) & (bsp['t276_BMI_BL']>= 20)&(bsp['t276_BMI_BL']< 25) & (bsp['t10_glucose_mg/dl_F1']>=55)&(bsp['t10_glucose_mg/dl_F1']<115)]

    #out.write("\n\nBSP/SPC Elderly > 80 and BMI and Glucose as constraints: %s individuals\n"%len(bspG5t))
    #bspG5t.to_csv(out,sep="\t")

    bspG5ta = bsp[(bsp['t37_sample_age_F1'].isin(range(80, 150))) & (bsp['t789_IBD_F1']==2) & (bsp['t790_chronic_diarrhea_F1']==2) & (bsp['t792_IBS_F1']==2) & (bsp['t903_diabetes_F1']==2)]
    print len(bsp[(bsp['t37_sample_age_F1'].isin(range(80, 150))) & (bsp['t789_IBD_F1']==2) & (bsp['t790_chronic_diarrhea_F1']==2) & (bsp['t792_IBS_F1']==2) & (bsp['t903_diabetes_F1']==2)])

    out.write("\n\nBSP/SPC Elderly > 80 without chronic_diarrhea, IBD, IBS and Diabetes: %s individuals\n"%len(bspG5ta))
    bspG5ta.to_csv(out,sep="\t")
    out.close()

    print 'Done!'

## Functions ##


if __name__ == '__main__':main()