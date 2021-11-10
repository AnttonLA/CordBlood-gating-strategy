#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  1 15:17:31 2020

This is the main python file where the full gating strategy is implemented!

@author: antton
"""

import aligater as ag
#import pandas as pd
#import numpy as np
from math import inf
#from styleframe import StyleFrame, utils 
#from random import randint
import datetime
#import re
import os
import sys


#Get file list
#Get repeats, store paths to them in 'repeats_filepaths' variable

path_to_files = "/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs"

def get_blacklist():
    #TODO: the plan was to edit the strings so that only file name is returned
    #but AliGater's mask feature seems to work fine with the full path name
    blacklist = ['/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/181203 CB/C7 98.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190401 CB/D3 333.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190610 CB/C7 529.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190610 CB/C8 530.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190701 CB/C6 616.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190715 CB/D6 651.fcs',
                 '/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs/190909 CB/E8 848.fcs']

    return blacklist

def get_filepaths(maxDate=None, minSampleID=None, maxSampleID=None, includeRepeats=False):
    """
    get_filepaths ceates a list with the exact filepath to each sample fcs file of interest.
    
    :maxDate: is a datetime date object that rejecs all samples OLDER than the specified date
    :maxSampleNum: is an int that limits the max CBID value allowed
    :includeRepeats: is a bool that determines wether or not to include additional
    .fcs files of samples that have been measured more than once.
    
    :return: this function returns a list of filepaths
    
    """
    if maxDate is None:
        maxdate = datetime.date(2021,6,24) #last date we are interested in
        print("No maximum date selected for the samples. Default cutoff is "+str(maxdate))
        
    relevant_folder_paths = []
    all_folders = sorted(os.listdir(path_to_files))  # Take ALL the folder names in 'fcs'
    if '.DS_Store' in all_folders: all_folders.remove('.DS_Store') #damn mac hidden files...
    
    for folder in all_folders:  # For each subfolder, check if we want it and if so append to 'relevabt_folder_paths'
        name_correction = False
        if folder.startswith('Natsumi-'): #Workaround for foldernames that start with 'Natsumi'
            name_correction = True
            og_foldername = folder
            folder = folder[8:]
        folder_date = datetime.date(int('20'+folder[:2]),int(folder[2:4]),int(folder[4:6]))  # Take date from folder name
        if folder_date > maxdate: # If we exceed the max date, stop
            continue
        if name_correction:
            folder = og_foldername
        #Add the folder name to the list of folders of interest
        relevant_folder_paths.append(path_to_files+"/"+folder)  # Add the complete path to the subfolders
    
    fcs_filepaths = []  # List of all fcs files.   
    for subfolder in relevant_folder_paths:
        for file in os.listdir(subfolder):
            if file.endswith('.fcs'):
                fcs_filepaths.append(subfolder+"/"+file)
    
    CBID_dict = dict()
    cbid_list = []
    for fcs in fcs_filepaths:
        cbid_dot_fcs = fcs.split('/')[-1].split(' ')[-1]
        cbid = cbid_dot_fcs.split('.')[0]
        if cbid.isdigit():
            if cbid in CBID_dict.keys(): #Append fcs filepath to the list corresponding to that key
                CBID_dict.setdefault(cbid, []).append(fcs) #append new value to list inside dict
    
            else: #if no key exists for that CBID yet
                CBID_dict[cbid] = [fcs]
    
            cbid_list.append(cbid) #Make list of cbids to use later as keys
    cbid_list = sorted(list(set(cbid_list)), key=int)
    if minSampleID:
        cbid_list = [x for x in cbid_list if int(x) >= int(minSampleID)]
    if maxSampleID:
        cbid_list = [x for x in cbid_list if int(x) <= int(maxSampleID)]
        
    final_filepaths = []
    repeat_filepaths = []
    unique_sample_fcs_filepaths = []
    for cbid in cbid_list:
        if includeRepeats and (len(CBID_dict[cbid]) >= 2): #If YES repeats and cbid has more than 2 samples associated
            final_filepaths.append(CBID_dict[cbid][-2]) #Second to last element of the corresponding list
            final_filepaths.append(CBID_dict[cbid][-1]) #Last element of the corresponding list
            
            #Add repeats only to additional list
            repeat_filepaths.append(CBID_dict[cbid][-2])
        else: #if no repeats or cbid only has single fcs file
            final_filepaths.append(CBID_dict[cbid][-1]) #Last element of the corresponding list
        unique_sample_fcs_filepaths.append(CBID_dict[cbid][-1])
    
    #NOTE: the latest fcs file for each CBID is considered to be the 'original', while older files are considered repeats
        
    print('Available samples:\n\tAll .fcs files: ', len(fcs_filepaths))
    print('\tUnique sample .fcs files included: ', len(unique_sample_fcs_filepaths))
    print('\tNumber of .fcs files in the final list of files to analyze:', len(final_filepaths))
    print('\tNumber of .fcs files belonging to repeats included: ', len(repeat_filepaths))

    return final_filepaths

def out_folder_list(image_path_prefix): # List of all output folders for the images
    
    folder_list = [image_path_prefix+"01-Viable",
        image_path_prefix+"02-Singlet",
        image_path_prefix+"03-CBMCs",
        image_path_prefix+"04-CD45high",        
        image_path_prefix+"05-CD45posCD34pos",
        image_path_prefix+"06-CD3negBACKGATE",
        image_path_prefix+"07-CD4posBACKGATE",
        image_path_prefix+"08-CD19negBACKGATE",       
        image_path_prefix+"09-CD14negBACKGATE",
        
        image_path_prefix+"10-CD16CD56/10A-CD16neg_56negBACKGATE",
        image_path_prefix+"10-CD16CD56/10B-CD16pospos_56posposBACKGATE",
        image_path_prefix+"10-CD16CD56/10C-CD16neg_56posposBACKGATE",
     
        image_path_prefix+"11-LinnegCD34pos",   
        image_path_prefix+"12-CD38neg",    
        image_path_prefix+"13-HSC_MLP_MPP",
        image_path_prefix+"14-B_NK",
        image_path_prefix+"15-CD10pos(MLP)",
        image_path_prefix+"16-CMP_GMP_MEP"]
    
    return folder_list

def gateFullDataset(my_sample, save_images=True, *args, **kwargs):
    """
    gateFullDataset describes the whole gating strategy for a single sample
    
    :my_sample: is the fcs file corresponding to a single blood sample
    :save_images: is a bool that determines wether or not images of the gates will be saved
    """
    #Setup
    ag.agconf.ag_verbose=False ##verbose ON
    sample_filepath = my_sample.filePath
    date_plate = ag.getFileName(ag.getParent(sample_filepath))
    sampleName = ag.getFileName(sample_filepath)
    
    event_cap = 3500000
    if len(my_sample()) > event_cap:  ## Cap samples with too many events. Take only the first 3.5M
        my_sample.fcsDF = my_sample.fcsDF.iloc[0:event_cap]
        sys.stderr.write("Event cap reached. Downsapmled to "+str(event_cap)+" events.\n")
    #Put inconsistent marker names into variables
    markers = my_sample().columns.tolist()
    for marker in markers[6:]:
        if marker.startswith('CD45'):
            CD45 = marker
        if marker.startswith('7AAD') or marker.startswith('CD235'):
            marker_7AAD = marker
    
    if save_images:
        fileName = 'notNone'  # Dirty coding
    else:
        fileName = None
    
    ## Gate 1:  Viable/non-viable separation through 7AAD
    livedead_step1 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=214000, parentGate=None, orientation='vertical', population='lower')
    livedead_step2 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=-800, parentGate=livedead_step1, orientation='horizontal', population='upper')   
    halfcut = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='bilog', thresh=150000, parentGate=livedead_step2, orientation='vertical', population='upper')    
    ylim_back = ag.densityDelimitation(my_sample, xCol= marker_7AAD, parentGate=halfcut, interval=[300,1100], limit_threshold=0.2, direction='right',scale='bilog',T=200)    
    if ylim_back == inf:  # Failsafe in case of no limit
        ylim_back = 900    
    if fileName:
        fileName=image_path_prefix+"01-Viable/"+date_plate+"-"+sampleName+"-Viable.jpeg"    
    livedead_final = ag.gateThreshold(my_sample, name="Viable", xCol='FSC 488/10-A', yCol= marker_7AAD, scale='linear', T=200, yscale='logicle', thresh=ylim_back, parentGate=livedead_step2, orientation='horizontal', population='lower', filePlot=fileName)
    

    ## Gate 2:  Singlets
    singlets_step1 = ag.gateThreshold(my_sample, name="remove_high_clutter", xCol='FSC 488/10-A', yCol='FSC 488/10-H', scale='linear', thresh=210000, parentGate=livedead_final,  orientation='vertical', population='lower')
    if fileName:
        fileName=image_path_prefix+"02-Singlet/"+date_plate+"-"+sampleName+"-Singlet.jpeg"
    singlets = ag.gatePC(my_sample, name="Singlets", xCol='FSC 488/10-A', yCol='FSC 488/10-H', center='centroid', adjustAngle=3,widthScale=2.5,heightScale=3.5, parentGate=singlets_step1, filePlot=fileName)


    ## Gate 3: CBMC cells  
    
    halfcut_middle = ag.gateThreshold(my_sample, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=120000, parentGate=singlets, orientation='vertical', population='upper')

    ylim_bot = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_middle, interval=[10000,20000], limit_threshold=0.05, direction='left',scale='linear')
    if ylim_bot == inf:
        ylim_bot = 10000
    CBMC_step1 = ag.gateCorner(my_sample, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=85000, yThresh=ylim_bot, xOrientation='lower', yOrientation='lower', Outer=True, parentGate=singlets)
    
    halfcut_tail = ag.gateThreshold(my_sample, name="right_tail",  xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', thresh=180000, parentGate=singlets, orientation='vertical', population='upper')
    ylim_top = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000], limit_threshold=0.2, direction='right',scale='linear') + 25000
    if ylim_top == inf:
        ylim_top = 140000

    CBMC = ag.horizontalPath(my_sample, name="CBMC",
                        xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                        startY=ylim_bot+30000, endY=ylim_top, xboundaries=[55000,140000],
                        yboundaries=[ylim_bot+25000,ylim_top+5000],
                        leftRight=True , direction='both',
                        maxStep=2, phi=0.1, bins=100, sigma=1,
                        scale='linear', parentGate=CBMC_step1)
    
    if fileName:
        fileName=image_path_prefix+"03-CBMCs/"+date_plate+"-"+sampleName+"-CBMCs.jpeg"
    
    ag.backGate(my_sample, population=CBMC, background_population=singlets, xCol='FSC 488/10-A', yCol='SSC 488/10-A', scale='linear', T=200, markersize=0.1,  filePlot=fileName)
        
    my_sample.update(ag.AGgate(CBMC, None,'FSC 488/10-A', 'SSC 488/10-A', "CBMC"), QC=True, scale='linear',xlim=[25000,180000], ylim=[0,200000])

    
    ## Gate 4: CD45+ AND CD34+ out of CBMC cells, "Declutter gate" to get rid of CD45- debree 
    xlim_arbitrary = 11000
    right_half = ag.gateThreshold(my_sample, name="right_half", xCol='CD34 PE-Cy7-A' , yCol=CD45, scale='bilog', T=200, thresh=xlim_arbitrary, parentGate=CBMC, orientation='vertical', population='upper')    
        # Get highest density point of CD34 blob to use as reference
    mean, median, sigma, maxVal = ag.axisStats(my_sample(), xCol=CD45, vI=right_half())
    upper_half = ag.gateThreshold(my_sample, name="upper_half", xCol='CD34 PE-Cy7-A', yCol=CD45, scale='bilog', T=200, thresh=maxVal, parentGate=CBMC, orientation='horizontal', population='upper')    
    lowr_halfof_upper_half = ag.gateThreshold(my_sample, name="lower_half_of_upper_half", xCol='CD34 PE-Cy7-A', yCol=CD45, scale='bilog', T=200, thresh=maxVal+3000, parentGate=upper_half, orientation='horizontal', population='lower')   
    xlim_middle = ag.valleySeek(my_sample, xCol='CD34 PE-Cy7-A', interval=[800, 11000], require_local_min=True, scale='bilog', T=200, parentGate=lowr_halfof_upper_half) -1000        
    if fileName:
        fileName=image_path_prefix+"04-CD45high/"+date_plate+"-"+sampleName+"-CD45high.jpeg"
    CD45_high = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal-2900), theta=-40, scale='bilog', T=200, population='upper', parentGate=CBMC, filePlot=fileName)
    
    
    ## Gate 5: CD34+ ratio relative to CD45+
    CD34pos = ag.gateTiltedLine(my_sample, name="CD34+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='lower', parentGate=CD45_high)
    if fileName:
        fileName=image_path_prefix+"05-CD45posCD34pos/"+date_plate+"-"+sampleName+"-CD45posCD34pos.jpeg"
    CD45pos = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='upper', parentGate=CD45_high, filePlot=fileName)

    my_sample.update(ag.AGgate(CD34pos, CD45pos,'CD34 PE-Cy7-A', CD45, "CD34+"), QC=True, scale='bilog', xscale='bilog', yscale='bilog', T=200, xlim=[-1000, 100000], ylim=[1000, 150000])
    ## Gate 6: CD3 out of CD45+
        ## Gate 6A: CD3-        
    ylim_middle = ag.valleySeek(my_sample, xCol='CD3 APC-H7-A', interval=[-1000, 7000], require_local_min=True, scale='bilog', T=200, parentGate=CD45_high)    
    CD3neg = ag.gateThreshold(my_sample, name="CD3-", xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD45_high, orientation='horizontal', population='lower')    
    if fileName:
        fileName=image_path_prefix+"06-CD3negBACKGATE/"+date_plate+"-"+sampleName+"-CD3negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD3neg, background_population=CBMC, xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)  
        ## Gate 6B: CD3+ 
    CD3pos = ag.gateThreshold(my_sample, name="CD3+", xCol=CD45 , yCol='CD3 APC-H7-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD45_high, orientation='horizontal', population='upper')
    
    my_sample.update(ag.AGgate(CD3pos, CD45pos, CD45, 'CD3 APC-H7-A', "CD3+"), QC=False)
    
    ## Gate 7: CD4, CD8 out of CD3+
        ## Gate 7A: CD4+
    xlim_middle = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True, scale='bilog', T=200, parentGate=CD3pos)
    cd4cd8_step1 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=xlim_middle, parentGate=CD3pos, orientation='vertical', population='upper')    
    CD4pos_final = ag.gatePC(my_sample, name="CD4+", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', center='centroid', adjustAngle=3,widthScale=2.5, scale='bilog', T=200, heightScale=3.5, parentGate=cd4cd8_step1)    
    if fileName:
        fileName=image_path_prefix+"07-CD4posBACKGATE/"+date_plate+"-"+sampleName+"-CD4posBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD4pos_final, background_population=CD3pos, xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName) 
   
    my_sample.update(ag.AGgate(CD4pos_final, CD3pos,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD4+_1"), QC=False)
        ##Gate 7B: CD8+
    ylim_middle = ag.valleySeek(my_sample, xCol='CD8 PerCP-Cy5.5-A', interval=[200, 2000], require_local_min=True, parentGate=CD3pos)
    cd4cd8_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=ylim_middle, parentGate=CD3pos, orientation='horizontal', population='upper')    
    xlim_middle_upper = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True, scale='bilog', T=200, parentGate=cd4cd8_step2)   
    cd4cd8_step3 = ag.gateThreshold(my_sample, name="separate_middle_x_axis", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=xlim_middle_upper, parentGate=cd4cd8_step2, orientation='vertical', population='lower')   
    CD8pos_final = ag.gatePC(my_sample, name="CD8+", xCol='CD4 (BV) 510-A' , yCol='CD8 PerCP-Cy5.5-A', center='centroid', adjustAngle=0,widthScale=2, scale='bilog', T=200, heightScale=3.5, parentGate=cd4cd8_step3)
    
    my_sample.update(ag.AGgate(CD8pos_final, CD3pos,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD8+"), QC=False)
    my_sample.update(ag.AGgate(CD4pos_final, CD8pos_final,'CD4 (BV) 510-A','CD8 PerCP-Cy5.5-A',"CD4+_2"), QC=False)

    ## Gate 8: CD19 out of CD3-
        ## Gate 8A: CD19-
    ylim_middle_cd19 = ag.valleySeek(my_sample, xCol='CD19 PE-Texas Red-A', interval=[200, 10000], require_local_min=True, scale='bilog', T=200, parentGate=CD3neg)
    CD19neg = ag.gateThreshold(my_sample, name="CD19-", xCol='CD3 APC-H7-A' , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horizontal', population='lower')
    if fileName:
        fileName=image_path_prefix+"08-CD19negBACKGATE/"+date_plate+"-"+sampleName+"-CD19negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD19neg, background_population=CD3neg, xCol='CD3 APC-H7-A'  , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName) 
        ## Gate 8B: CD19+
    CD19pos = ag.gateThreshold(my_sample, name="CD19+", xCol='CD3 APC-H7-A' , yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horizontal', population='upper', filePlot=fileName)
            
    my_sample.update(ag.AGgate(CD19pos, CD45pos, 'CD3 APC-H7-A', 'CD19 PE-Texas Red-A', "CD19+"), QC=False)

    ## Gate 9: CD14 out of CD3-/CD19-    
        ## Gate 7A: CD14- out of CD19-
    ylim_middle_cd14 = ag.valleySeek(my_sample, xCol='CD14 (BV) 605-A', interval=[0, 3000], require_local_min=True, parentGate=CD19neg, scale='bilog')
    CD14neg = ag.gateThreshold(my_sample, name="CD14-", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A', scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg, orientation='horizontal', population='lower')
    if fileName:
        fileName=image_path_prefix+"09-CD14negBACKGATE/"+date_plate+"-"+sampleName+"-CD14negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD14neg, background_population=CD19neg, xCol='CD16 (BV) 786-A' , yCol='CD14 (BV) 605-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)     
        ## Gate 7B: CD14+ out of CD19-        
    CD14pos = ag.gateThreshold(my_sample, name="CD14+", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A', scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg, orientation='horizontal', population='upper')
        
    my_sample.update(ag.AGgate(CD14pos, CD45pos, 'CD19 PE-Texas Red-A', 'CD14 (BV) 605-A', "CD14+"), QC=False)

    ## Gate 10: CD16 and CD56 values out of CD14-

    aprox_xlim_middle = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True, scale='bilog', T=200, parentGate=CD14neg)    
    left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle, parentGate=CD14neg, orientation='vertical', population='lower')    
    aprox_ylim_middle = ag.valleySeek(my_sample, xCol='CD56 (BV) 650-A', interval=[700, 1300], require_local_min=True, scale='bilog', T=200, parentGate=CD14neg)    
    lower_left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_ylim_middle, parentGate=left_half, orientation='horizontal', population='lower')    
    main_ylim_middle = ag.densityDelimitation(my_sample, xCol='CD56 (BV) 650-A', parentGate=lower_left_half, interval=[0,2000], limit_threshold=0.07, direction='right',scale='linear')    
    if main_ylim_middle == inf:
        main_ylim_middle = aprox_ylim_middle    
    lower_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg, orientation='horizontal', population='lower')    
    upper_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg, orientation='horizontal', population='upper')    
    aprox_xlim_middle2 = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True, scale='bilog', T=200, parentGate=lower_half)    
    CD16neg_56neg_step1 = ag.gateThreshold(my_sample, name="lower_left_side", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2, parentGate=lower_half, orientation='vertical', population='lower')    
    main_xlim_middle = ag.densityDelimitation(my_sample, xCol='CD16 (BV) 786-A', parentGate=CD16neg_56neg_step1, interval=[0,2000], limit_threshold=0.07, direction='right',scale='linear')
    if main_xlim_middle == inf:
        main_xlim_middle = aprox_xlim_middle    
    CD16neg_56neg_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle, parentGate=CD16neg_56neg_step1, orientation='vertical', population='lower')
        ## CD16-CD56-, a.k.a. LIN NEGATIVE
    CD16neg_56neg_final = ag.gatePC(my_sample, name="CD16-CD56-", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', center='centroid', adjustAngle=0,widthScale=2, scale='bilog', T=200, heightScale=3.5, parentGate=CD16neg_56neg_step2)
    if fileName:
        fileName=image_path_prefix+"10-CD16CD56/10A-CD16neg_56negBACKGATE/"+date_plate+"-"+sampleName+"-CD16neg_56negBACKGATE.jpeg" 
    ag.backGate(my_sample, population=CD16neg_56neg_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
    central_section = ag.horizontalPath(my_sample, name="hor_path", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', population='lower', startY=3000, endY=12000,  xboundaries=[0,500000], yboundaries=[2000,15000], leftRight=True , direction='both', maxStep=2, phi=0.1, bins=100, sigma=1, scale='bilog', T=200, parentGate=upper_half)   
         ## CD16++CD56++
    CD16pospos_56pospos_final = ag.gateThreshold(my_sample, name="CD16+CD56+", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle, parentGate=central_section, orientation='vertical', population='upper')
    if fileName:
        fileName=image_path_prefix+"10-CD16CD56/10B-CD16pospos_56posposBACKGATE/"+date_plate+"-"+sampleName+"-CD16pospos_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16pospos_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
    
    my_sample.update(ag.AGgate(CD16pospos_56pospos_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16++CD56++"), QC=False)
        ## CD16-CD56++    
    CD16neg_56pospos_final = ag.horizontalPath(my_sample, name="CD16-CD56++", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', population='upper', startY=3000, endY=12000, xboundaries=[0,500000], yboundaries=[2000,15000], leftRight=True , direction='both', maxStep=2, phi=0.1, bins=100, sigma=1, scale='bilog', T=200, parentGate=upper_half)
    if fileName:
        fileName=image_path_prefix+"10-CD16CD56/10C-CD16neg_56posposBACKGATE/"+date_plate+"-"+sampleName+"-CD16neg_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16neg_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=fileName)
        
    my_sample.update(ag.AGgate(CD16neg_56pospos_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16-CD56++"), QC=False)
        ## CD16++CD56-
    CD16pospos_56neg_final = ag.gateThreshold(my_sample, name="CD16++CD56-", xCol='CD16 (BV) 786-A' , yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2, parentGate=lower_half, orientation='vertical', population='upper')
 
    my_sample.update(ag.AGgate(CD16pospos_56neg_final, CD45pos,'CD16 (BV) 786-A' ,'CD56 (BV) 650-A',"CD16++CD56-"), QC=False)
    
    ## Gate 11: Lin-CD34+ out of CD16-CD56-/CD14-/CD19-/CD3-/CD45+
    
    if fileName:
        fileName=image_path_prefix+"11-LinnegCD34pos/"+date_plate+"-"+sampleName+"-LinnegCD34pos.jpeg"
    CD34pos_step1 = ag.gateTiltedLine(my_sample, name="tilted_gate_cd34", xCol='CD34 PE-Cy7-A' , yCol=CD45, startPoint=(xlim_middle,maxVal), endLimits=(None, maxVal+7500), theta=40, scale='bilog', T=200, population='lower', parentGate=CD16neg_56neg_final, filePlot=fileName)
    linnegCD34pos = ag.gatePC(my_sample, name="Lin-CD34+", xCol='CD34 PE-Cy7-A' , yCol=CD45, center='centroid', adjustAngle=3,widthScale=2.5, heightScale=3.5, parentGate=CD34pos_step1)
    
    my_sample.update(ag.AGgate(linnegCD34pos, CD45pos,'CD34 PE-Cy7-A' , CD45,"linneg_cd34pos"), QC=True,  scale='bilog', xscale='bilog', yscale='bilog', T=200, xlim=[-1000, 80000], ylim=[0, 150000])

    ## Gate 12: CD38 out of Lin-CD34+ (PLACEHOLDER until new CD38 gate is developed)
    
    gated_df = my_sample.fcsDF.loc[CD34pos()].copy()
    cd38_values_list = sorted(list(gated_df['CD38 (BV) 421-A']))
    if cd38_values_list:  #If the list is not empty
        index = int(len(cd38_values_list)*0.3)
        ylim_middle = cd38_values_list[index]
    else:
        return my_sample

        ## Gate 11A: CD38-
    if fileName:
        fileName=image_path_prefix+"12-CD38neg/"+date_plate+"-"+sampleName+"-CD38neg.jpeg"
    CD38neg= ag.gateThreshold(my_sample, name="CD38-", xCol='CD34 PE-Cy7-A' , yCol='CD38 (BV) 421-A',  scale='bilog', T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horizontal', population='lower', filePlot=fileName)
        ## Gate 11B: CD38+
    CD38pos = ag.gateThreshold(my_sample, name="CD38+", xCol='CD34 PE-Cy7-A' , yCol='CD38 (BV) 421-A',  scale='bilog', T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horizontal', population='upper')
    
    ## Gate 13: HSCs, MPPs and MLPs out of CD38-   
    ylim = ag.densityDelimitation(my_sample, xCol='CD90 PE (R-phycoerythrin)-A', interval=[200,5000], limit_threshold=0.5, direction='right',scale='bilog',T=200, parentGate=CD38neg)    
    xlim = ag.densityDelimitation(my_sample, xCol='CD45RA FITC-A', interval=[-100,300], limit_threshold=0.05, direction='right',scale='bilog',T=200, parentGate=CD38neg)   
    if ylim == inf:  # Plan B... RATHER PROBLEMATIC BECAUSE IT'S INCONSISTENT!!!
        right_hand= ag.gateThreshold(my_sample, name="remove_clutter_3", xCol='CD45RA FITC-A', yCol='CD90 PE (R-phycoerythrin)-A', scale='bilog', T=200, thresh=xlim, parentGate=CD38neg, orientation='vertical', population='upper')        
        ylim = ag.densityDelimitation(my_sample, xCol='CD90 PE (R-phycoerythrin)-A', interval=[200,5000], limit_threshold=0.2, direction='right',scale='bilog',T=200, parentGate=right_hand)
        if ylim == inf:
            ylim = 5000  # Plan C...
    if fileName:
        fileName=image_path_prefix+"13-HSC_MLP_MPP/"+date_plate+"-"+sampleName+"-HSC_MLP_MPP.jpeg"
    HSC, irrelevant, MLP, MPP = ag.quadGate(my_sample, names=['HSC', 'NA', 'MLP', 'MPP'], xCol='CD45RA FITC-A', yCol='CD90 PE (R-phycoerythrin)-A', xThresh=xlim, yThresh=ylim, scale='bilog', T=200, parentGate=CD38neg, filePlot=fileName)

    my_sample.update(ag.AGgate(HSC, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"HSC_1"), QC=True, scale='bilog', xscale='bilog', yscale='bilog', T=200, xlim=[0, 1000], ylim=[300, 11000])
    my_sample.update(ag.AGgate(HSC, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"HSC_2"), QC=False)
    my_sample.update(ag.AGgate(MLP, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MLP_1"), QC=False)
    my_sample.update(ag.AGgate(MLP, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MLP_2"), QC=False)
    my_sample.update(ag.AGgate(MPP, CD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MPP_1"), QC=False)
    my_sample.update(ag.AGgate(MPP, linnegCD34pos,'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',"MPP_2"), QC=False)
    
    ## Gate 14: non-BNK (CD10neg) and B-NK out of CD38+
    ylim_cd10 = ag.densityDelimitation(my_sample, xCol='CD10 APC (Allophycocyanin)-A', parentGate=halfcut, interval=[-100,800], limit_threshold=0.2, direction='right',scale='bilog',T=200)
    
        ## B-NK progenitors, CD10+?
    if fileName:
        fileName=image_path_prefix+"14-B_NK/"+date_plate+"-"+sampleName+"-B_NK.jpeg"        
    BNK_prog = ag.gateCorner(my_sample, name="BNK", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', xThresh=200, yThresh=ylim_cd10, xOrientation='upper', yOrientation='upper', Outer=False, scale='bilog', T=200, parentGate=CD38pos, filePlot=fileName)
    
    my_sample.update(ag.AGgate(BNK_prog, CD38pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_1"), QC=False) 
    my_sample.update(ag.AGgate(BNK_prog, CD34pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_2"), QC=False)
    my_sample.update(ag.AGgate(BNK_prog, linnegCD34pos,'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A',"B-NK_3"), QC=False)

        ## CD10-
    nonBNK= ag.gateThreshold(my_sample, name="CD10-", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=CD38pos, orientation='horizontal', population='lower')

    ## Gate 15: CD10 out of MLP
    if fileName:
        fileName=image_path_prefix+"15-CD10pos(MLP)/"+date_plate+"-"+sampleName+"-CD10pos(MLP).jpeg"     
    CD10pos_MLP = ag.gateThreshold(my_sample, name="cd10+_MLP", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=MLP, orientation='horizontal', population='upper',filePlot=fileName)
    CD10neg_MLP = ag.gateThreshold(my_sample, name="cd10-_MLP", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=ylim_cd10, parentGate=MLP, orientation='horizontal', population='lower')

    my_sample.update(ag.AGgate(CD10pos_MLP, MLP,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_1"), QC=False)
    my_sample.update(ag.AGgate(CD10pos_MLP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_2"), QC=False)
    my_sample.update(ag.AGgate(CD10pos_MLP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10+(MLP)_3"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, MLP,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_1"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_2"), QC=False)
    my_sample.update(ag.AGgate(CD10neg_MLP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CD10-(MLP)_3"), QC=False)    
    
    ## Gate 16: CNPs, GMPs and MEPs out of non-BNK CD38+ cells
    ylim = ag.densityDelimitation(my_sample, xCol='CD135 (BV) 711-A', interval=[-100,800], limit_threshold=0.5, direction='left',scale='bilog',T=200, parentGate=nonBNK)
    
    xlim_cut_tail = ag.densityDelimitation(my_sample, xCol='CD45RA FITC-A', interval=[-100,800], limit_threshold=0.2, direction='right',scale='bilog',T=200, parentGate=nonBNK)
    if fileName:
        fileName=image_path_prefix+"16-CMP_GMP_MEP/"+date_plate+"-"+sampleName+"-CMP_GMP_MEP.jpeg"
    CMP, GMP, irrelevant, MEP = ag.quadGate(my_sample, names=['CMP', 'GMP', 'NA', 'MEP'], xCol='CD45RA FITC-A', yCol='CD135 (BV) 711-A', xThresh=xlim_cut_tail, yThresh=ylim, scale='bilog', T=200,  parentGate=nonBNK, filePlot=fileName)
    
    my_sample.update(ag.AGgate(CMP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CMP_1"), QC=False)
    my_sample.update(ag.AGgate(CMP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"CMP_2"), QC=False)
    my_sample.update(ag.AGgate(MEP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"MEP_1"), QC=False)
    my_sample.update(ag.AGgate(MEP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"MEP_2"), QC=False)
    my_sample.update(ag.AGgate(GMP, CD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"GMP_1"), QC=False)
    my_sample.update(ag.AGgate(GMP, linnegCD34pos,'CD45RA FITC-A', 'CD135 (BV) 711-A',"GMP_2"), QC=False)
    
    return my_sample


if __name__ == '__main__':
    
    ## Folder structure to save images at. Maybe should be a function
        #Note that 'image_path_prefix' is a global variable and 
        #is used by 'gateFullDataset'
    date_string=str(datetime.date.today()) #2021-10-25 , for example
    week_string="w"+str(datetime.date.today().isocalendar()[1]) #w43 , for example
    week_n_year_string=week_string+"_"+str(datetime.date.today().isocalendar()[0]) #w43-2021 , for example
    image_path_prefix="/home/antton/Projects/CB_Data_Analysis/output/gating/aligater_images/images_"+week_n_year_string+"/"
    
    if not os.path.exists(image_path_prefix):    
        os.mkdir(image_path_prefix) #Create parent folder for all images at this date
    out_folders = out_folder_list(image_path_prefix) #List containing desired folder structure to save images in
    if not os.path.exists(image_path_prefix+"10-CD16CD56/"): #create subfolder so sub-sub folders can be created next step
            os.mkdir(image_path_prefix+"10-CD16CD56/")
    for folder in out_folders: #Actual creation of the folders
        if not os.path.exists(folder):
            os.mkdir(folder)
    
    batch_size = 500
    lowest_sampleID=3002
    highest_sampleID=3002#batch_size
    
    #fcs files to gate        
     
    while True:#lowest_sampleID < 3400:
        print("Gating all samples from CBID"+str(lowest_sampleID)+" to CBID"+str(highest_sampleID))
        filepaths = get_filepaths(minSampleID=lowest_sampleID,maxSampleID=highest_sampleID) #Get list of paths to each fcs file of interest

        #Define experiment object
        CB_exp=ag.AGExperiment(filepaths, filters=['fcs'], mask=['30min','45min','Neg','test']+get_blacklist(),\
                               experiment_name="cord_blood_experiment-"+week_n_year_string+"-samples_"+str(lowest_sampleID)+"_to_"+str(highest_sampleID),\
                               flourochrome_area_filter=True, QC=True, QCbins=128)
        CB_exp.apply(gateFullDataset) #Apply this function to every sample
        
        CB_exp.printExperiment("/home/antton/Projects/CB_Data_Analysis/output/gating/aligater_output/cblood_phenotypes-"+week_n_year_string+"-samples_"\
                               +str(lowest_sampleID)+"_to_"+str(highest_sampleID)+".txt")
        break
        lowest_sampleID=highest_sampleID+1 #Update limits to measure next batch
        highest_sampleID=highest_sampleID+batch_size
