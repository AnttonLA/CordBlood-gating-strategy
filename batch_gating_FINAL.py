#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script is used to run the AliGater gating strategy on all the samples in the dataset.
It describes and implements the gating strategy and either calculates or reads the parameters used for each gate.

AliGater takes a list of .fcs file names as input to create an AGExperiment object. You can then call the .apply()
method with a function describing a gating strategy, and that strategy will be run for any number of .fcs files stored
by the object. In this case, the function 'gate_full_dataset' is used to describe the gating strategy.

The gating strategy contains 18 steps, which are divided into pre-CD38 and post-CD38 steps. The pre-CD38 steps are run
only once per each sample, while the post-CD38 steps are run four times per sample, changing the CD38 gate each time.




@author: Antton Lamarka
"""

import aligater as ag
from math import inf
import datetime
import os
import sys

# Get file list
# Get repeats, store paths to them in 'repeats_filepaths' variable
import numpy as np

path_to_files = "/home/antton/cbio3/data/BloodVariome/Cord Blood/fcs"


def get_blacklist():
    """
    :return: get_blacklist returns a list of file names that are to  be used with the "mask"
    feature of AliGater, so that these files are ignored. This helps avoid crashes.
    The file names are either written in the script or read from files.
    """
    blacklist = []

    # Add list of samples that make gating strategy crash
    with open('../data/interim/variety_text_file_lists/gating_blacklist.txt', 'r') as blacklist_file:
        for line in blacklist_file:
            blacklist.append(line.rstrip())

    # Add samples with incorrect markers/without any markers to avoid crashes
    # TODO: edit this file
    with open("../data/interim/variety_text_file_lists/missing_markers_Aurelie_SOLVED/samples_without_markers.txt",
              'r') as infile:
        for line in infile:
            blacklist.append(line.rstrip())

    # Add faulty samples (confirmed as faulty by Aitzkoa) to the blacklist
    with open("../output/AliGater_gating_QC/w49_2021/confirmed_faulty_samples_w49.txt", 'r') as infile:
        for line in infile:
            blacklist.append('/home/antton/cbio3/data/BloodVariome/Cord Blood/' + line.rstrip())

    return blacklist


def get_filepaths(max_date=None, min_sample_id=None, max_sample_id=None, include_repeats=False):
    """
    get_filepaths creates a list with the global filepath to each fcs file corresponding to a sample of interest.
    The files on the list are then intended to be used in Flow Cytometry gating using AliGater.

    :param max_date: a datetime date object that rejects all samples OLDER than the specified date
    :param min_sample_id: an int or number str that limits the min CBID value allowed
    :param max_sample_id: an int or number str that limits the max CBID value allowed
    :param include_repeats: a bool that determines whether to include additional .fcs files of samples that have been
    measured more than once.
    
    :return: this function returns a list of filepaths
    
    """
    if max_date is None:
        max_date = datetime.date(2022, 6, 24)  # last date we are interested in
        print("No maximum date selected for the samples. Default cutoff is " + str(max_date))

    relevant_folder_paths = []
    all_folders = sorted(os.listdir(path_to_files))  # Take ALL the folder names in 'fcs'
    for folder_to_ignore in ['.DS_Store', 'Aurelie-to-Antton']:  # list of folders not to take from
        if folder_to_ignore in all_folders:
            all_folders.remove(folder_to_ignore)

    for folder in all_folders:  # For each sub-folder, check if we want it and if so append to 'relevant_folder_paths'
        name_correction_needed = False
        og_folder_name = folder
        if og_folder_name.startswith('Natsumi-'):  # Workaround for folder-names that start with 'Natsumi'
            name_correction_needed = True
            folder = folder[8:]  # Remove 'Natsumi-' from the folder name
        # Get the date from folder name:
        folder_date = datetime.date(int('20' + folder[:2]), int(folder[2:4]), int(folder[4:6]))
        if folder_date > max_date:  # If we exceed the max date, stop
            continue
        if name_correction_needed:
            folder = og_folder_name
        # Add the folder name to the list of folders of interest
        relevant_folder_paths.append(path_to_files + "/" + folder)  # Add the complete path to the folders

    fcs_filepaths = []  # List of all fcs files.   
    for subfolder in relevant_folder_paths:
        for file in os.listdir(subfolder):
            if file.endswith('.fcs'):
                fcs_filepaths.append(subfolder + "/" + file)

    cbid_dict = dict()
    cbid_list = []
    for fcs in fcs_filepaths:
        cbid_dot_fcs = fcs.split('/')[-1].split(' ')[-1]
        cbid = cbid_dot_fcs.split('.')[0]
        if cbid.isdigit():
            if cbid in cbid_dict.keys():  # Append fcs filepath to the list corresponding to that key
                cbid_dict.setdefault(cbid, []).append(fcs)  # append new value to list inside dict

            else:  # if no key exists for that CBID yet
                cbid_dict[cbid] = [fcs]

            cbid_list.append(cbid)  # Make list of CBIDs to use later as keys
    cbid_list = sorted(list(set(cbid_list)), key=int)
    if min_sample_id:
        cbid_list = [x for x in cbid_list if int(x) >= int(min_sample_id)]
    if max_sample_id:
        cbid_list = [x for x in cbid_list if int(x) <= int(max_sample_id)]

    final_filepaths = []
    repeat_filepaths = []
    unique_sample_fcs_filepaths = []
    for cbid in cbid_list:
        if include_repeats and (len(cbid_dict[cbid]) >= 2):  # If repeats are included and this CBID has more than 2
            # samples associated with it

            final_filepaths.append(cbid_dict[cbid][-2])  # Second to last element of the corresponding list
            final_filepaths.append(cbid_dict[cbid][-1])  # Last element of the corresponding list

            # Add repeats only to additional list
            repeat_filepaths.append(cbid_dict[cbid][-2])
        else:  # if no repeats or CBID only has single fcs file
            final_filepaths.append(cbid_dict[cbid][-1])  # Last element of the corresponding list
        unique_sample_fcs_filepaths.append(cbid_dict[cbid][-1])

    # NOTE: the latest fcs file for each CBID is considered to be the 'original', while older files are considered
    # repeats

    print('Available samples:\n\tAll .fcs files: ', len(fcs_filepaths))
    print('\tUnique sample .fcs files included: ', len(unique_sample_fcs_filepaths))
    print('\tNumber of .fcs files in the final list of files to analyze:', len(final_filepaths))
    print('\tNumber of .fcs files belonging to repeats included: ', len(repeat_filepaths))

    return final_filepaths


def get_out_folder_list(input_image_path_prefix: str) -> list:
    """
    get_out_folder_list returns a list of names of the folders where we intend to save
    the images resulting from the gating. The parent path of all of these folders
    can be edited.

    :param input_image_path_prefix: the path to the directory where our folder structure will be stored.
    
    :return: List of all output folders for the images
    """

    folder_list = [input_image_path_prefix + "01-Viable",
                   input_image_path_prefix + "02-Singlet",
                   input_image_path_prefix + "03-CBMCs",
                   input_image_path_prefix + "04-CD45high",
                   input_image_path_prefix + "05-CD45posCD34pos",
                   input_image_path_prefix + "06-CD3negBACKGATE",
                   input_image_path_prefix + "07-CD4CD8/07A-CD4posBACKGATE/",
                   input_image_path_prefix + "07-CD4CD8/07B-CD8posBACKGATE/",
                   input_image_path_prefix + "08-CD19negBACKGATE",

                   input_image_path_prefix + "09-CD14/09A-CD14negBACKGATE",
                   input_image_path_prefix + "09-CD14/09B-CD14posMonocytes",

                   input_image_path_prefix + "10-CD16CD56/10A-CD16neg_56negBACKGATE",
                   input_image_path_prefix + "10-CD16CD56/10B-CD16pospos_56posposBACKGATE",
                   input_image_path_prefix + "10-CD16CD56/10C-CD16neg_56posposBACKGATE",

                   input_image_path_prefix + "11-LinnegCD34pos",

                   input_image_path_prefix + "CD45RA_threshold",
                   input_image_path_prefix + "CD90_threshold",
                   input_image_path_prefix + "CD10_threshold",
                   input_image_path_prefix + "CD135_threshold",

                   input_image_path_prefix + "CD38-0.3/12-CD38neg",
                   input_image_path_prefix + "CD38-0.3/13-HSC_MLP_MPP",
                   input_image_path_prefix + "CD38-0.3/14-B_NK",
                   input_image_path_prefix + "CD38-0.3/15-CD10pos(MLP)",
                   input_image_path_prefix + "CD38-0.3/16-CMP_GMP_MEP",
                   input_image_path_prefix + "CD38-0.3/17-trueHSC",

                   input_image_path_prefix + "CD38-0.2/12-CD38neg",
                   input_image_path_prefix + "CD38-0.2/13-HSC_MLP_MPP",
                   input_image_path_prefix + "CD38-0.2/14-B_NK",
                   input_image_path_prefix + "CD38-0.2/15-CD10pos(MLP)",
                   input_image_path_prefix + "CD38-0.2/16-CMP_GMP_MEP",
                   input_image_path_prefix + "CD38-0.2/17-trueHSC",

                   input_image_path_prefix + "CD38-0.1/12-CD38neg",
                   input_image_path_prefix + "CD38-0.1/13-HSC_MLP_MPP",
                   input_image_path_prefix + "CD38-0.1/14-B_NK",
                   input_image_path_prefix + "CD38-0.1/15-CD10pos(MLP)",
                   input_image_path_prefix + "CD38-0.1/16-CMP_GMP_MEP",
                   input_image_path_prefix + "CD38-0.1/17-trueHSC",

                   input_image_path_prefix + "CD38-0.05/12-CD38neg",
                   input_image_path_prefix + "CD38-0.05/13-HSC_MLP_MPP",
                   input_image_path_prefix + "CD38-0.05/14-B_NK",
                   input_image_path_prefix + "CD38-0.05/15-CD10pos(MLP)",
                   input_image_path_prefix + "CD38-0.05/16-CMP_GMP_MEP",
                   input_image_path_prefix + "CD38-0.05/17-trueHSC"]

    return folder_list


def gate_full_dataset(my_sample, save_images=True):
    """
    gateFullDataset describes the whole gating strategy for a single sample. The AliGater experiment object
    uses its method 'apply' to apply this gating strategy to each of the desired fcs files. 
    
    :param my_sample: AliGater object built from the fcs file corresponding to a single sample
    :param save_images: a bool that determines if images of the gates will be saved
    
    :return: my_sample (AliGater object, now updated with the results from the gating)
    """

    # Setup
    ag.agconf.ag_verbose = False
    sample_filepath = my_sample.filePath
    date_plate = ag.getFileName(ag.getParent(sample_filepath))
    sample_name = ag.getFileName(sample_filepath)

    # Event cap.  We cap samples with too many events. Take only the first 3.5M
    event_cap = 3200000  # Used to be 3500000, but it kept crashing because the computer struggles with so many events
    if len(my_sample()) > event_cap:
        my_sample.fcsDF = my_sample.fcsDF.iloc[0:event_cap]
        sys.stderr.write("Event cap reached. Down sampled to " + str(event_cap) + " events.\n")

    # Marker check
    necessary_markers = ['FSC 488/10-H', 'FSC 488/10-A', 'SSC 488/10-A',
                         'CD45RA FITC-A', 'CD8 PerCP-Cy5.5-A', 'CD34 PE-Cy7-A',
                         'CD90 PE (R-phycoerythrin)-A', 'CD19 PE-Texas Red-A',
                         'CD56 (BV) 650-A', 'CD135 (BV) 711-A', 'CD16 (BV) 786-A',
                         'CD38 (BV) 421-A', 'CD14 (BV) 605-A', 'CD4 (BV) 510-A',
                         'CD3 APC-H7-A', 'CD10 APC (Allophycocyanin)-A']

    # Put inconsistent marker names into variables
    markers = my_sample().columns.tolist()

    for marker in necessary_markers:
        if marker not in markers:  # check all necessary markers are present in file
            print(f"WARNING: {marker} couldn't be found in file {sample_name}. Skipping.")
            return my_sample

    # A few of the marker names are inconsistent across files. We'll fix them here.
    marker_cd45 = ''
    marker_7aad = ''
    for marker in markers[6:]:  # For marker names that tend to change between samples, store the name in a variable
        if marker.startswith('CD45'):  # it's a miracle this works...
            marker_cd45 = marker
        if marker.startswith('7AAD') or marker.startswith('CD235'):
            marker_7aad = marker

    # List of antibodies. Used when getting MFI values
    antibody_list = [x for x in my_sample().columns.to_list() if not x.startswith('FSC') if not x.startswith('SSC') if
                     x.startswith('CD')]

    # Image save switch
    if save_images:
        file_name = 'notNone'  # dummy variable to prevent errors
    else:
        file_name = None

    # Alternate gating strategy switch. If the sample is in the 'missgated_gate4' list, apply alternative strategy
    alternate_gate4_strategy = False
    with open('../output/AliGater_gating_QC/w49_2021/missgated_gate4_xlim_tooFarBack_w49.txt', 'r') as in_file:
        for line in in_file:
            if my_sample.filePath.endswith(line.rstrip() + '.fcs'):
                alternate_gate4_strategy = True
                print(f"Sample {line.rstrip()} .fcs belongs to missgate list for gate 4 (CD34+ gate) with incorrect "
                      f"xlim. Applying alternate gate 4 strategy.")

    # Read the lists of samples that require adjusted values
    path_to_adjust_files = '/home/antton/Projects/CB_Data_Analysis/data/processed/AliGater_gating_corrections/'
    CD45RA_BLOB_list = []
    with open(path_to_adjust_files + 'CD45RA_BLOB_list.txt', 'r') as in_file:
        for line in in_file:
            CD45RA_BLOB_list.append(line.rstrip())
    gate1_ylim_400_list = []
    with open(path_to_adjust_files + 'gate1_ylim_400.txt', 'r') as in_file:
        for line in in_file:
            gate1_ylim_400_list.append(line.rstrip())
    gate1_ylim_600_list = []
    with open(path_to_adjust_files + 'gate1_ylim_600.txt', 'r') as in_file:
        for line in in_file:
            gate1_ylim_600_list.append(line.rstrip())
    gate1_ylim_1000_list = []
    with open(path_to_adjust_files + 'gate1_ylim_1000.txt', 'r') as in_file:
        for line in in_file:
            gate1_ylim_1000_list.append(line.rstrip())
    gate1_ylim_2000_list = []
    with open(path_to_adjust_files + 'gate1_ylim_2000.txt', 'r') as in_file:
        for line in in_file:
            gate1_ylim_2000_list.append(line.rstrip())
    gate2_witdthScale_5_list = []
    with open(path_to_adjust_files + 'gate2-widthScale_5.txt', 'r') as in_file:
        for line in in_file:
            gate2_witdthScale_5_list.append(line.rstrip())
    gate3_xboundary_45000_list = []
    with open(path_to_adjust_files + 'gate3-xboundary_45000.txt', 'r') as in_file:
        for line in in_file:
            gate3_xboundary_45000_list.append(line.rstrip())
    gate3_ylim_bot_8000_list = []
    with open(path_to_adjust_files + 'gate3-ylim_bot_8000.txt', 'r') as in_file:
        for line in in_file:
            gate3_ylim_bot_8000_list.append(line.rstrip())
    gate4_maxVal_3000_list = []
    with open(path_to_adjust_files + 'gate4-maxVal_3000.txt', 'r') as in_file:
        for line in in_file:
            gate4_maxVal_3000_list.append(line.rstrip())
    gate4_maxVal_7000_list = []
    with open(path_to_adjust_files + 'gate4-maxVal_7000.txt', 'r') as in_file:
        for line in in_file:
            gate4_maxVal_7000_list.append(line.rstrip())
    gate4_xlim_1000_list = []
    with open(path_to_adjust_files + 'gate4-xlim_1000.txt', 'r') as in_file:
        for line in in_file:
            gate4_xlim_1000_list.append(line.rstrip())
    gate4_xlim_4000_list = []
    with open(path_to_adjust_files + 'gate4-xlim_4000.txt', 'r') as in_file:
        for line in in_file:
            gate4_xlim_4000_list.append(line.rstrip())
    gate6_ylim_1000_list = []
    with open(path_to_adjust_files + 'gate6-ylim_1000.txt', 'r') as in_file:
        for line in in_file:
            gate6_ylim_1000_list.append(line.rstrip())
    gate7_xlim_800_list = []
    with open(path_to_adjust_files + 'gate7-xlim_800.txt', 'r') as in_file:
        for line in in_file:
            gate7_xlim_800_list.append(line.rstrip())
    gate7_xlim_1000_list = []
    with open(path_to_adjust_files + 'gate7-xlim_1000.txt', 'r') as in_file:
        for line in in_file:
            gate7_xlim_1000_list.append(line.rstrip())
    gate7_ylim_400_list = []
    with open(path_to_adjust_files + 'gate7-ylim_400.txt', 'r') as in_file:
        for line in in_file:
            gate7_ylim_400_list.append(line.rstrip())
    gate7_ylim_1000_list = []
    with open(path_to_adjust_files + 'gate7-ylim_1000.txt', 'r') as in_file:
        for line in in_file:
            gate7_ylim_1000_list.append(line.rstrip())
    gate9_ylim_500_list = []
    with open(path_to_adjust_files + 'gate9-ylim_500.txt', 'r') as in_file:
        for line in in_file:
            gate9_ylim_500_list.append(line.rstrip())
    gate9_ylim_800_list = []
    with open(path_to_adjust_files + 'gate9-ylim_800.txt', 'r') as in_file:
        for line in in_file:
            gate9_ylim_800_list.append(line.rstrip())
    gate9_ylim_1100_list = []
    with open(path_to_adjust_files + 'gate9-ylim_1100.txt', 'r') as in_file:
        for line in in_file:
            gate9_ylim_1100_list.append(line.rstrip())
    gate9_ylim_2000_list = []
    with open(path_to_adjust_files + 'gate9-ylim_2000.txt', 'r') as in_file:
        for line in in_file:
            gate9_ylim_2000_list.append(line.rstrip())
    gate10_xlim_1500_list = []
    with open(path_to_adjust_files + 'gate10-xlim_1500.txt', 'r') as in_file:
        for line in in_file:
            gate10_xlim_1500_list.append(line.rstrip())
    gate10_xlim_3000_list = []
    with open(path_to_adjust_files + 'gate10-xlim_3000.txt', 'r') as in_file:
        for line in in_file:
            gate10_xlim_3000_list.append(line.rstrip())
    gate10_ylim_500_list = []
    with open(path_to_adjust_files + 'gate10-ylim_500.txt', 'r') as in_file:
        for line in in_file:
            gate10_ylim_500_list.append(line.rstrip())
    CD45RA_threshold_100_list = []
    with open(path_to_adjust_files + 'CD45RA_threshold_100.txt', 'r') as in_file:
        for line in in_file:
            CD45RA_threshold_100_list.append(line.rstrip())
    CD45RA_threshold_120_list = []
    with open(path_to_adjust_files + 'CD45RA_threshold_120.txt', 'r') as in_file:
        for line in in_file:
            CD45RA_threshold_120_list.append(line.rstrip())
    CD45RA_threshold_200_list = []
    with open(path_to_adjust_files + 'CD45RA_threshold_200.txt', 'r') as in_file:
        for line in in_file:
            CD45RA_threshold_200_list.append(line.rstrip())
    CD45RA_threshold_280_list = []
    with open(path_to_adjust_files + 'CD45RA_threshold_280.txt', 'r') as in_file:
        for line in in_file:
            CD45RA_threshold_280_list.append(line.rstrip())
    CD90_threshold_1000_list = []
    with open(path_to_adjust_files + 'CD90_threshold_1000.txt', 'r') as in_file:
        for line in in_file:
            CD90_threshold_1000_list.append(line.rstrip())
    CD90_threshold_2000_list = []
    with open(path_to_adjust_files + 'CD90_threshold_2000.txt', 'r') as in_file:
        for line in in_file:
            CD90_threshold_2000_list.append(line.rstrip())
    CD90_threshold_5000_list = []
    with open(path_to_adjust_files + 'CD90_threshold_5000.txt', 'r') as in_file:
        for line in in_file:
            CD90_threshold_5000_list.append(line.rstrip())

    # GATING STARTS HERE

    ####################################################################################################################
    # Gate 1:  Viable/non-viable separation through 7AAD
    ####################################################################################################################

    livedead_step1 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol=marker_7aad,
                                      scale='linear', T=200, yscale='bilog', thresh=214000, parentGate=None,
                                      orientation='vertical', population='lower')
    livedead_step2 = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol=marker_7aad,
                                      scale='linear', T=200, yscale='bilog', thresh=-800, parentGate=livedead_step1,
                                      orientation='horizontal', population='upper')
    halfcut = ag.gateThreshold(my_sample, name="remove_clutter", xCol='FSC 488/10-A', yCol=marker_7aad, scale='linear',
                               T=200, yscale='bilog', thresh=150000, parentGate=livedead_step2, orientation='vertical',
                               population='upper')

    if f'fcs/{date_plate}/{sample_name}' in gate1_ylim_400_list:
        ylim_back = 400
    elif f'fcs/{date_plate}/{sample_name}' in gate1_ylim_600_list:
        ylim_back = 600
    elif f'fcs/{date_plate}/{sample_name}' in gate1_ylim_1000_list:
        ylim_back = 1000
    elif f'fcs/{date_plate}/{sample_name}' in gate1_ylim_2000_list:
        ylim_back = 2000
    else:
        ylim_back = ag.densityDelimitation(my_sample, xCol=marker_7aad, parentGate=halfcut, interval=[300, 1100],
                                           limit_threshold=0.2, direction='right', scale='bilog', T=200)
    if ylim_back == inf:  # Failsafe in case of no limit
        ylim_back = 900
    if file_name:
        file_name = image_path_prefix + "01-Viable/" + date_plate + "-" + sample_name + "-Viable.jpeg"
    livedead_final = ag.gateThreshold(my_sample, name="Viable", xCol='FSC 488/10-A', yCol=marker_7aad, scale='linear',
                                      T=200, yscale='logicle', thresh=ylim_back, parentGate=livedead_step2,
                                      orientation='horizontal', population='lower', filePlot=file_name)

    ####################################################################################################################
    # Gate 2:  Singlets
    ####################################################################################################################

    if f'fcs/{date_plate}/{sample_name}' in gate2_witdthScale_5_list:
        widthScale = 5
    else:
        widthScale = 2.5

    singlets_step1 = ag.gateThreshold(my_sample, name="remove_high_clutter", xCol='FSC 488/10-A', yCol='FSC 488/10-H',
                                      scale='linear', thresh=210000, parentGate=livedead_final, orientation='vertical',
                                      population='lower')
    if file_name:
        file_name = image_path_prefix + "02-Singlet/" + date_plate + "-" + sample_name + "-Singlet.jpeg"
    singlets = ag.gatePC(my_sample, name="Singlets", xCol='FSC 488/10-A', yCol='FSC 488/10-H', center='centroid',
                         adjustAngle=1, widthScale=widthScale, heightScale=3.5, parentGate=singlets_step1,
                         filePlot=file_name)

    ####################################################################################################################
    # Gate 3: CBMC cells
    ####################################################################################################################

    halfcut_middle = ag.gateThreshold(my_sample, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A',
                                      scale='linear', thresh=120000, parentGate=singlets, orientation='vertical',
                                      population='upper')

    ylim_bot = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_middle,
                                      interval=[10000, 20000], limit_threshold=0.05, direction='left', scale='linear')
    if ylim_bot == inf:
        ylim_bot = 10000

    if f'fcs/{date_plate}/{sample_name}' in gate3_ylim_bot_8000_list:
        ylim_bot = 8000

    CBMC_step1 = ag.gateCorner(my_sample, name="cut_corner", xCol='FSC 488/10-A', yCol='SSC 488/10-A', xThresh=85000,
                               yThresh=ylim_bot, xOrientation='lower', yOrientation='lower', Outer=True,
                               parentGate=singlets)

    halfcut_tail = ag.gateThreshold(my_sample, name="right_tail", xCol='FSC 488/10-A', yCol='SSC 488/10-A',
                                    scale='linear', thresh=180000, parentGate=singlets, orientation='vertical',
                                    population='upper')
    ylim_top = ag.densityDelimitation(my_sample, xCol='SSC 488/10-A', parentGate=halfcut_tail, interval=[50000, 125000],
                                      limit_threshold=0.2, direction='right', scale='linear') + 25000
    if ylim_top == inf:
        ylim_top = 140000

    if f'fcs/{date_plate}/{sample_name}' in gate3_xboundary_45000_list:
        xboundary_low = 45000
    else:
        xboundary_low = 55000

    cbmc = ag.horizontalPath(my_sample, name="CBMC",
                             xCol='FSC 488/10-A', yCol='SSC 488/10-A', population='lower',
                             startY=ylim_bot + 30000, endY=ylim_top, xboundaries=[xboundary_low, 140000],
                             yboundaries=[ylim_bot + 25000, ylim_top + 5000],
                             leftRight=True, direction='both',
                             maxStep=2, phi=0.1, bins=100, sigma=1,
                             scale='linear', parentGate=CBMC_step1)

    if file_name:
        file_name = image_path_prefix + "03-CBMCs/" + date_plate + "-" + sample_name + "-CBMCs.jpeg"

    ag.backGate(my_sample, population=cbmc, background_population=singlets, xCol='FSC 488/10-A', yCol='SSC 488/10-A',
                scale='linear', T=200, markersize=0.1, filePlot=file_name)

    my_sample.update(ag.AGgate(cbmc, None, 'FSC 488/10-A', 'SSC 488/10-A', "CBMC"), QC=True, scale='linear',
                     xlim=[25000, 180000], ylim=[0, 200000])

    ####################################################################################################################
    # Gate 4: CD45+ AND CD34+ out of CBMC cells, "Declutter gate" to get rid of CD45- debree
    ####################################################################################################################

    xlim_arbitrary = 11000
    ylim_arbitrary = 20000
    # Take right side and calculate maximum density point in the y-axis (the location of the CD34 cluster)
    right_half = ag.gateCorner(my_sample, name="right_half", xCol='CD34 PE-Cy7-A', yCol=marker_cd45, scale='bilog',
                               T=200,
                               xThresh=xlim_arbitrary, yThresh=ylim_arbitrary, parentGate=cbmc, xOrientation='upper',
                               yOrientation='lower')
    # Get the highest density point of CD34 blob to use as reference
    if f'fcs/{date_plate}/{sample_name}' in gate4_maxVal_3000_list:
        maxVal = 3000
    elif f'fcs/{date_plate}/{sample_name}' in gate4_maxVal_7000_list:
        maxVal = 7000
    else:
        mean, median, sigma, maxVal = ag.axisStats(my_sample(), xCol=marker_cd45, vI=right_half(),
                                                   bins=10000)  # It is important to have a high number of bins for the
    # estimation to work correctly
    # Take a cut of the upper half that includes the CD34 cluster, so as tu find where they separate in the x-axis
    upper_half = ag.gateThreshold(my_sample, name="upper_half", xCol='CD34 PE-Cy7-A', yCol=marker_cd45, scale='bilog',
                                  T=200,
                                  thresh=maxVal, parentGate=cbmc, orientation='horizontal', population='upper')
    y_correction_amount = 3000
    lower_half_of_upper_half = ag.gateThreshold(my_sample, name="lower_half_of_upper_half", xCol='CD34 PE-Cy7-A',
                                                yCol=marker_cd45, scale='bilog', T=200,
                                                thresh=maxVal + y_correction_amount,
                                                parentGate=upper_half, orientation='horizontal', population='lower')
    x_correction_amount = 2000  # value used to offset the calculated point
    xlim_middle0 = ag.valleySeek(my_sample, xCol='CD34 PE-Cy7-A', interval=[800, 11000], require_local_min=True,
                                 scale='bilog', T=200,
                                 parentGate=lower_half_of_upper_half)
    if alternate_gate4_strategy:
        xlim_middle0 += x_correction_amount  # add the correction back up, since these samples don't seem to need it
    if f'fcs/{date_plate}/{sample_name}' in gate4_xlim_1000_list:
        xlim_middle0 = 1000
    elif f'fcs/{date_plate}/{sample_name}' in gate4_xlim_4000_list:
        xlim_middle0 = 4000

    if file_name:
        file_name = image_path_prefix + "04-CD45high/" + date_plate + "-" + sample_name + "-CD45high.jpeg"
    CD45_high = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A', yCol=marker_cd45,
                                  startPoint=(xlim_middle0, maxVal),
                                  endLimits=(None, maxVal - 1.5 * y_correction_amount), theta=-65, scale='bilog',
                                  T=200, population='upper', parentGate=cbmc, filePlot=file_name)

    ####################################################################################################################
    # Gate 5: CD34+ ratio relative to CD45+
    ####################################################################################################################

    CD34pos = ag.gateTiltedLine(my_sample, name="CD34+", xCol='CD34 PE-Cy7-A', yCol=marker_cd45,
                                startPoint=(xlim_middle0, maxVal), endLimits=(None, maxVal + 7500), theta=40,
                                scale='bilog', T=200, population='lower', parentGate=CD45_high)
    if file_name:
        file_name = image_path_prefix + "05-CD45posCD34pos/" + date_plate + "-" + sample_name + "-CD45posCD34pos.jpeg"
    CD45pos = ag.gateTiltedLine(my_sample, name="CD45+", xCol='CD34 PE-Cy7-A', yCol=marker_cd45,
                                startPoint=(xlim_middle0, maxVal), endLimits=(None, maxVal + 7500), theta=40,
                                scale='bilog', T=200, population='upper', parentGate=CD45_high, filePlot=file_name)

    my_sample.update(ag.AGgate(CD34pos, CD45pos, 'CD34 PE-Cy7-A', marker_cd45, "CD34+"), MFI=True,
                     extra_MFI=antibody_list,
                     QC=True, scale='bilog', xscale='bilog', yscale='bilog', T=200, xlim=[-1000, 100000],
                     ylim=[1000, 150000])

    ####################################################################################################################
    # Gate 6: CD3 out of CD45high
    ####################################################################################################################

    # Gate 6a: CD3-
    if f'fcs/{date_plate}/{sample_name}' in gate6_ylim_1000_list:
        ylim_middle = 1000
    else:
        ylim_middle = ag.valleySeek(my_sample, xCol='CD3 APC-H7-A', interval=[-1000, 7000], require_local_min=True,
                                    scale='bilog', T=200, parentGate=CD45_high)
    if ylim_middle == np.Infinity:
        ylim_middle = 1000
    CD3neg = ag.gateThreshold(my_sample, name="CD3-", xCol=marker_cd45, yCol='CD3 APC-H7-A', scale='bilog', T=200,
                              thresh=ylim_middle, parentGate=CD45_high, orientation='horizontal', population='lower')
    if file_name:
        file_name = image_path_prefix + "06-CD3negBACKGATE/" + date_plate + "-" + sample_name + "-CD3negBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD3neg, background_population=cbmc, xCol=marker_cd45, yCol='CD3 APC-H7-A',
                scale='bilog',
                T=200, markersize=0.1, filePlot=file_name)
    # Gate 6b: CD3+
    CD3pos = ag.gateThreshold(my_sample, name="CD3+", xCol=marker_cd45, yCol='CD3 APC-H7-A', scale='bilog', T=200,
                              thresh=ylim_middle, parentGate=CD45_high, orientation='horizontal', population='upper')

    my_sample.update(ag.AGgate(CD3pos, CD45pos, marker_cd45, 'CD3 APC-H7-A', "CD3+"), MFI=True, extra_MFI=antibody_list,
                     QC=False)

    ####################################################################################################################
    # Gate 7: CD4, CD8 out of CD3+
    ####################################################################################################################

    # Gate 7a: CD4+
    xlim_middle = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True,
                                scale='bilog', T=200, parentGate=CD3pos)
    if xlim_middle == np.Infinity:
        xlim_middle = 2000

    if f'fcs/{date_plate}/{sample_name}' in gate7_xlim_800_list:
        xlim_middle = 800
    elif f'fcs/{date_plate}/{sample_name}' in gate7_xlim_1000_list:
        xlim_middle = 1000

    cd4cd8_step1 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A', yCol='CD8 PerCP-Cy5.5-A',
                                    scale='bilog', T=200, thresh=xlim_middle, parentGate=CD3pos, orientation='vertical',
                                    population='upper')
    CD4pos_final = ag.gatePC(my_sample, name="CD4+", xCol='CD4 (BV) 510-A', yCol='CD8 PerCP-Cy5.5-A', center='centroid',
                             adjustAngle=3, widthScale=2.5, scale='bilog', T=200, heightScale=3.5,
                             parentGate=cd4cd8_step1)
    if file_name:
        file_name = image_path_prefix + "07-CD4CD8/07A-CD4posBACKGATE/" + date_plate + "-" + sample_name + "-CD4posBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD4pos_final, background_population=CD3pos, xCol='CD4 (BV) 510-A',
                yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)

    my_sample.update(ag.AGgate(CD4pos_final, CD3pos, 'CD4 (BV) 510-A', 'CD8 PerCP-Cy5.5-A', "CD4+_1"), MFI=True,
                     extra_MFI=antibody_list, QC=False)

    # Gate 7b: CD8+
    ylim_middle = ag.valleySeek(my_sample, xCol='CD8 PerCP-Cy5.5-A', interval=[200, 2000], require_local_min=True,
                                parentGate=CD3pos)
    if ylim_middle == np.Infinity:
        ylim_middle = 1000

    if f'fcs/{date_plate}/{sample_name}' in gate7_ylim_400_list:
        ylim_middle = 400
    elif f'fcs/{date_plate}/{sample_name}' in gate7_ylim_1000_list:
        ylim_middle = 1000

    cd4cd8_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD4 (BV) 510-A', yCol='CD8 PerCP-Cy5.5-A',
                                    scale='bilog', T=200, thresh=ylim_middle, parentGate=CD3pos,
                                    orientation='horizontal', population='upper')
    xlim_middle_upper = ag.valleySeek(my_sample, xCol='CD4 (BV) 510-A', interval=[1000, 5000], require_local_min=True,
                                      scale='bilog', T=200, parentGate=cd4cd8_step2)
    if xlim_middle_upper == np.Infinity:
        xlim_middle_upper = 2000
    cd4cd8_step3 = ag.gateThreshold(my_sample, name="separate_middle_x_axis", xCol='CD4 (BV) 510-A',
                                    yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, thresh=xlim_middle_upper,
                                    parentGate=cd4cd8_step2, orientation='vertical', population='lower')
    CD8pos_final = ag.gatePC(my_sample, name="CD8+", xCol='CD4 (BV) 510-A', yCol='CD8 PerCP-Cy5.5-A', center='centroid',
                             adjustAngle=0, widthScale=2, scale='bilog', T=200, heightScale=3.5,
                             parentGate=cd4cd8_step3)

    if file_name:
        file_name = image_path_prefix + "07-CD4CD8/07B-CD8posBACKGATE/" + date_plate + "-" + sample_name + "-CD8posBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD8pos_final, background_population=CD3pos, xCol='CD4 (BV) 510-A',
                yCol='CD8 PerCP-Cy5.5-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)

    my_sample.update(ag.AGgate(CD8pos_final, CD3pos, 'CD4 (BV) 510-A', 'CD8 PerCP-Cy5.5-A', "CD8+"), MFI=True,
                     extra_MFI=antibody_list, QC=False)
    my_sample.update(ag.AGgate(CD4pos_final, CD8pos_final, 'CD4 (BV) 510-A', 'CD8 PerCP-Cy5.5-A', "CD4+_2"), MFI=False,
                     QC=False)

    ####################################################################################################################
    # Gate 8: CD19 out of CD3-
    ####################################################################################################################

    # Gate 8a: CD19-
    ylim_middle_cd19 = ag.valleySeek(my_sample, xCol='CD19 PE-Texas Red-A', interval=[200, 10000],
                                     require_local_min=True, scale='bilog', T=200, parentGate=CD3neg)
    if ylim_middle_cd19 == np.Infinity:
        ylim_middle_cd19 = 1000
    CD19neg = ag.gateThreshold(my_sample, name="CD19-", xCol='CD3 APC-H7-A', yCol='CD19 PE-Texas Red-A', scale='bilog',
                               T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horizontal',
                               population='lower')
    if file_name:
        file_name = image_path_prefix + "08-CD19negBACKGATE/" + date_plate + "-" + sample_name + "-CD19negBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD19neg, background_population=CD3neg, xCol='CD3 APC-H7-A',
                yCol='CD19 PE-Texas Red-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)
    # Gate 8b: CD19+
    CD19pos = ag.gateThreshold(my_sample, name="CD19+", xCol='CD3 APC-H7-A', yCol='CD19 PE-Texas Red-A', scale='bilog',
                               T=200, thresh=ylim_middle_cd19, parentGate=CD3neg, orientation='horizontal',
                               population='upper', filePlot=file_name)

    my_sample.update(ag.AGgate(CD19pos, CD3neg, 'CD3 APC-H7-A', 'CD19 PE-Texas Red-A', "CD19+_1"), MFI=True,
                     extra_MFI=antibody_list, QC=False)

    my_sample.update(ag.AGgate(CD19pos, CD45pos, 'CD3 APC-H7-A', 'CD19 PE-Texas Red-A', "CD19+_2"), MFI=True,
                     extra_MFI=antibody_list, QC=False)

    ####################################################################################################################
    # Gate 9: CD14 out of CD3-/CD19- (Monocytes)
    ####################################################################################################################

    # Gate 9a: CD14- out of CD19-
    ylim_middle_cd14 = ag.valleySeek(my_sample, xCol='CD14 (BV) 605-A', interval=[0, 3000], require_local_min=True,
                                     parentGate=CD19neg, scale='bilog')
    if ylim_middle_cd14 == np.Infinity:
        ylim_middle_cd14 = 800

    if f'fcs/{date_plate}/{sample_name}/' in gate9_ylim_500_list:
        ylim_middle_cd14 = 500
    if f'fcs/{date_plate}/{sample_name}/' in gate9_ylim_800_list:
        ylim_middle_cd14 = 800
    elif f'fcs/{date_plate}/{sample_name}/' in gate9_ylim_1100_list:
        ylim_middle_cd14 = 1100
    elif f'fcs/{date_plate}/{sample_name}/' in gate9_ylim_2000_list:
        ylim_middle_cd14 = 2000
    CD14neg = ag.gateThreshold(my_sample, name="CD14-", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A',
                               scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg,
                               orientation='horizontal', population='lower')
    if file_name:
        file_name = image_path_prefix + "09-CD14/09A-CD14negBACKGATE/" + date_plate + "-" + sample_name + \
                    "-CD14negBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD14neg, background_population=CD19neg, xCol='CD16 (BV) 786-A',
                yCol='CD14 (BV) 605-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)
    # Gate 9b: CD14+ out of CD19- (All monocytes)
    CD14pos = ag.gateThreshold(my_sample, name="CD14+", xCol='CD19 PE-Texas Red-A', yCol='CD14 (BV) 605-A',
                               scale='bilog', T=200, thresh=ylim_middle_cd14, parentGate=CD19neg,
                               orientation='horizontal', population='upper')

    # out of parent
    my_sample.update(ag.AGgate(CD14pos, CD19neg, 'CD19 PE-Texas Red-A', 'CD14 (BV) 605-A', "CD14+_1"), MFI=True,
                     extra_MFI=antibody_list, QC=False)

    my_sample.update(ag.AGgate(CD14pos, CD45pos, 'CD19 PE-Texas Red-A', 'CD14 (BV) 605-A', "CD14+_2"), MFI=True,
                     extra_MFI=antibody_list, QC=False)

    ####################################################################################################################
    # Gate 10: CD16 and CD56 values out of CD14-
    ####################################################################################################################

    aprox_xlim_middle = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True,
                                      scale='bilog', T=200, parentGate=CD14neg)
    if aprox_xlim_middle == np.Infinity:
        aprox_xlim_middle = 2000
    left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A', yCol='CD56 (BV) 650-A',
                                 scale='bilog', T=200, thresh=aprox_xlim_middle, parentGate=CD14neg,
                                 orientation='vertical', population='lower')
    aprox_ylim_middle = ag.valleySeek(my_sample, xCol='CD56 (BV) 650-A', interval=[700, 1300], require_local_min=True,
                                      scale='bilog', T=200, parentGate=CD14neg)
    if aprox_ylim_middle == np.Infinity:
        aprox_ylim_middle = 1000
    lower_left_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A',
                                       yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_ylim_middle,
                                       parentGate=left_half, orientation='horizontal', population='lower')
    main_ylim_middle = ag.densityDelimitation(my_sample, xCol='CD56 (BV) 650-A', parentGate=lower_left_half,
                                              interval=[0, 2000], limit_threshold=0.07, direction='right',
                                              scale='linear')
    if main_ylim_middle == inf:
        main_ylim_middle = aprox_ylim_middle
    if f'fcs/{date_plate}/{sample_name}' in gate10_ylim_500_list:
        main_ylim_middle = 500

    lower_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A', yCol='CD56 (BV) 650-A',
                                  scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg,
                                  orientation='horizontal', population='lower')
    upper_half = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A', yCol='CD56 (BV) 650-A',
                                  scale='bilog', T=200, thresh=main_ylim_middle, parentGate=CD14neg,
                                  orientation='horizontal', population='upper')
    aprox_xlim_middle2 = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A', interval=[900, 18000], require_local_min=True,
                                       scale='bilog', T=200, parentGate=lower_half)
    if aprox_xlim_middle2 == np.Infinity:
        aprox_xlim_middle2 = 2000
    CD16neg_56neg_step1 = ag.gateThreshold(my_sample, name="lower_left_side", xCol='CD16 (BV) 786-A',
                                           yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2,
                                           parentGate=lower_half, orientation='vertical', population='lower')
    main_xlim_middle = ag.densityDelimitation(my_sample, xCol='CD16 (BV) 786-A', parentGate=CD16neg_56neg_step1,
                                              interval=[0, 2000], limit_threshold=0.07, direction='right',
                                              scale='linear')
    if main_xlim_middle == inf:
        main_xlim_middle = aprox_xlim_middle
    if f'fcs/{date_plate}/{sample_name}' in gate10_xlim_1500_list:
        main_xlim_middle = 1500
    elif f'fcs/{date_plate}/{sample_name}' in gate10_xlim_3000_list:
        main_xlim_middle = 3000

    CD16neg_56neg_step2 = ag.gateThreshold(my_sample, name="separate_middle", xCol='CD16 (BV) 786-A',
                                           yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle,
                                           parentGate=CD16neg_56neg_step1, orientation='vertical', population='lower')
    # CD16-CD56-, a.k.a. LIN NEGATIVE
    CD16neg_56neg_final = ag.gatePC(my_sample, name="CD16-CD56-", xCol='CD16 (BV) 786-A', yCol='CD56 (BV) 650-A',
                                    center='centroid', adjustAngle=0, widthScale=2, scale='bilog', T=200,
                                    heightScale=3.5, parentGate=CD16neg_56neg_step2)
    if file_name:
        file_name = image_path_prefix + "10-CD16CD56/10A-CD16neg_56negBACKGATE/" + date_plate + "-" + sample_name \
                    + "-CD16neg_56negBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16neg_56neg_final, background_population=CD14neg, xCol='CD16 (BV) 786-A',
                yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)
    central_section = ag.horizontalPath(my_sample, name="hor_path", xCol='CD16 (BV) 786-A', yCol='CD56 (BV) 650-A',
                                        population='lower', startY=3000, endY=12000, xboundaries=[0, 500000],
                                        yboundaries=[2000, 15000], leftRight=True, direction='both', maxStep=2, phi=0.1,
                                        bins=100, sigma=1, scale='bilog', T=200, parentGate=upper_half)
    # CD16++CD56++
    CD16pospos_56pospos_final = ag.gateThreshold(my_sample, name="CD16+CD56+", xCol='CD16 (BV) 786-A',
                                                 yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=main_xlim_middle,
                                                 parentGate=central_section, orientation='vertical', population='upper')
    if file_name:
        file_name = image_path_prefix + "10-CD16CD56/10B-CD16pospos_56posposBACKGATE/" + date_plate + "-" \
                    + sample_name + "-CD16pospos_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16pospos_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A',
                yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)

    my_sample.update(
        ag.AGgate(CD16pospos_56pospos_final, CD14neg, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16++CD56++_1"), MFI=True,
        extra_MFI=antibody_list, QC=False)

    my_sample.update(
        ag.AGgate(CD16pospos_56pospos_final, CD45pos, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16++CD56++_2"), MFI=True,
        extra_MFI=antibody_list, QC=False)

    # CD16-CD56++
    CD16neg_56pospos_final = ag.horizontalPath(my_sample, name="CD16-CD56++", xCol='CD16 (BV) 786-A',
                                               yCol='CD56 (BV) 650-A', population='upper', startY=3000, endY=12000,
                                               xboundaries=[0, 500000], yboundaries=[2000, 15000], leftRight=True,
                                               direction='both', maxStep=2, phi=0.1, bins=100, sigma=1, scale='bilog',
                                               T=200, parentGate=upper_half)
    if file_name:
        file_name = image_path_prefix + "10-CD16CD56/10C-CD16neg_56posposBACKGATE/" + date_plate + "-" + sample_name \
                    + "-CD16neg_56posposBACKGATE.jpeg"
    ag.backGate(my_sample, population=CD16neg_56pospos_final, background_population=CD14neg, xCol='CD16 (BV) 786-A',
                yCol='CD56 (BV) 650-A', scale='bilog', T=200, markersize=0.1, filePlot=file_name)

    my_sample.update(ag.AGgate(CD16neg_56pospos_final, CD14neg, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16-CD56++_1"),
                     MFI=True, extra_MFI=antibody_list, QC=False)
    my_sample.update(ag.AGgate(CD16neg_56pospos_final, CD45pos, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16-CD56++_2"),
                     MFI=True, extra_MFI=antibody_list, QC=False)

    # CD16++CD56-
    CD16pospos_56neg_final = ag.gateThreshold(my_sample, name="CD16++CD56-", xCol='CD16 (BV) 786-A',
                                              yCol='CD56 (BV) 650-A', scale='bilog', T=200, thresh=aprox_xlim_middle2,
                                              parentGate=lower_half, orientation='vertical', population='upper')

    my_sample.update(ag.AGgate(CD16pospos_56neg_final, CD14neg, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16++CD56-_1"),
                     MFI=True, extra_MFI=antibody_list, QC=False)
    my_sample.update(ag.AGgate(CD16pospos_56neg_final, CD45pos, 'CD16 (BV) 786-A', 'CD56 (BV) 650-A', "CD16++CD56-_2"),
                     MFI=True, extra_MFI=antibody_list, QC=False)

    ####################################################################################################################
    # Gate10b : Classical vs. Non-classical Monocytes

    if file_name:
        file_name = image_path_prefix + "09-CD14/09B-CD14posMonocytes/" + date_plate + "-" + sample_name + \
                    "-CD14posMonocytes.jpeg"

    CD16_threshold = ag.valleySeek(my_sample, xCol='CD16 (BV) 786-A',
                                   interval=[main_xlim_middle, main_xlim_middle + 10000], require_local_min=True,
                                   scale='bilog', T=200, parentGate=CD14pos)
    if CD16_threshold == np.Infinity:
        CD16_threshold = main_xlim_middle + 5000
    # Classical Monocytes
    CD14posCD16neg = ag.gateThreshold(my_sample, name="CD14+CD16-", xCol='CD16 (BV) 786-A', yCol='CD14 (BV) 605-A',
                                      scale='bilog', T=200, thresh=CD16_threshold, parentGate=CD14pos,
                                      orientation='vertical', population='lower', filePlot=file_name)
    my_sample.update(ag.AGgate(CD14posCD16neg, CD14pos, 'CD16 (BV) 786-A', 'CD14 (BV) 605-A', "CD14+CD16-_1"),
                     MFI=True, extra_MFI=antibody_list, QC=False)
    my_sample.update(ag.AGgate(CD14posCD16neg, CD45pos, 'CD16 (BV) 786-A', 'CD14 (BV) 605-A', "CD14+CD16-_2"),
                     MFI=False, QC=False)

    # Non-classical Monocytes
    CD14posCD16pos = ag.gateThreshold(my_sample, name="CD14+CD16-", xCol='CD16 (BV) 786-A', yCol='CD14 (BV) 605-A',
                                      scale='bilog', T=200, thresh=CD16_threshold, parentGate=CD14pos,
                                      orientation='vertical', population='upper')
    my_sample.update(ag.AGgate(CD14posCD16pos, CD14pos, 'CD16 (BV) 786-A', 'CD14 (BV) 605-A', "CD14+CD16+_1"),
                     MFI=True, extra_MFI=antibody_list, QC=False)
    my_sample.update(ag.AGgate(CD14posCD16pos, CD45pos, 'CD16 (BV) 786-A', 'CD14 (BV) 605-A', "CD14+CD16+_2"),
                     MFI=False, QC=False)

    ####################################################################################################################
    # Gate 11: Lin-CD34+ out of CD16-CD56-/CD14-/CD19-/CD3-/CD45+
    ####################################################################################################################

    if file_name:
        file_name = image_path_prefix + "11-LinnegCD34pos/" + date_plate + "-" + sample_name + "-LinnegCD34pos.jpeg"
    CD34pos_step1 = ag.gateTiltedLine(my_sample, name="tilted_gate_cd34", xCol='CD34 PE-Cy7-A', yCol=marker_cd45,
                                      startPoint=(xlim_middle0, maxVal), endLimits=(None, maxVal + 7500), theta=40,
                                      scale='bilog', T=200, population='lower', parentGate=CD16neg_56neg_final,
                                      filePlot=file_name)
    linnegCD34pos = ag.gatePC(my_sample, name="Lin-CD34+", xCol='CD34 PE-Cy7-A', yCol=marker_cd45, center='centroid',
                              adjustAngle=3, widthScale=2.5, heightScale=3.5, parentGate=CD34pos_step1)

    my_sample.update(ag.AGgate(linnegCD34pos, CD45pos, 'CD34 PE-Cy7-A', marker_cd45, "Lin-CD34+_1"), MFI=True,
                     extra_MFI=antibody_list, QC=True, scale='bilog', xscale='bilog', yscale='bilog', T=200,
                     xlim=[-1000, 80000], ylim=[0, 150000])
    # Out of parent
    my_sample.update(ag.AGgate(linnegCD34pos, CD16neg_56neg_final, 'CD34 PE-Cy7-A', marker_cd45, "Lin-CD34+_2"),
                     MFI=False, QC=False)
    # Ratio of Lin-CD34+ over total CD34+
    my_sample.update(ag.AGgate(linnegCD34pos, CD34pos, 'CD34 PE-Cy7-A', marker_cd45, "Lin-CD34+_3"),
                     MFI=False, QC=False)

    ####################################################################################################################
    # Calculate all the relevant thresholds, and produce images as a sanity check
    ####################################################################################################################

    # CD45RA - Using all lineage negative CD34+ events
    if f'fcs/{date_plate}/{sample_name}' in CD45RA_BLOB_list:  # Use fixed gate for the weird blob samples
        CD45RA_threshold = 200
    elif f'fcs/{date_plate}/{sample_name}' in CD45RA_threshold_100_list:
        CD45RA_threshold = 100
    elif f'fcs/{date_plate}/{sample_name}' in CD45RA_threshold_120_list:
        CD45RA_threshold = 120
    elif f'fcs/{date_plate}/{sample_name}' in CD45RA_threshold_200_list:
        CD45RA_threshold = 200
    elif f'fcs/{date_plate}/{sample_name}' in CD45RA_threshold_280_list:
        CD45RA_threshold = 280
    else:
        CD45RA_threshold = ag.densityDelimitation(my_sample, xCol='CD45RA FITC-A', interval=[50, 200],
                                                  limit_threshold=0.1, direction='right', scale='bilog', T=200,
                                                  parentGate=linnegCD34pos)

    if CD45RA_threshold == inf:
        CD45RA_threshold = 200
    # Save image for QC
    if file_name:  # Just a control to check that the CD90 limit is being calculated correctly. Gate produces
        # image but gate output is not stored
        file_name = image_path_prefix + "CD45RA_threshold/" + date_plate + "-" + sample_name + "-CD45RA_CONTROL.jpeg"
    ag.gateThreshold(my_sample, name="CD45RA_threshold", xCol='CD45RA FITC-A', yCol='CD90 PE (R-phycoerythrin)-A',
                     scale='bilog', T=200, thresh=CD45RA_threshold, orientation='vertical', population='lower',
                     parentGate=linnegCD34pos, filePlot=file_name)

    ####################################################################################################################
    # CD90 -  Using total CD45+ events to get the CD90 limit right (take only the right side of the CD45 cells)
    CD45pos_right_side = ag.gateThreshold(my_sample, name="right_blob", xCol='CD45RA FITC-A',
                                          yCol='CD90 PE (R-phycoerythrin)-A', scale='bilog', T=200, thresh=6500,
                                          parentGate=CD45pos, orientation='vertical', population='upper')

    if f'fcs/{date_plate}/{sample_name}' in CD90_threshold_1000_list:
        CD90_threshold = 1000
    elif f'fcs/{date_plate}/{sample_name}' in CD90_threshold_2000_list:
        CD90_threshold = 2000
    elif f'fcs/{date_plate}/{sample_name}' in CD90_threshold_5000_list:
        CD90_threshold = 5000
    else:
        CD90_threshold = ag.densityDelimitation(my_sample, xCol='CD90 PE (R-phycoerythrin)-A', interval=[200, 5000],
                                                limit_threshold=0.04, direction='right', scale='bilog', T=200,
                                                parentGate=CD45pos_right_side)
    if CD90_threshold == inf:  # Plan B: fixed gate
        CD90_threshold = 1500

    # Save image for QC
    if file_name:  # Just a control to check that the CD90 limit is being calculated correctly. Gate produces
        # image but gate output is not stored
        file_name = image_path_prefix + "CD90_threshold/" + date_plate + "-" + sample_name + "-CD90_CONTROL.jpeg"
    ag.gateThreshold(my_sample, name="CD90_threshold", xCol='CD45RA FITC-A',
                     yCol='CD90 PE (R-phycoerythrin)-A', scale='bilog', T=200,
                     thresh=CD90_threshold, orientation='horizontal', population='lower',
                     parentGate=CD45pos, filePlot=file_name)

    ####################################################################################################################
    # CD10 - Using 'CD16neg_56neg_final'
    CD10_threshold = ag.densityDelimitation(my_sample, xCol='CD10 APC (Allophycocyanin)-A',
                                            interval=[-100, 800], limit_threshold=0.2, direction='right',
                                            scale='bilog', T=200, parentGate=CD16neg_56neg_final)

    # Save image for QC
    if file_name:  # Just a control to check that the CD90 limit is being calculated correctly. Gate produces
        # image but gate output is not stored
        file_name = image_path_prefix + "CD10_threshold/" + date_plate + "-" + sample_name + "-CD10_CONTROL.jpeg"
    ag.gateThreshold(my_sample, name="CD10_threshold", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A',
                     scale='bilog', T=200, thresh=CD10_threshold, orientation='horizontal', population='lower',
                     parentGate=CD16neg_56neg_final, filePlot=file_name)

    ####################################################################################################################
    # CD135 - Using all CD45 positive events
    CD135_threshold = ag.densityDelimitation(my_sample, xCol='CD135 (BV) 711-A', interval=[-100, 800],
                                             limit_threshold=0.7, direction='left', scale='bilog', T=200,
                                             parentGate=CD16neg_56neg_final)

    # Save image for QC
    if file_name:  # Just a control to check that the CD90 limit is being calculated correctly. Gate produces
        # image but gate output is not stored
        file_name = image_path_prefix + "CD135_threshold/" + date_plate + "-" + sample_name + "-CD135_CONTROL.jpeg"
    ag.gateThreshold(my_sample, name="CD135_threshold", xCol='CD45RA FITC-A', yCol='CD135 (BV) 711-A',
                     scale='bilog', T=200, thresh=CD135_threshold, orientation='horizontal', population='lower',
                     parentGate=CD45pos, filePlot=file_name)

    ####################################################################################################################

    def gate_beyond_cd38(my_sample, CD38_threshold, CD45RA_threshold, CD90_threshold, CD10_threshold, CD135_threshold,
                         file_name=None):
        """
        This function implements the remaining gates, 12 to 17. It is a separate function so that it can be looped with
        different CD38 threshold values, as the more stringent the threshold is, the more enriched the population is
        supposed to be for HSCs.

        :param my_sample: AliGater sample object
        :param CD38_threshold: float showing the percentage of events considered CD38 negative
        :param CD45RA_threshold: CD45RA threshold, calculated with densityDelimitation
        :param CD90_threshold: CD90 threshold, calculated with densityDelimitation
        :param CD10_threshold: CD10 threshold, calculated with densityDelimitation
        :param CD135_threshold: CD135 threshold, calculated with densityDelimitation
        :param file_name: file name for the image to be saved. Default is None, in which case no image is saved.
        :return:
        """
        subfolder_name = f'CD38-{CD38_threshold}'
        name_suffix = f'_CD38_{CD38_threshold}'

        # TODO: Add minimum # of events requirement at the beginning of each loop?
        ################################################################################################################
        # Gate 12: CD38 out of Lin-CD34+
        ################################################################################################################

        gated_df = my_sample.fcsDF.loc[linnegCD34pos()].copy()
        cd38_values_list = sorted(list(gated_df['CD38 (BV) 421-A']))
        if cd38_values_list:  # If the list is not empty
            index = int(len(cd38_values_list) * CD38_threshold)  # Find index corresponding to event at threshold
            ylim_middle = cd38_values_list[index]
        else:
            print(f"WARNING: CD38 value list is empty for sample {sample_name}. Skipping.")
            return my_sample

        # Gate 12a: CD38-
        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/12-CD38neg/{date_plate}-{sample_name}-CD38neg{name_suffix}.jpeg'
        CD38neg = ag.gateThreshold(my_sample, name="CD38-", xCol='CD34 PE-Cy7-A', yCol='CD38 (BV) 421-A', scale='bilog',
                                   T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horizontal',
                                   population='lower', filePlot=file_name)
        # Gate 12b: CD38+
        CD38pos = ag.gateThreshold(my_sample, name="CD38+", xCol='CD34 PE-Cy7-A', yCol='CD38 (BV) 421-A', scale='bilog',
                                   T=200, thresh=ylim_middle, parentGate=linnegCD34pos, orientation='horizontal',
                                   population='upper')

        ################################################################################################################
        # Gate 13b: CD90+ out of CD38+
        CD90_onCD38pos = ag.gateCorner(my_sample, name="CD90posCD38pos", xCol='CD45RA FITC-A',
                                       yCol='CD90 PE (R-phycoerythrin)-A', scale='bilog', T=200,
                                       xThresh=CD45RA_threshold, yThresh=CD90_threshold,
                                       xOrientation='lower', yOrientation='upper', Outer=False,
                                       parentGate=CD38pos)

        my_sample.update(ag.AGgate(CD90_onCD38pos, CD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',
                                   f'CD90+CD38+{name_suffix}_1'), MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(ag.AGgate(CD90_onCD38pos, linnegCD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',
                                   f'CD90+CD38+{name_suffix}_2'), MFI=False, QC=False)

        ################################################################################################################
        # Gate 13: HSCs, MPPs and MLPs out of CD38-
        ################################################################################################################

        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/13-HSC_MLP_MPP/{date_plate}-{sample_name}-HSC_MLP_MPP{name_suffix}.jpeg'
        HSC, CD90posCD45RApos, MLP, MPP = ag.quadGate(my_sample, names=['HSC', 'NA', 'MLP', 'MPP'],
                                                      xCol='CD45RA FITC-A',
                                                      yCol='CD90 PE (R-phycoerythrin)-A', xThresh=CD45RA_threshold,
                                                      yThresh=CD90_threshold, scale='bilog', T=200, parentGate=CD38neg,
                                                      filePlot=file_name)

        my_sample.update(ag.AGgate(HSC, CD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'HSC{name_suffix}_1'),
                         MFI=True, extra_MFI=antibody_list, QC=True, scale='bilog', xscale='bilog', yscale='bilog',
                         T=200, xlim=[0, 1000], ylim=[300, 11000])
        my_sample.update(
            ag.AGgate(HSC, linnegCD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'HSC{name_suffix}_2'),
            MFI=False, QC=False)
        my_sample.update(ag.AGgate(MLP, CD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'MLP{name_suffix}_1'),
                         MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(
            ag.AGgate(MLP, linnegCD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'MLP{name_suffix}_2'),
            MFI=False, QC=False)
        my_sample.update(ag.AGgate(MPP, CD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'MPP{name_suffix}_1'),
                         MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(
            ag.AGgate(MPP, linnegCD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A', f'MPP{name_suffix}_2'),
            MFI=False, QC=False)
        my_sample.update(ag.AGgate(CD90posCD45RApos, CD34pos, 'CD45RA FITC-A', 'CD90 PE (R-phycoerythrin)-A',
                                   f'CD38-CD90+CD45RA+{name_suffix}'), MFI=True, extra_MFI=antibody_list, QC=False)

        ################################################################################################################
        # Gate 14: non-BNK (CD10neg) and B-NK out of CD38+
        ################################################################################################################

        # CD10+, B-NK progenitors
        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/14-B_NK/{date_plate}-{sample_name}-B_NK{name_suffix}.jpeg'
        BNK_prog = ag.gateCorner(my_sample, name="BNK", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A',
                                 xThresh=CD45RA_threshold, yThresh=CD10_threshold,
                                 xOrientation='upper', yOrientation='upper',
                                 Outer=False, scale='bilog', T=200, parentGate=CD38pos, filePlot=file_name)

        my_sample.update(
            ag.AGgate(BNK_prog, CD38pos, 'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A', f'B-NK{name_suffix}_1'),
            MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(
            ag.AGgate(BNK_prog, CD34pos, 'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A', f'B-NK{name_suffix}_2'),
            MFI=False, QC=False)
        my_sample.update(
            ag.AGgate(BNK_prog, linnegCD34pos, 'CD45RA FITC-A', 'CD10 APC (Allophycocyanin)-A', f'B-NK{name_suffix}_3'),
            MFI=False, QC=False)

        # CD10-, later divided into CMPs, GMPs and MEPs
        nonBNK = ag.gateThreshold(my_sample, name="CD10-", xCol='CD45RA FITC-A', yCol='CD10 APC (Allophycocyanin)-A',
                                  scale='bilog', T=200, thresh=CD10_threshold, orientation='horizontal',
                                  population='lower', parentGate=CD38pos)

        ################################################################################################################
        # Gate 15: CD10 out of MLP
        ################################################################################################################

        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/15-CD10pos(MLP)/{date_plate}-{sample_name}-CD10pos(MLP){name_suffix}.jpeg'
        CD10pos_MLP = ag.gateThreshold(my_sample, name="cd10+_MLP", xCol='CD45RA FITC-A',
                                       yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=CD10_threshold,
                                       parentGate=MLP, orientation='horizontal', population='upper', filePlot=file_name)
        CD10neg_MLP = ag.gateThreshold(my_sample, name="cd10-_MLP", xCol='CD45RA FITC-A',
                                       yCol='CD10 APC (Allophycocyanin)-A', scale='bilog', T=200, thresh=CD10_threshold,
                                       parentGate=MLP, orientation='horizontal', population='lower')

        my_sample.update(ag.AGgate(CD10pos_MLP, MLP, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10+(MLP){name_suffix}_1',
                                   IgnoreCellLimit=True),
                         MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(
            ag.AGgate(CD10pos_MLP, CD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10+(MLP){name_suffix}_2',
                      IgnoreCellLimit=True),
            MFI=False, QC=False)
        my_sample.update(
            ag.AGgate(CD10pos_MLP, linnegCD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10+(MLP){name_suffix}_3',
                      IgnoreCellLimit=True),
            MFI=False, QC=False)
        my_sample.update(ag.AGgate(CD10neg_MLP, MLP, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10-(MLP){name_suffix}_1',
                                   IgnoreCellLimit=True),
                         MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(
            ag.AGgate(CD10neg_MLP, CD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10-(MLP){name_suffix}_2',
                      IgnoreCellLimit=True),
            MFI=False, QC=False)
        my_sample.update(
            ag.AGgate(CD10neg_MLP, linnegCD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CD10-(MLP){name_suffix}_3',
                      IgnoreCellLimit=True),
            MFI=False, QC=False)

        ################################################################################################################
        # Gate 16: CMPs, GMPs and MEPs out of non-BNK (CD10-) CD38+ cells
        ################################################################################################################

        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/16-CMP_GMP_MEP/{date_plate}-{sample_name}-CMP_GMP_MEP{name_suffix}.jpeg'
        CMP, GMP, irrelevant, MEP = ag.quadGate(my_sample, names=['CMP', 'GMP', 'NA', 'MEP'],
                                                xCol='CD45RA FITC-A', yCol='CD135 (BV) 711-A',
                                                xThresh=CD45RA_threshold, yThresh=CD135_threshold,
                                                scale='bilog', T=200, parentGate=nonBNK, filePlot=file_name)

        my_sample.update(ag.AGgate(CMP, CD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CMP{name_suffix}_1'), MFI=True,
                         extra_MFI=antibody_list, QC=False)
        my_sample.update(ag.AGgate(CMP, linnegCD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CMP{name_suffix}_2'),
                         MFI=False,
                         QC=False)
        my_sample.update(ag.AGgate(CMP, nonBNK, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'CMP{name_suffix}_3'), MFI=False,
                         QC=False)
        my_sample.update(ag.AGgate(MEP, CD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'MEP{name_suffix}_1'), MFI=True,
                         extra_MFI=antibody_list, QC=False)
        my_sample.update(ag.AGgate(MEP, linnegCD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'MEP{name_suffix}_2'),
                         MFI=False,
                         QC=False)
        my_sample.update(ag.AGgate(MEP, nonBNK, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'MEP{name_suffix}_3'), MFI=False,
                         QC=False)
        my_sample.update(ag.AGgate(GMP, CD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'GMP{name_suffix}_1'), MFI=True,
                         extra_MFI=antibody_list, QC=False)
        my_sample.update(ag.AGgate(GMP, linnegCD34pos, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'GMP{name_suffix}_2'),
                         MFI=False,
                         QC=False)
        my_sample.update(ag.AGgate(GMP, nonBNK, 'CD45RA FITC-A', 'CD135 (BV) 711-A', f'GMP{name_suffix}_3'), MFI=False,
                         QC=False)
        ################################################################################################################
        # Gate 17: CD10- and CD135- from "HSC"s
        ################################################################################################################

        if file_name:
            file_name = f'{image_path_prefix}{subfolder_name}/17-trueHSC/{date_plate}-{sample_name}-trueHSC{name_suffix}.jpeg'
        true_HSC = ag.gateCorner(my_sample, name="trueHSC", xCol='CD10 APC (Allophycocyanin)-A',
                                 yCol='CD135 (BV) 711-A',
                                 xThresh=CD10_threshold, yThresh=CD135_threshold,
                                 scale='bilog', T=200, xOrientation='lower',
                                 yOrientation='lower', Outer=False, parentGate=HSC, filePlot=file_name)

        my_sample.update(ag.AGgate(true_HSC, CD34pos, 'CD10 APC (Allophycocyanin)-A', 'CD135 (BV) 711-A',
                                   f'trueHSC{name_suffix}_1'), MFI=True, extra_MFI=antibody_list, QC=False)
        my_sample.update(ag.AGgate(true_HSC, linnegCD34pos, 'CD10 APC (Allophycocyanin)-A', 'CD135 (BV) 711-A',
                                   f'trueHSC{name_suffix}_2'), MFI=False, QC=False)
        my_sample.update(ag.AGgate(true_HSC, HSC, 'CD10 APC (Allophycocyanin)-A', 'CD135 (BV) 711-A',
                                   f'trueHSC{name_suffix}_3'), MFI=False, QC=False)

        return my_sample

    # Loop over CD38 threshold values and store all the outputs into 'my_sample'
    for CD38_threshold in [0.3, 0.2, 0.1, 0.05]:
        my_sample = gate_beyond_cd38(my_sample, CD38_threshold, CD45RA_threshold, CD90_threshold, CD10_threshold,
                                     CD135_threshold, file_name)

    return my_sample


if __name__ == '__main__':

    # Folder structure to save images at. Maybe should be a function
    # Note that 'image_path_prefix' is a global variable and is used by 'gate_full_dataset'
    date_string = str(datetime.date.today())  # 2021-10-25 , for example
    week_string = "w" + str(datetime.date.today().isocalendar()[1])  # w43 , for example
    week_and_year_string = week_string + "_" + str(datetime.date.today().isocalendar()[0])  # 'w43_2021' , for example

    week_and_year_string = 'w33_2022_FINAL_2'  # TODO: delete this

    image_path_prefix = "../output/gating/aligater_images/images_" + week_and_year_string + "/"
    if not os.path.exists(image_path_prefix):
        os.mkdir(image_path_prefix)  # Create parent folder for all images at this date

    out_folders = get_out_folder_list(image_path_prefix)  # List containing desired folder structure to save images in
    # if not os.path.exists(
    #        image_path_prefix + "10-CD16CD56/"):  # create subfolder so sub-sub folders can be created next step
    #    os.mkdir(image_path_prefix + "10-CD16CD56/")
    for individual_folder in out_folders:  # Actual creation of the folders
        if not os.path.exists(individual_folder):
            os.makedirs(individual_folder)

    batch_size = 500
    lowest_sampleID = 1
    highest_sampleID = lowest_sampleID + batch_size - 1
    maximum_date = datetime.date(2022, 6, 24)
    test_bool = False  # TODO: find a cleaner way to do this

    while lowest_sampleID < 3900:
        print("Gating all samples from CBID" + str(lowest_sampleID) + " to CBID" + str(highest_sampleID))
        filepaths = get_filepaths(max_date=maximum_date, min_sample_id=lowest_sampleID, max_sample_id=highest_sampleID,
                                  include_repeats=False)  # Get list of paths to each fcs file of interest

        # TESTING
        if test_bool:
            print("RUNNING ALIGATER IN TEST MODE")
            filepaths = filepaths[2:5]  # Take only 3 files per batch for testing
        # Define experiment object
        exp_folder_name = "cord_blood_experiment-" + week_and_year_string + "-samples_" + str(
            lowest_sampleID) + "_to_" + str(highest_sampleID)
        CB_exp = ag.AGExperiment(filepaths, filters=['fcs'], mask=['30min', '45min', 'Neg', 'test'] + get_blacklist(),
                                 experiment_name=exp_folder_name, flourochrome_area_filter=True, QC=True,
                                 QCbins=128)  # Set up experiment (including folder for output)
        CB_exp.apply(gate_full_dataset)  # Apply this function to every sample

        CB_exp.printExperiment(
            "../output/gating/aligater_output/" + exp_folder_name + "/cblood_phenotypes-" + week_and_year_string
            + "-samples_" + str(lowest_sampleID) + "_to_" + str(highest_sampleID) + ".txt")  # save main output file

        # Update limits to measure next batch
        lowest_sampleID = highest_sampleID + 1
        highest_sampleID = highest_sampleID + batch_size
