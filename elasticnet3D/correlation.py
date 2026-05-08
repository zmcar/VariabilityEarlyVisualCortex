# add correlation values for a subject to one .csv file. 5x5 matrix
# find greatest correlation value for that subject, add to matrix with 'subject_id', 'b1', 'b2', 'pearson/spearman', 'azimuth/altitude/polar angle/eccentricity/circular correlation', 'correlation val'
#find max correlation val, plot scatter for this simulation vs empirical
# could return count for all maximum correlation values, determine which b1/b2 is most common

# check for patient id in directory
# if none, check for patient_id (done)
# if none, return "no 'results' folder"

# check 2d vs flat

# want to compare all correlation scores, finding the highest score for each patient, add this to a list, and then find the highest score of all the patients.
# preserve the b1 and b2 values to determine if there is a trend for a specific b1, or b2, or b1 + b2 combination
# add condition that checks for the patient id in the directory, if there is no matching folder, check for patient_id (done)
# additional condition to return "no 'result' folder"
# compute the correlation values for 2d vs flat as well

import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import pandas as pd
import csv
import argparse
import numpy as np
import functions.dstools as dst
from scipy.stats import spearmanr
from pathlib import Path
import scipy.io as sio
import matplotlib.pyplot as plt
os.chdir(r'c:\Users\61439\Variability_Early_Visual_Cortex_new\VariabilityEarlyVisualCortex') 
loadDir = r'\\storage.erc.monash.edu\shares\MNHS-dshi0006\VariabilityEarlyVisualCortex';

#iterating through each subject 

#add in a check or skip if the .mat file doesn't exist

#iterating through each combination of b1 and b2 values to load corresponding .mat files, extract correlation values, and save them to CSV files
#use commented code for all subjects once folder structure confirmed
# list_file = os.path.join(os.path.dirname(__file__), 'list_subj.txt')
# with open('list_subj.txt', 'r', encoding='utf-8') as f:
#     subject_id = [line.strip() for line in f if line.strip()]
subject_id = ['100610', '114823', '192439']

b1_b2_vals = [10, 20, 40, 80, 160]

for x in subject_id:
    n==0
    thisDir = os.path.join(loadDir, x+'(done)')

    #load original retinotopy data for the subject
    arealBorders = sio.loadmat(os.path.join(thisDir, 'arealBorder_' + x + '.mat'))
    maltitude = arealBorders['grid_altitude_i'] 
    mazimuth = arealBorders['grid_azimuth_i']
    shape2d = arealBorders['areaMatrix'][0][0].shape
    retinotopy = np.column_stack((mazimuth.ravel(order='F'), maltitude.ravel(order='F')))
    
    #check against main code
    nAreas = len(arealBorders['areaMatrix'][0])
    mask_idx = [None] * nAreas
    for iarea in range(nAreas):
        mask_idx[iarea] = np.where(arealBorders['areaMatrix'][0][iarea].ravel(order='F') == 1)[0]
    varIdx = [1, 2]
    mask_var_idx = np.concatenate([mask_idx[index] for index in varIdx])

    summary = sio.loadmat(os.path.join(thisDir, f'summary_correlation_{x}.mat'))

    # set up matrices
    # need to add corresponding b1/b2 values as titles
    spear_corr_altitude = np.zeros(5,5)
    spear_corr_azimuth = np.zeros(5,5)
    corr_pa = np.zeros(5,5)
    corr_ecc = np.zeros(5,5)

    for i in range(len(b1_b2_vals)):
        for j in range(len(b1_b2_vals)):
            n =+1
            b1 = b1_b2_vals[i]
            b2 = b1_b2_vals[j]
            
            print(f"Processing {x}, b1={b1}, b2={b2}")
            # Load the .mat file
            mat_file_path = os.path.join(thisDir, f'summary_{x}V+D_{x}_b1_{b1}_b2_{b2}.mat')
            if not os.path.exists(mat_file_path):
                print(f"MAT file not found for b1={b1}, b2={b2}, skipping.")
                continue
            mat_contents = sio.loadmat(mat_file_path)

            # Extract the results and pearson correlation coefficient values of the simulation
            result = mat_contents['result'][0,0]
            result = np.asarray(result)
            result = np.squeeze(result)
            
            # result = mat_contents['result'][0]
            # result = np.squeeze(result)

            # result2d = mat_contents['result2d']
            # result2d_flat = mat_contents['result2d_flat']

            # will create six 5x5 matrices for each patient_id

            #pearson correlation (cartesian)
            pear_corr_altitude = summary['corr_altitude'][i,j]
            pear_corr_azimuth = summary['corr_azimuth'][i,j]

            # do we need this circular corelation?
            # pear_corr_pa = summary['corr_pa'][i,j]

            #spearman correlation (cartesian)
            spear_corr_altitude[i,j] = spearmanr(retinotopy[mask_var_idx,1], result[:,1])[0]
            spear_corr_azimuth[i,j] = spearmanr(retinotopy[mask_var_idx,0], result[:,0])[0]
            #can't do circular spearman correlation for PA, so just report the pearson correlation for PA

            #convert data to polar
            #simulation data
            retinotopy_ecc, retinotopy_pa = dst.cartesian_to_polar(retinotopy[:,0],retinotopy[:,1])
            retinotopy_pol = np.column_stack((retinotopy_ecc,retinotopy_pa)) #[ecc, pa] in [deg]
            #empirical data
            result_ecc, result_pa = dst.cartesian_to_polar(result[:,0],result[:,1])

            #spearman correlation (polar)
            corr_pa[i,j] = spearmanr(np.pi/180*retinotopy_pol[mask_var_idx, 1], np.pi/180*result_pa)
            corr_ecc[i,j] = spearmanr(retinotopy_pol[mask_var_idx, 0], result_ecc)

            # on final iteration, we want to create several .csv files that tells us all correlation values, then determine maximum correlation value, and plot that
            if i==4 and j==4:
                # Save the correlation values to a CSV file
                os.makedirs('./results', exist_ok=True)
                row_col_labels = ["10", "20", "40", "80", "160"]
                df.index.name = "B2/B1"
                # spearman correlation (altitude)
                csv_file_path = f'./results/spearcorr_{x}_altitude.csv'
                df = pd.DataFrame(spear_corr_altitude, index=row_col_labels, columns=row_col_labels)
                df.to_csv(csv_file_path)
                print(f"Saved correlation values to {csv_file_path}")

                # spearman correlation (azimuth)
                csv_file_path = f'./results/spearcorr_{x}_azimuth.csv'
                df = pd.DataFrame(spear_corr_azimuth, index=row_col_labels, columns=row_col_labels)
                df.to_csv(csv_file_path)
                print(f"Saved correlation values to {csv_file_path}")

                # spearman correlation (polar angle)
                csv_file_path = f'./results/spearcorr_{x}_pa.csv'
                df = pd.DataFrame(corr_pa, index=row_col_labels, columns=row_col_labels)
                df.to_csv(csv_file_path)
                print(f"Saved correlation values to {csv_file_path}")

                # spearman correlation (eccentricity)
                csv_file_path = f'./results/spearcorr_{x}_ecc.csv'
                df = pd.DataFrame(corr_ecc, index=row_col_labels, columns=row_col_labels)
                df.to_csv(csv_file_path)
                print(f"Saved correlation values to {csv_file_path}")

                #determine max correlation value for all b1/b2 combinations, then plot this one.

                
                
                #plotting empirical vs simulated data for visual representation of the correlation
                plt.figure(figsize=(12, 5))
                #eccentricity
                plt.subplot(1, 2, 1)
                plt.scatter(retinotopy[mask_var_idx,0], result[:,0], alpha=0.5)
                plt.xlabel('Empirical Azimuth')
                plt.ylabel('Simulated Azimuth')
                plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')

                #polar angle
                plt.subplot(1, 2, 2)
                plt.scatter(retinotopy[mask_var_idx,1], result[:,1], alpha=0.5)
                plt.xlabel('Empirical altitude')
                plt.ylabel('Simulated altitude')
                plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')
                plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')

                # #plotting empirical vs simulated data for visual representation of the correlation
                # plt.figure(figsize=(12, 5))
                # #eccentricity
                # plt.subplot(1, 2, 1)
                # plt.scatter(retinotopy_pol[mask_var_idx,0], result_ecc, alpha=0.5)
                # plt.xlabel('Empirical Eccentricity (deg)')
                # plt.ylabel('Simulated Eccentricity (deg)')
                # plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')

                # #polar angle
                # plt.subplot(1, 2, 2)
                # plt.scatter(retinotopy_pol[mask_var_idx,1], result_pa, alpha=0.5)
                # plt.xlabel('Empirical Polar Angle (deg)')
                # plt.ylabel('Simulated Polar Angle (deg)')
                # plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')
                # plt.title(f'Subject {x} - b1: {b1}, b2: {b2}')

                plt.savefig(f'./results/scatter_{x}_b1_{b1}_b2_{b2}.png', dpi=150)
                plt.close()
                
