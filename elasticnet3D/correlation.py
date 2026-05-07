import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
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

            pear_corr_altitude = summary['corr_altitude'][i,j]
            pear_corr_azimuth = summary['corr_azimuth'][i,j]
            # pear_corr_pa = summary['corr_pa'][i,j]

            # pear_corr_altitude = np.corrcoef(retinotopy[mask_var_idx,1], result[:,1])[0,1]
            # pear_corr_azimuth = np.corrcoef(retinotopy[mask_var_idx,0], result[:,0])[0,1]

            #determine spearman correlation values
            spear_corr_altitude = spearmanr(retinotopy[mask_var_idx,1], result[:,1])[0]
            spear_corr_azimuth = spearmanr(retinotopy[mask_var_idx,0], result[:,0])[0]
            #can't do circular spearman correlation for PA, so just report the pearson correlation for PA

            #now for polar coordinates
            #simulation data
            retinotopy_ecc, retinotopy_pa = dst.cartesian_to_polar(retinotopy[:,0],retinotopy[:,1])
            retinotopy_pol = np.column_stack((retinotopy_ecc,retinotopy_pa)) #[ecc, pa] in [deg]
            #empirical data
            result_ecc, result_pa = dst.cartesian_to_polar(result[:,0],result[:,1])
            corr_pa = spearmanr(np.pi/180*retinotopy_pol[mask_var_idx, 1], np.pi/180*result_pa)
            corr_ecc = spearmanr(retinotopy_pol[mask_var_idx, 0], result_ecc)

            if n == 25:
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
                # Save the correlation values to a CSV file
                os.makedirs('./results', exist_ok=True)
                csv_file_path = f'./results/corr_{x}_b1_{b1}_b2_{b2}.csv'
                with open(csv_file_path, mode='w', newline='') as csv_file:
                    writer = csv.writer(csv_file)
                    writer.writerow(['Correlation Type', 'Value'])
                    writer.writerow(['Pearson Correlation Azimuth', pear_corr_azimuth])
                    writer.writerow(['Pearson Correlation Altitude', pear_corr_altitude])
                    # writer.writerow(['Circular Correlation Polar Angle', pear_corr_pa])
                    writer.writerow(['Spearman Correlation Azimuth', spear_corr_azimuth])
                    writer.writerow(['Spearman Correlation Altitude', spear_corr_altitude])
                    writer.writerow(['Spearman Correlation Polar Angle', corr_pa.correlation])
                    writer.writerow(['Spearman Correlation Eccentricity', corr_ecc.correlation])

                print(f"Saved correlation values to {csv_file_path}")
