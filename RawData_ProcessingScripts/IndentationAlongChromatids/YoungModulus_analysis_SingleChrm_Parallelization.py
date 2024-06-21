# import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Nanoscope_converter import nanoscope_converter
from nanoscope_ChromosomeCurveFit import area_viscoelasticity
import glob  # this is for listing all the files in a given directory
# import cv2
import os
from os.path import isfile
import math
# import seaborn as sns
from scipy.stats import sem
from YoungModulus_acrossAxis_function import YMAcrossAxis_function
import sys


# ---- Custom function for binning all the data points ----
def lin_bin(x, y, n_bins):
    boundaries = np.linspace(min(x), max(x), n_bins+1)
    x_m = []
    y_m = []
    x_er = []
    y_er = []
    for i in range(n_bins):
        sel = np.array((x > boundaries[i]) & (x <= boundaries[i+1]))
        x_m.append(np.average([x[i] for i in range(len(sel)) if sel[i]]))
        y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]]))
        x_er.append(sem([x[i] for i in range(len(sel)) if sel[i]]))
        y_er.append(sem([y[i] for i in range(len(sel)) if sel[i]]))
    return np.array(x_m), np.array(y_m), np.array(x_er), np.array(y_er)

# ---- parsing the input parameter as a script variable ----
j = int(sys.argv[1])
# j = 0
print(j)

# if the script is running on the office computer,
# then change the folder path to the local folder, otherwise use the USB
# if os.getenv('COMPUTERNAME') == 'DESKTOP-EGCQ9GC':
#     all_chrms_path = r"C:\Users\Andrea Ridolfi\Desktop\VU material\Lab\20221219_Chrms_Indentation_different_chrms_areas\Chromosome_indentations"
# else:
#     all_chrms_path = r"E:\20221219_Chrms_Indentation_different_chrms_areas\Chromosome_indentations"

all_chrms_path = "Chromosome_indentations"

print('Computing...')
single_chrms_folders = [f for f in glob.glob(all_chrms_path + '**/*', recursive=True)]
# Debugging
# single_chrms_folders = [f for f in glob.glob(all_chrms_path + '**/*', recursive=True) if "PS130520" in f]
# print(single_chrms_folders)
# end Debugging
# path to the folder containing the indentation files
telomere_dists = []
centromere_dists = []
arms_dists = []
telomere_YM = []
telomere_YM_area = []
centromere_YM = []
centromere_YM_area = []
arms_YM = []
arms_YM_area = []
telomere_cp_x = []
centromere_cp_x = []
arms_cp_x = []
telomere_R2 = []
centromere_R2 = []
arms_R2 = []

folder = single_chrms_folders[j]
# folder = single_chrms_folders
chrm_name = folder[-8:]
print(chrm_name)
# storing the directories of the indentation folders
indent_folders = [f for f in glob.glob(single_chrms_folders[j] + '**/*', recursive=True)
                  if not(isfile(f))]
# for each indentation area, process the curves with the viscoelastic index function
for indent_area_folder in indent_folders:
    # print(indent_area_folder)
    distances, YM_vals, R2_vals, cp_x_vals, YM_area_vals = YMAcrossAxis_function(indent_area_folder)
    if 'Telomere' in indent_area_folder:
        telomere_dists.extend(distances)
        telomere_YM.extend(YM_vals)
        telomere_R2.extend(R2_vals)
        telomere_YM_area.extend(YM_area_vals)
        telomere_cp_x.extend(cp_x_vals)
    if 'Centromere' in indent_area_folder:
        centromere_dists.extend(distances)
        centromere_YM.extend(YM_vals)
        centromere_R2.extend(R2_vals)
        centromere_YM_area.extend(YM_area_vals)
        centromere_cp_x.extend(cp_x_vals)
    if 'Arms' in indent_area_folder:
        arms_dists.extend(distances)
        arms_YM.extend(YM_vals)
        arms_R2.extend(R2_vals)
        arms_YM_area.extend(YM_area_vals)
        arms_cp_x.extend(cp_x_vals)

df_telomere = pd.DataFrame({'dist': telomere_dists,
                           'YM': telomere_YM,
                           'R2': telomere_R2,
                           'YM_area': telomere_YM_area,
                           'area': ['telomere']*len(telomere_dists),
                           'cpx': telomere_cp_x,
                           'chrm': chrm_name})
df_centromere = pd.DataFrame({'dist': centromere_dists,
                           'YM': centromere_YM,
                           'R2': centromere_R2,
                           'YM_area': centromere_YM_area,
                           'area': ['centromere']*len(centromere_dists),
                           'cpx': centromere_cp_x,
                           'chrm': chrm_name})
df_arms = pd.DataFrame({'dist': arms_dists,
                           'YM': arms_YM,
                           'R2': arms_R2,
                           'YM_area': arms_YM_area,
                           'area': ['arms']*len(arms_dists),
                           'cpx': arms_cp_x,
                           'chrm': chrm_name})

ALL_df = pd.concat([df_telomere, df_centromere, df_arms])

# if the script is running on the office computer,
# then change the folder path to the local folder, otherwise use the USB
# if os.getenv('COMPUTERNAME') == 'DESKTOP-EGCQ9GC':
#     dst_path = r"C:\Users\Andrea Ridolfi\Desktop\YM_Indentation_parts_values_Contact150nm.csv"
# else:
#     dst_path = r"C:\Users\Andrea\Desktop\YM_Indentation_parts_values_Contact150nm.csv"

dst_path = "YM_Indentation_parts_values_Contact150nm.txt"
# print(ALL_df)
ALL_df.to_csv("20240424//" + chrm_name + "_" + dst_path, sep='\t')
print('Finished...')


#
# import matplotlib.pyplot as plt
# import pandas as pd
# import numpy as np
# import os
# import seaborn as sns
# from scipy.stats import sem
#
# # ---- Custom function for binning all the data points ----
# def lin_bin(x, y, n_bins, R2):
#     boundaries = np.linspace(min(x), max(x), n_bins+1)
#     x_m = []
#     y_m = []
#
#     x_er = []
#     y_er = []
#     print('R2 from fit: ', np.average(R2))
#     for i in range(n_bins):
#         sel = np.array((x > boundaries[i]) & (x <= boundaries[i+1]))
#         # R2_val = [R2[i] for i in range(len(sel)) if sel[i]]
#         x_m.append(np.average([x[i] for i in range(len(sel)) if sel[i]]))
#         # y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]], weights=R2_val))
#         y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]]))
#         x_er.append(sem([x[i] for i in range(len(sel)) if sel[i]]))
#         y_er.append(sem([y[i] for i in range(len(sel)) if sel[i]]))
#     return np.array(x_m), np.array(y_m), np.array(x_er), np.array(y_er)
#
# # if the script is running on the office computer,
# # then change the folder path to the local folder, otherwise use the USB
# if os.getenv('COMPUTERNAME') == 'DESKTOP-EGCQ9GC':
#     path = r"C:\Users\Andrea Ridolfi\Desktop\YM_Indentation_parts_values_Contact150nm.csv"
# else:
#     path = r"C:\Users\Andrea\Desktop\YM_Indentation_parts_values_Contact150nm.csv"
#
# df = pd.read_csv(path)
# df_telomere = df[df['area'] == 'telomere']
# df_centromere = df[df['area'] == 'centromere']
# df_arms = df[df['area'] == 'arms']
#
# # df_telomere = df_telomere[df_telomere['cpx'] > 150]
# # df_centromere = df_centromere[df_centromere['cpx'] > 150]
# # df_arms = df_arms[df_arms['cpx'] > 150]
# df_telomere = df_telomere[df_telomere['YM'] > 0]
# df_centromere = df_centromere[df_centromere['YM'] > 0]
# df_arms = df_arms[df_arms['YM'] > 0]
# df_telomere['YM_area'] = df_telomere['YM_area'] * 1e3
# df_centromere['YM_area'] = df_centromere['YM_area'] * 1e3
# df_arms['YM_area'] = df_arms['YM_area'] * 1e3
# df_telomere = df_telomere[df_telomere['R2'] >= 0]
# df_centromere = df_centromere[df_centromere['R2'] >= 0]
# df_arms = df_arms[df_arms['R2'] >= 0]
# print(len(df_telomere['YM']))
# print(len(df_centromere['YM']))
# print(len(df_arms['YM']))
#
# # ---- Plotting ----
# n_bins = 20
# fig, ax = plt.subplots(1, 3)
# c_palette = sns.color_palette('rainbow')
# palette = 'rainbow'
#
# # ---- TELOMERES ----
# ax[0].set_title('Telomeres')
# ax[0].set_xlabel('Distance from centre (nm)')
# ax[0].set_ylabel('YM (kPa)')
# ax[0].set_ylim(0, 2000)
# sns.scatterplot(x=df_telomere['dist'], y=df_telomere['YM'], alpha=0.15, ax=ax[0], c=c_palette[1])
# tel_dist_bin, tel_eta_bin, tel_dist_err, tel_eta_err = \
#     lin_bin(np.array(df_telomere['dist']), np.array(df_telomere['YM']), n_bins, np.array(df_telomere['R2']))
# ax[0].errorbar(tel_dist_bin, tel_eta_bin, xerr=tel_dist_err, yerr=tel_eta_err,
#                fmt='o', c=c_palette[1])
#
# # check for thickness distribution
# # ax0 = ax[0].twinx()
# # ax0.set_ylim(150, 800)
# # ax0.scatter(df_telomere['dist'], df_telomere['cpx'], alpha=0.25, c='peachpuff')
# # tel_dist_bin, tel_eta_bin, tel_dist_err, tel_eta_err = \
# #     lin_bin(np.array(df_telomere['dist']), np.array(df_telomere['cpx']), n_bins, np.array(df_telomere['R2']))
# # ax0.errorbar(tel_dist_bin, tel_eta_bin, xerr=tel_dist_err, yerr=tel_eta_err,
# #                fmt='o', c='orange')
#
# print('Telomere YM (mean, std) in the central part (kPa): ',
#       tel_eta_bin[8:-6].mean(), tel_eta_bin[8:-6].std())
#
# # ---- CENTROMERE ----
# ax[1].set_title('Centromere')
# ax[1].set_xlabel('Distance from centre (nm)')
# ax[1].set_ylabel('YM (kPa)')
# ax[1].set_ylim(0, 2000)
# sns.scatterplot(x=df_centromere['dist'], y=df_centromere['YM'], alpha=0.15, ax=ax[1], c=c_palette[3])
# cent_dist_bin, cent_eta_bin, cent_dist_err, cent_eta_err =\
#     lin_bin(np.array(df_centromere['dist']), np.array(df_centromere['YM']), n_bins, np.array(df_centromere['R2']))
# ax[1].errorbar(cent_dist_bin, cent_eta_bin, xerr=cent_dist_err, yerr=cent_eta_err,
#                fmt='o', c=c_palette[3])
#
# # ----check for thickness distribution ----
# # ax1 = ax[1].twinx()
# # ax1.set_ylim(150, 800)
# # ax1.scatter(df_centromere['dist'], df_centromere['cpx'], alpha=0.25, c='peachpuff')
# # cent_dist_bin, cent_eta_bin, cent_dist_err, cent_eta_err = \
# #     lin_bin(np.array(df_centromere['dist']), np.array(df_centromere['cpx']), n_bins, np.array(df_centromere['R2']))
# # ax1.errorbar(cent_dist_bin, cent_eta_bin, xerr=cent_dist_err, yerr=cent_eta_err,
# #                fmt='o', c='orange')
#
# print('Centromere YM (mean, std) in the central part (kPa): ',
#       cent_eta_bin[8:-6].mean(), cent_eta_bin[8:-6].std())
#
# # ---- ARMS ----
# ax[2].set_title('Arms')
# ax[2].set_xlabel('Distance from centre (nm)')
# ax[2].set_ylabel('YM (kPa)')
# ax[2].set_ylim(0, 2000)
# sns.scatterplot(x=df_arms['dist'], y=df_arms['YM'], alpha=0.15, ax=ax[2], c=c_palette[5])
# arm_dist_bin, arm_eta_bin, arm_dist_err, arm_eta_err =\
#     lin_bin(np.array(df_arms['dist']), np.array(df_arms['YM']), n_bins, np.array(df_arms['R2']))
# ax[2].errorbar(arm_dist_bin, arm_eta_bin, xerr=arm_dist_err, yerr=arm_eta_err,
#                fmt='o', c=c_palette[5])
#
# # check for thickness distribution
# # ax2 = ax[2].twinx()
# # ax2.set_ylim(150, 800)
# # ax2.scatter(df_arms['dist'], df_arms['cpx'], alpha=0.25, c='peachpuff')
# # arm_dist_bin, arm_eta_bin, arm_dist_err, arm_eta_err = \
# #     lin_bin(np.array(df_arms['dist']), np.array(df_arms['cpx']), n_bins,  np.array(df_arms['R2']))
# # ax2.errorbar(arm_dist_bin, arm_eta_bin, xerr=arm_dist_err, yerr=arm_eta_err,
# #                fmt='o', c='orange')
# #
# print('Arms YM (mean, std) in the central part (kPa): ',
#       arm_eta_bin[8:-6].mean(), arm_eta_bin[8:-6].std())
#
# plt.show()
#
#
# fig, ax = plt.subplots()
# # ax.set_ylim(0, 1000)
# ax.set_ylabel('Young Modulus (kPa)')
# ax.set_xlabel('Distance from centre (nm)')
# ax.plot(tel_dist_bin, tel_eta_bin, color=c_palette[1], linestyle='solid', marker='o', label='Telomeres')
# ax.errorbar(tel_dist_bin, tel_eta_bin, xerr=tel_dist_err, yerr=tel_eta_err,
#                fmt='o', c=c_palette[1])
# ax.plot(cent_dist_bin, cent_eta_bin, color=c_palette[3], linestyle='solid', marker='o', label='Centromere')
# ax.errorbar(cent_dist_bin, cent_eta_bin, xerr=cent_dist_err, yerr=cent_eta_err,
#                fmt='o', c=c_palette[3])
# ax.plot(arm_dist_bin, arm_eta_bin, color=c_palette[5], linestyle='solid', marker='o', label='Arms')
# ax.errorbar(arm_dist_bin, arm_eta_bin, xerr=arm_dist_err, yerr=arm_eta_err,
#                fmt='o', c=c_palette[5])
# plt.legend()
# plt.show()
#
