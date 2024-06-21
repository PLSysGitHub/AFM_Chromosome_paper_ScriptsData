import pandas as pd
import numpy as np
from Nanoscope_converter import nanoscope_converter
from nanoscope_ChromosomeCurveFit import indentation_fit
import glob  # this is for listing all the files in a given directory
# import cv2
import os
import math
# import seaborn as sns
from scipy.stats import sem
import sys
# import matplotlib.pyplot as plt


# # function for calculating the distance between two points
def dist(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)


def fit_func_linear(x, a, b):  # sensitivity is corrected by applying a linear fit to the last part of the curve
   return a * x + b

# parameter to run the analysis within a for loop in the terminal
j = int(sys.argv[1])



all_path = "AFM_CurvesAcrossChromatids_AllData10"
folders = [f for f in glob.glob(all_path + '**/*', recursive=True)]
ALL_ym = []
ALL_err = []
ALL_dist = []
ALL_area = []
ALL_E_bis = []
ALL_cp_x = []

folder = folders[j]
chrm_name = folder[-8:]

# storing the directories of all the files
all_files = [f for f in glob.glob(folders[j] + '**/*', recursive=True)]
# image of the chromosome is the only .tif file in the folder
img_path = [item for item in all_files if '.tif' in item]
#chrm_img = cv2.imread(img_path[0])
# plt.imshow(chrm_img, cmap='afmhot')
# storing the path of all the curves (the filetype is .000)
curves_path = sorted([item for item in all_files if '.000' in item])
# ---- if the total number of curves is not a multiple of n,
n = 10  # number of curves for each line
if len(curves_path) % n != 0:
    print("Wrong number of curves in the folder")
    exit()
# If everything's good --> process the curves in chunks of length equal to n
curves_name = []
Young_modulus = []
Young_modulus_TipC = []
all_curves = []
cpoints = []
x_coords_tot = []
y_coords_tot = []
total_ym_values = []
total_area = []
total_distances = []
total_ym_err = []
total_E_bis = []
total_cp_x_values = []

for i in range(0, len(curves_path), n):
    sub_list = curves_path[i: i+n]
    # print(sub_list)
    x_coords = []
    y_coords = []
    ym_values = []
    area_val = []
    distances = []
    err_values = []
    E_bis_values = []
    cp_x_values = []
    valid = []
    for path in sub_list:
        curve_data = indentation_fit(path)
        x_pos, y_pos, ymodulus, ymodulus_err, cp_x, ind_area, E_bis = \
            curve_data[0], curve_data[1], curve_data[2], curve_data[3], curve_data[4], curve_data[5], curve_data[6]
        # Selecting only those curves with contact point higher than 150 nm
        # (the rest is regarded as surface)
        if cp_x > 150:
            x_coords.append(x_pos)
            y_coords.append(y_pos)
            area_val.append(ind_area)
            ym_values.append(ymodulus)
            err_values.append(ymodulus_err)
            cp_x_values.append(cp_x)
            E_bis_values.append(E_bis*1000)
        else:
            continue
    # get the coordinates of the middle point of the line
    xm = (x_coords[0] + x_coords[-1]) / 2
    ym = (y_coords[0] + y_coords[-1]) / 2
    #  if the chromosome is laying horizontal, rotate all the lines and
    #  process them as if all the chromosomes were vertical
    if np.abs((x_coords[0] - x_coords[-1])) < np.abs((y_coords[0] - y_coords[-1])):
        xm, ym = ym, xm
        x_coords, y_coords = y_coords, x_coords
    # computing the distances of each curve from the middle
    for j in range(len(x_coords)):
        if xm >= x_coords[j]:
            distance = - dist(x_coords[j], y_coords[j], xm, ym)
        else:
            distance = dist(x_coords[j], y_coords[j], xm, ym)
        # distance = dist(x_coords[j], y_coords[j], xm, ym)
        distances.append(distance)
    total_ym_values.extend(ym_values)
    total_distances.extend(distances)
    total_ym_err.extend(err_values)
    total_area.extend(area_val)
    total_E_bis.extend(E_bis_values)
    total_cp_x_values.extend(cp_x_values)
    # plt.scatter(distances, ym_values)
    # plt.show()
    # plt.scatter(x_coords, y_coords)
    # plt.scatter(xm, ym, marker='x')
ALL_dist.extend(total_distances)
ALL_ym.extend(total_ym_values)
ALL_err.extend(total_ym_err)
ALL_area.extend(total_area)
ALL_E_bis.extend(total_E_bis)
ALL_cp_x.extend(total_cp_x_values)
# plt.show()
# fig, ax = plt.subplots()
# ax.set_xlim(-700, 700)
# ax.set_ylim(0, 0.15)
# ax.errorbar(total_distances, total_ym_values, yerr=total_ym_err, fmt='o')
# plt.scatter(total_distances, total_ym_values)
# plt.show()

#  ---- Saving the obtained data as a dataframe ----
df = pd.DataFrame(
    {'ALL_ym': ALL_ym,
     'ALL_err': ALL_err,
     'ALL_dist': ALL_dist,
     'ALL_Area': ALL_area,
     'ALL_E_bis': ALL_E_bis,
     'ALL_cp_x': ALL_cp_x,
     'chrm': chrm_name})

df.to_csv("parallelization_results_20240424//" + chrm_name + "_YM_AcrossChromatids_AllData10.txt", sep='\t')
