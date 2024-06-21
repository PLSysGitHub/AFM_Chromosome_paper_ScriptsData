# import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Nanoscope_converter import nanoscope_converter
from nanoscope_ChromosomeCurveFit import indentation_fit
import glob  # this is for listing all the files in a given directory
# import cv2
import os
from os.path import isfile
import math
# import seaborn as sns
from scipy.stats import sem


# # function for calculating the distance between two points
def dist(x1, y1, x2, y2):
    return np.sqrt((x2-x1)**2 + (y2-y1)**2)


def fit_func_linear(x, a, b):  # sensitivity is corrected by applying a linear fit to the last part of the curve
    return a * x + b


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


def YMAcrossAxis_function(path):
    curves_path = sorted([f for f in glob.glob(path + '**/*', recursive=True)])
    # ---- if the total number of curves is not a multiple of n,
    # curves are not ordered in lines of n --> manually check! ----
    n = 10  # number of curves for each line
    if len(curves_path) % n != 0:
        print("Wrong number of curves in the folder")
        exit()
    # If everything's good --> process the curves in chunks of length equal to n
    curves_name = []
    all_curves = []
    cpoints = []
    x_coords_tot = []
    y_coords_tot = []
    total_distances = []
    total_YM_values = []
    total_YM_area_values = []
    total_cp_x_values = []
    total_R2 = []
    ALL_dist = []
    ALL_YM = []
    ALL_YM_area = []
    ALL_cp_x = []
    ALL_R2 = []
    for i in range(0, len(curves_path), n):
        sub_list = curves_path[i: i+n]
        x_coords = []
        y_coords = []
        YM_values = []
        YM_area_values = []
        cp_x_values = []
        R2_values = []
        valid = []
        for path in sub_list:
            curve_data = indentation_fit(path)
            x_pos, y_pos, YM, R2, cp_x, YM_area = curve_data[0], curve_data[1],\
                                                  curve_data[2], curve_data[3],\
                                                  curve_data[4], curve_data[6]
            # Selecting only those curves with contact point higher than 200 or 100 nm
            # (the rest is regarded as surface)
            if cp_x > 150:
                x_coords.append(x_pos)
                y_coords.append(y_pos)
                YM_values.append(YM)
                R2_values.append(R2)
                YM_area_values.append(YM_area)
                cp_x_values.append(cp_x)
            else:
                continue
        #  if the chromocome is laying horizontal, rotate all the lines and
        #  process them as if all the chromosomes were vertical
            # x_coords.append(x_pos)
            # y_coords.append(y_pos)
            # YM_values.append(YM)
            # R2_values.append(R2)
            # YM_area_values.append(YM_area)
            # cp_x_values.append(cp_x)
        # get the coordinates of the middle point of the line
        xm = (max(x_coords) + min(x_coords)) / 2
        ym = (max(y_coords) + min(y_coords)) / 2
        '''try this new approach (and delete the if statement above)'''
        if np.abs((x_coords[0] - x_coords[-1])) < np.abs((y_coords[0] - y_coords[-1])):
            xm, ym = ym, xm
            x_coords, y_coords = y_coords, x_coords
        # computing the distances of each curve from the middle
        distances = []
        for j in range(len(x_coords)):
            if xm >= x_coords[j]:
                distance = - dist(x_coords[j], y_coords[j], xm, ym)
            else:
                distance = dist(x_coords[j], y_coords[j], xm, ym)
            # distance = dist(x_coords[j], y_coords[j], xm, ym)
            distances.append(distance)
        total_distances.extend(distances)
        total_YM_values.extend(YM_values)
        total_R2.extend(R2_values)
        total_YM_area_values.extend(YM_area_values)
        total_cp_x_values.extend(cp_x_values)
        # plt.scatter(distances, eta_values)
        # plt.show()
        # plt.scatter(x_coords, y_coords)
        # plt.scatter(xm, ym, marker='x')
    # plt.show()
    ALL_dist.extend(total_distances)
    ALL_YM.extend(total_YM_values)
    ALL_R2.extend(total_R2)
    ALL_YM_area.extend(total_YM_area_values)
    ALL_cp_x.extend(total_cp_x_values)
    # plt.show()
    # fig, ax = plt.subplots()
    # ax.set_xlim(-700, 700)
    # ax.set_ylim(0, 0.15)
    # ax.scatter(ALL_dist, ALL_eta, s=4)
    # plt.show()
    return ALL_dist, ALL_YM, ALL_R2, ALL_cp_x, ALL_YM_area