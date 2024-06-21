import pandas as pd
pd.options.display.max_columns = 999
import numpy as np
import matplotlib.pyplot as plt
import glob

''' Plotting and saving a representative Force-Time curve for the 
    AFM-rheological measurements'''

repr_curve_folder = "Representative_ForceTime_curves"
curves_path = sorted([c for c in glob.glob(repr_curve_folder + '**/*', recursive=True)], reverse=True)
curve_path_2 = curves_path[0]
curve_path_20 = curves_path[1]
curve_path_200 = curves_path[2]

all_curve_paths = [curve_path_2, curve_path_20, curve_path_200]
# add to the two following list the frequencies and paths of curves you want displayed (2, 20, 200 Hz)
osc_freq = ["20 Hz"]
curve_paths = [curve_path_20]
# Plotting
fig, ax = plt.subplot_mosaic([['a']])
for i, curve_path in enumerate(curve_paths):
    # Stripping the path in order to get the chromosome number and its respective curve number
    chrm_name = curve_path.split('/')[-1].split('-')[0]
    curve_index = curve_path.split('/')[-1].split('_')[-1][:2]
    with open(curve_path, 'r') as doc:
        curve = doc.read()
        #  Each section of the experiment is divided by the expression '\n#\n'
        curve = curve.split("\n#\n")
        # print(len(curve))
        # Extracting the useful parts of the curve data
        # Approach phase
        extend_data = curve[2].split('\n\n')[0].split("\n")
        # Relaxation phase
        pause_data = curve[4].split('\n\n')[0].split("\n")
        # Oscillation phase
        modulation_data = curve[6].split('\n\n')[0].split("\n")
        # Retraction phase
        retraction_data = curve[8].split('\n\n')[0].split("\n")
    # Header for the Approach (and Retraction) phase
    extend_header = ["Vertical Tip Position", "Vertical Deflection",
                       "Height", "Fast Scanner Height (measured & smoothed)",
                       "Fast Scanner Height (measured)", "Height (measured & smoothed)",
                       "Height (measured)", "Series Time", "Segment Time"]
    # Header for the Relaxation (and Oscillation) phase
    pause_header = ["Vertical Tip Position", "Vertical Deflection",
                 "Height", "Height (measured & smoothed)",
                 "Height (measured)", "Series Time", "Segment Time"]
    # Load the parts of the curve as .txt and store them as DataFrames
    extend_data = np.loadtxt(extend_data, delimiter=' ')
    pause_data = np.loadtxt(pause_data, delimiter=' ')
    modulation_data = np.loadtxt(modulation_data, delimiter=' ')
    retraction_data = np.loadtxt(retraction_data, delimiter=' ')

    df_extend = pd.DataFrame(extend_data)
    df_extend.columns = extend_header
    # additional column for identification and later plots
    df_extend["Section"] = "Approach"
    df_retraction = pd.DataFrame(retraction_data)
    df_retraction.columns = extend_header
    df_retraction["Section"] = "Retraction"
    # saving a subset with only the meaningful info
    df_extend = df_extend[["Vertical Tip Position", "Vertical Deflection",
                 "Height", "Height (measured & smoothed)",
                 "Height (measured)", "Series Time", "Segment Time", "Section"]]
    df_retraction = df_retraction[["Vertical Tip Position", "Vertical Deflection",
                 "Height", "Height (measured & smoothed)",
                 "Height (measured)", "Series Time", "Segment Time", "Section"]]

    df_pause = pd.DataFrame(pause_data)
    df_pause.columns = pause_header
    df_pause["Section"] = "Relaxation"
    df_modulation = pd.DataFrame(modulation_data)
    df_modulation.columns = pause_header
    df_modulation["Section"] = "Oscillation"
    # Concatenating all the different parts of the curve in a single DataFrame
    curve_df = pd.concat([df_extend, df_pause, df_modulation, df_retraction])
    # Saving the DataFrame as a .txt file
    # curve_df.to_csv("Representative_AFM-Rheo_Curve.txt", sep='\t')

    # PLotting (Vertical Deflection is multiplied by the cantilever spring constant)
    ax['a'].plot(curve_df["Series Time"].astype(float), curve_df["Vertical Deflection"].astype(float) * 0.137e12,
                 label=osc_freq[i])
ax['a'].set_xlabel("Time (s)")
ax['a'].set_ylabel("Force (pN)")
ax['a'].legend()
plt.title("AFM microrheology measurement")
plt.show()
