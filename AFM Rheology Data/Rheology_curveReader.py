import pandas as pd
pd.options.display.max_columns = 999
import numpy as np
import matplotlib.pyplot as plt
import glob
from scipy.optimize import curve_fit


'''Extract the relaxation part from each force curve performed on a chromosome, compute the average force relaxation 
    curve for each chromosome and performs a double exponential fit to extract characteristic relaxation times.'''

path_2 = "Probed_chromosomes/Chromosome3_2Hz-data-2023.06.06-14.43.45.054_processed-2023.07.07-11.24.11/curves"
path_20 = "Probed_chromosomes/Chromosome3_20Hz-data-2023.06.06-14.41.55.666_processed-2023.07.07-11.24.12/curves"
path_200 = "Probed_chromosomes/Chromosome3_200Hz-data-2023.06.06-14.47.21.668_processed-2023.07.07-11.24.12/curves"

# Double exponential function for fitting the curves
def double_expFit(time, Af, tf, As, ti, Fs):
    F = Af * np.exp(-time / tf) + As * np.exp(-time / (ti)) + Fs
    return F

# Wrapper function to force ti to be always larger than tf
def wrapper_func(time, Af, tf, As, ti, Fs):
    if ti <= tf:
        return 1e18  # Return a very large value to discourage the fit
    return double_expFit(time, Af, tf, As, ti, Fs)

# Paths to the curves folders
curves_20 = sorted([c for c in glob.glob(path_20 + '**/*', recursive=True)])
curves_2 = sorted([c for c in glob.glob(path_2 + '**/*', recursive=True)])
curves_200 = sorted([c for c in glob.glob(path_200 + '**/*', recursive=True)])

def chrm_powerLawExponent(curves):
    """Takes as input the folder with the curves on a single chromosome, extract the Relaxation part, compute the
    average curve."""
    curves_list = sorted([c for c in glob.glob(curves + '**/*', recursive=True)])
    # DataFrame for storing data on the chromosome
    chrm_df = pd.DataFrame()
    for curve_path in curves_list[1:]:
        # print(curve_path)
        # Stripping the path in order to get the chromosome number and its respective curve number
        chrm_name = curve_path.split('/')[-1].split('-')[0]
        curve_index = curve_path.split('/')[-1].split('_')[-1][:2]
        with open(curve_path, 'r') as doc:
            curve = doc.read()
            #  Each section is divided by the pattern '\\r\\n#\\r\\n'
            curve = curve.split("\n#\n")
            # The fourth section contains the Relaxation part + some text
            # keeping the info and getting rid of the remaining text
            pause_data = curve[4].split('\n\n')[0].split("\n")
        # Storing pause_data into a Dataframe
        df_header = ["Vertical Tip Position", "Vertical Deflection",
                     "Height", "Height (measured & smoothed)",
                     "Height (measured)", "Series Time", "Segment Time"]
        pause_data = np.loadtxt(pause_data, delimiter=' ')
        df = pd.DataFrame(pause_data)
        df.columns = df_header
        df["Chromosome"] = chrm_name
        df["Curve_#"] = curve_index
        #  Excluding curves with contact points < 150nm (more representative of the surface than chromosome)
        if df["Height (measured & smoothed)"].astype(float)[0] * -1e9 > 150:
            # calculate the force (deflection * cantilever spring constant) for curves on chromosome
            df["Force (pN)"] = df["Vertical Deflection"].astype(float) * 0.137e12
            # store all the curves in a single DataFrame
            chrm_df = pd.concat([chrm_df, df], ignore_index=True)

    # Building the average force-relaxation curve for the chromosome
    avg_force = chrm_df.groupby("Segment Time", as_index=True)["Force (pN)"].mean().to_numpy()
    sem_force = chrm_df.groupby("Segment Time", as_index=True)["Force (pN)"].sem().to_numpy()
    delta_f = avg_force[0] - avg_force[-1]
    time = chrm_df["Segment Time"].unique()
    df_avgCurve = pd.DataFrame({
        "Time (s)": time,
        "Force (pN)": avg_force,
        "s.e.m": sem_force,
        "Delta F": delta_f,
        "Chromosome": chrm_name
    })
    return df_avgCurve

'''Iterating the function on all the chromosomes (only 20Hz oscillation exps have been considered) and applying the 
    fit'''
all_chrms = sorted([c for c in glob.glob("Probed_chromosomes" + '**/*', recursive=True) if "20Hz" in c])
all_avg_df = pd.DataFrame()
for chrm in all_chrms:
    path = str(chrm) + "//curves//"
    df_avg = chrm_powerLawExponent(path)
    # Storing all the average curves for all the chromsomes in one DataFrame
    all_avg_df = pd.concat([all_avg_df, df_avg], ignore_index=True)

fig, ax = plt.subplots(1, 2)
# Initializing arrays where the fitting parameters will be stored
fit_params = np.zeros((len(all_avg_df["Chromosome"].unique()), 5))
fit_param_std = np.zeros((len(all_avg_df["Chromosome"].unique()), 5))
DeltaForce = np.zeros(len(all_avg_df["Chromosome"].unique()))
# DataFrame to store the final results of the fit
df_fit = pd.DataFrame()
for i, chrm in enumerate(all_avg_df["Chromosome"].unique()):
    time = all_avg_df[all_avg_df["Chromosome"] == chrm]["Time (s)"].to_numpy()
    force = all_avg_df[all_avg_df["Chromosome"] == chrm]["Force (pN)"].to_numpy()
    sem = all_avg_df[all_avg_df["Chromosome"] == chrm]["s.e.m"].to_numpy()
    # Fitting the average curve
    [params, pcov] = curve_fit(wrapper_func, time, force, p0=[5, 0.05, 5, 0.05, 1],
                               bounds=([0, 0, 0, 0, 0], [np.inf, np.inf, np.inf, np.inf, np.inf]),
                               max_nfev=10000000)
    # standard deviation errors on the parameters
    p_err = np.sqrt(np.diag(pcov))
    fit = double_expFit(time, *params)
    # temporary dataframe for storing the current fit results
    df_t = pd.DataFrame({"fit": fit,
                        "Chromosome": chrm})
    # DataFrame containing all the fits for all the chromosomes
    df_fit = pd.concat([df_fit, df_t])
    fit_params[i] = params
    fit_param_std[i] = p_err
    DeltaForce[i] = all_avg_df[all_avg_df["Chromosome"] == chrm]["Delta F"].unique()
    # print(fit_params[i][1])

    # Plotting
    ax[0].plot(time, force)
    ax[0].plot(time, fit, c='r')
    ax[0].set_xlabel("Time (s)")
    ax[0].set_ylabel("Force (pN)")

# Saving the average curves and the fits in .txt files
# all_avg_df.to_csv("Average_ForceRelax_curves.txt", sep='\t')
# df_fit.to_csv("Fit_ForceRelax_curves.txt", sep='\t')

# Analysis of the fitting parameters
t1 = [fit_params[i][1] for i in range(fit_params.shape[0])]
t2 = [fit_params[i][3] for i in range(fit_params.shape[0])]
t1_std = [fit_param_std[i][1] for i in range(fit_param_std.shape[0])]
t2_std = [fit_param_std[i][3] for i in range(fit_param_std.shape[0])]
Fs = [fit_params[i][4] for i in range(fit_params.shape[0])]
fit_results = pd.DataFrame({
    "t1": t1,
    "t1_std": t1_std,
    "t2": t2,
    "t2_std": t2_std,
    "Delta Force (pN)": Fs})

# Saving the fitting parameters to a .txt file
# fit_results.to_csv("Fitting_parameters.txt", sep='\t')

print("Average t1 = " + str(np.mean(t1)) + " ± " + str(np.std(t1)))
print("Average t2 = " + str(np.mean(t2)) + " ± " + str(np.std(t2)))

# Plotting
ax[1].errorbar(x=DeltaForce, y=t1, yerr=t1_std, marker="o", linestyle='', label="$\\tau_1$")
ax[1].errorbar(x=DeltaForce, y=t2, yerr=t2_std, marker="o", linestyle='', label="$\\tau_2$")
ax[1].set_yscale("log")
ax[1].set_ylabel("Time (s)")
ax[1].set_xlabel("Delta F (pN)")
ax[1].legend()
plt.tight_layout()
plt.show()
