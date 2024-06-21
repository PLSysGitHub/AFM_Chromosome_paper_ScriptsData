import pandas as pd
pd.options.display.max_columns = 999
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from scipy import stats

'''Script for generating Figure of the manuscript.'''

# ---- Representative AFM curve ----
curve_path = "Representative_AFM-Rheo_Curve.txt"
df_curve = pd.read_csv(curve_path, sep='\t')
# cantilever spring constant
cantilever_k = 0.137e12
df_curve["Force (pN)"] = df_curve["Vertical Deflection"].astype(float) * cantilever_k

print("Curve header:", df_curve.columns)

# Results from the microrheology analysis from JPK-Bruker software
df_path = "Rheology_results_allChrms.txt"
df = pd.read_csv(df_path, sep='\t')
# Calculating mean and s.e.m of G' (Storage) and G'' (Loss)
G_Loss_avg = df.groupby("Oscillation Frequency [Hz]")["Shear Loss G'' [Pa]"].mean().values
G_Loss_sem = df.groupby("Oscillation Frequency [Hz]")["Shear Loss G'' [Pa]"].sem().values
G_Storage_avg = df.groupby("Oscillation Frequency [Hz]")["Shear Storage G' [Pa]"].mean().values
G_Storage_sem = df.groupby("Oscillation Frequency [Hz]")["Shear Storage G' [Pa]"].sem().values
# Importing data of Loss Tangent from OT experiments (Meijering et al., Nature)
f_fast = np.loadtxt("OT_NaturePaper_data/f_fast.dat", unpack=True)
f_slow = np.loadtxt("OT_NaturePaper_data/f_slow.dat", unpack=True)
k_im_fast = np.loadtxt("OT_NaturePaper_data/k_im_fast.dat", unpack=True)
k_im = np.loadtxt("OT_NaturePaper_data/k_im.dat", unpack=True)
k_re_fast = np.loadtxt("OT_NaturePaper_data/k_re_fast.dat", unpack=True)
k_re = np.loadtxt("OT_NaturePaper_data/k_re.dat", unpack=True)
std_im_fast = np.loadtxt("OT_NaturePaper_data/std_im_fast.dat", unpack=True)
std_re_fast = np.loadtxt("OT_NaturePaper_data/std_re_fast.dat", unpack=True)
std_im = np.loadtxt("OT_NaturePaper_data/std_im.dat", unpack=True)
std_re = np.loadtxt("OT_NaturePaper_data/std_re.dat", unpack=True)
lt_freq = np.loadtxt("OT_NaturePaper_data/lt_freq.dat", unpack=True)
lt = np.loadtxt("OT_NaturePaper_data/lt.dat", unpack=True)
std_lt = np.loadtxt("OT_NaturePaper_data/std_lt.dat", unpack=True)
sem_lt = np.loadtxt("OT_NaturePaper_data/sem_lt.dat", unpack=True)

#  ---- Plotting ----
# single column width = 8,9 cm, 2-column width = 18,2 cm
cm = 1/2.54  # centimeters in inches
fig = plt.figure(figsize=(18*cm, 18*cm))
ax = fig.subplot_mosaic([['a', 'a', 'a'],
                         ['b', 'c', 'd'],
                         ['e', 'f', 'g']])
for label, a in ax.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(1/72, 3/72, fig.dpi_scale_trans)
    a.text(0.0, 1.0, label, transform=a.transAxes + trans,
            fontsize='large', weight='bold', va='bottom', fontfamily='sans-serif')
# Representative Force-Time curve
ax['a'].plot(df_curve["Series Time"], df_curve["Force (pN)"])
ax['a'].set_xlabel("Time (s)")
ax['a'].set_ylabel("Force (pN)")
ax['a'].set_xlim(0, 2.7)
# add labels for the different parts of the curve
sections = df_curve.groupby("Section")["Series Time"].max().to_numpy()
ax['a'].vlines(x=sections[:-1], ymin=-20, ymax=170, colors='gray', ls='--')
ax['a'].text(.4, 150, r'Approach', fontsize=10)
ax['a'].text(0.95, 150, r'Relaxation', fontsize=10)
ax['a'].text(1.45, 150, r'Oscillation', fontsize=10)
ax['a'].text(2.05, 150, r'Retraction', fontsize=10)

freq = [2, 20, 200]
pal_b = sns.color_palette("Blues", 20)
pal_c = sns.color_palette("Reds", 20)

# Plotting G' and G'' values in errorbars for all the chromosomes
for i, chrm in enumerate(df["Filename"].unique()):
    avg_g1 = df[df["Filename"] == chrm].groupby("Oscillation Frequency [Hz]")["Shear Storage G' [Pa]"].mean().to_numpy()
    sem_g1 = df[df["Filename"] == chrm].groupby("Oscillation Frequency [Hz]")["Shear Storage G' [Pa]"].sem().to_numpy()
    ax['b'].errorbar(x=freq, y=avg_g1, yerr=sem_g1, label="G'", c=pal_b[i],
                     marker="o",markersize=3, linestyle='-', linewidth=1,
                     capsize=2, capthick=1)
    avg_g2 = df[df["Filename"] == chrm].groupby("Oscillation Frequency [Hz]")["Shear Loss G'' [Pa]"].mean().to_numpy()
    sem_g2 = df[df["Filename"] == chrm].groupby("Oscillation Frequency [Hz]")["Shear Loss G'' [Pa]"].sem().to_numpy()
    ax['c'].errorbar(x=freq, y=avg_g2, yerr=sem_g2, label="G''", c=pal_c[i],
                     marker="o", markersize=3, linestyle='-', linewidth=1,
                     capsize=2, capthick=1)

ax['b'].set_xlabel("Oscillation Frequency (Hz)")
ax['b'].set_ylabel("G' (Pa)")
ax['b'].set_xscale('log')
ax['b'].set_yscale('log')
ax['c'].set_xlabel("Oscillation Frequency (Hz)")
ax['c'].set_ylabel("G'' (Pa)")
ax['c'].set_xscale('log')
ax['c'].set_yscale('log')

# Average G' and G''
ax['d'].errorbar(x=freq, y=G_Storage_avg, yerr=G_Storage_sem,
                 label="G'", c=pal_b[15], marker="o",
                 markersize=3, linestyle='-', linewidth=1,
                 capsize=2, capthick=1)
ax['d'].errorbar(x=freq, y=G_Loss_avg, yerr=G_Loss_sem,
                 label="G''", c=pal_c[15], marker="o",
                 markersize=3, linestyle='-', linewidth=1,
                 capsize=2, capthick=1)
exp_linex = np.linspace(2, 200, 100)
exp_liney = 60*exp_linex ** 0.5
ax['d'].plot(exp_linex, exp_liney, c='grey', linestyle='--')
ax['d'].set_xscale('log')
ax['d'].set_yscale('log')
ax['d'].set_xlabel("Oscillation Frequency (Hz)")
ax['d'].set_ylabel("$G'_{avg}$ , $G''_{avg}$ (Pa)")
ax['d'].legend(loc="lower right", frameon=False)

# Loss Tangent for AFM
lt_avg = df.groupby("Oscillation Frequency [Hz]")["Loss Tan"].mean().values
lt_sem = df.groupby("Oscillation Frequency [Hz]")["Loss Tan"].sem().values
ax['e'].errorbar(x=freq, y=lt_avg, yerr=lt_sem, c="k", label="AFM",
                 marker="o", markersize=3, linestyle='-',
                 linewidth=1, capsize=2, capthick=1)
ax['e'].errorbar(x=lt_freq, y=lt, yerr=sem_lt, c="gray", label="OT",
                 marker="o", markersize=3, linestyle='-',
                 linewidth=1, capsize=2, capthick=1)
ax['e'].set_xscale('log')
ax['e'].set_yscale('log')
ax['e'].set_xlabel("Oscillation Frequency (Hz)")
ax['e'].set_ylabel("Loss Tangent")
ax['e'].legend(loc="center", frameon=False)

# Loading the average force-relaxation curves for the chromosomes + fits
avg_curve_df = pd.read_csv("Average_ForceRelax_curves.txt", sep='\t')
df_fit = pd.read_csv("Fit_ForceRelax_curves.txt", sep='\t')
# Plotting
pal_f = sns.color_palette("colorblind", 20)
for i, chrm in enumerate(avg_curve_df["Chromosome"].unique()):
    time = avg_curve_df[avg_curve_df["Chromosome"] == chrm]["Time (s)"].to_numpy()
    force = avg_curve_df[avg_curve_df["Chromosome"] == chrm]["Force (pN)"].to_numpy()
    fit = df_fit[df_fit["Chromosome"] == chrm]["fit"].to_numpy()
    ax['f'].plot(time, force / np.max(force), c=pal_f[i], alpha=.5)
    ax['f'].plot(time, fit / np.max(force), '--', c=pal_f[i], linewidth=.7)
ax['f'].set_xlabel("Time (s)")
ax['f'].set_ylabel("$F(t=0) / F(t\\rightarrow\\infty)$")
ax['f'].set_xlim(0, 0.5)

# Loading the fitting parameters
fit_params_df = pd.read_csv("Fitting_parameters.txt", sep='\t')
t1 = fit_params_df["t1"].to_numpy()
t1_std = fit_params_df["t1_std"].to_numpy()
t2 = fit_params_df["t2"].to_numpy()
print("Fs: ", fit_params_df["Delta Force (pN)"].mean(), fit_params_df["Delta Force (pN)"].sem())
t2_std = fit_params_df["t2_std"].to_numpy()
print(f"Average $\\tau_1$: {np.round(np.mean(t1), decimals=3)} ± {np.round(stats.sem(t1), decimals=3)}")
print(f"Average $\\tau_2$: {np.round(np.mean(t2), decimals=3)} ± {np.round(stats.sem(t2), decimals=3)}")
DeltaForce = []
F0 = []
# Plotting characteristic relaxation times vs Force drops
for i, chrm in enumerate(avg_curve_df["Chromosome"].unique()):
    force = avg_curve_df[avg_curve_df["Chromosome"] == chrm]["Force (pN)"].to_numpy()
    F0.append(force[0])
    DeltaForce.append(force[0] - fit_params_df["Delta Force (pN)"].to_numpy()[i])
fit_params_df["F0 (pN)"] = F0
fit_params_df["Force_drop (pN)"] = DeltaForce
ax['g'].errorbar(x=DeltaForce, y=t2, yerr=t2_std, marker="o",
                 markersize=3, linestyle='', linewidth=1,
                 capsize=2, capthick=1, label="$\\tau_2$", c='skyblue')
ax['g'].errorbar(x=DeltaForce, y=t1, yerr=t1_std, marker="o",
                 markersize=3, linestyle='', linewidth=1,
                 capsize=2, capthick=1, label="$\\tau_1$", c='teal')
ax['g'].set_yscale("log")
ax['g'].set_ylabel("Time (s)")
ax['g'].set_xlabel("$\\Delta$Force (pN)")
ax['g'].legend(frameon=False)

# Saving the results of the fitting parameters to a .csv for Table in SI
# fit_params_df.to_csv("TableResults_from_fit_and_figure.csv", sep=',')

plt.tight_layout()
# saving the figure
# plt.savefig("RheologyAFM.svg", dpi=500)
plt.show()

