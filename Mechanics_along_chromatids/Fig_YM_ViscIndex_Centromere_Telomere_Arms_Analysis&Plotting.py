import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import pandas as pd
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
import numpy as np
import seaborn as sns
from scipy.stats import sem, ttest_ind

'''This script generates the plots displayed in the figures regarding chromosome mechanics along the chromatids'''

# ---- Young Modulus (Y.M.) data analysis ----

n_bins = 15
def lin_bin_YM(x, y, n_bins):
    '''Bins the data in x and y, i.e., distance from CCA and YM, and compute the avg and s.e.m for each bin.'''
    boundaries = np.linspace(min(x), max(x), n_bins+1)
    x_m = []
    y_m = []
    x_er = []
    y_er = []
    for i in range(n_bins):
        sel = np.array((x >= boundaries[i]) & (x <= boundaries[i+1]))
        x_m.append(np.average([x[i] for i in range(len(sel)) if sel[i]]))
        y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]]))
        x_er.append(sem([x[i] for i in range(len(sel)) if sel[i]]))
        y_er.append(sem([y[i] for i in range(len(sel)) if sel[i]]))
    return np.array(x_m), np.array(y_m), np.array(x_er), np.array(y_er)


# Path for the results of all the Young Modulus fits
path = "20240424_YM_Indentation_parts_values_Contact150nm.txt"
df = pd.read_csv(path, sep='\t')
# Separate the results of telomeres, centromeres and arms
df_telomere = df[df['area'] == 'telomere']
df_centromere = df[df['area'] == 'centromere']
df_arms = df[df['area'] == 'arms']
# Discard values of YM < 0 or coming from fits with R2 <= 0
df_telomere = df_telomere[df_telomere['YM'] > 0]
df_centromere = df_centromere[df_centromere['YM'] > 0]
df_arms = df_arms[df_arms['YM'] > 0]
df_telomere = df_telomere[df_telomere['R2'] > 0]
df_centromere = df_centromere[df_centromere['R2'] > 0]
df_arms = df_arms[df_arms['R2'] > 0]
# Convert the YM calculated from the area under the curve into kPa
df_telomere['YM_area'] = df_telomere['YM_area'] * 1e3
df_centromere['YM_area'] = df_centromere['YM_area'] * 1e3
df_arms['YM_area'] = df_arms['YM_area'] * 1e3

# Analysis of the results for Telomeres
tel_dist_bin_ym, tel_eta_bin_ym, tel_dist_err_ym, tel_eta_err_ym = \
    lin_bin_YM(np.array(df_telomere['dist']), np.array(df_telomere['YM']), n_bins)
print('Telomere Y.M. (mean, s.e.m) in the central part (kPa): ',
      tel_eta_bin_ym[4:11].mean(), tel_eta_bin_ym[4:-6].std() / len(tel_eta_bin_ym[4:-6]))

# Analysis of the results for Centromeres
cent_dist_bin_ym, cent_eta_bin_ym, cent_dist_err_ym, cent_eta_err_ym =\
    lin_bin_YM(np.array(df_centromere['dist']), np.array(df_centromere['YM']), n_bins)
print('Centromere Y.M. (mean, s.e.m) in the central part (kPa): ',
      cent_eta_bin_ym[4:11].mean(), cent_eta_bin_ym[4:11].std() / len(cent_eta_bin_ym[4:11]))

# Analysis of the results for Arms
arm_dist_bin_ym, arm_eta_bin_ym, arm_dist_err_ym, arm_eta_err_ym =\
    lin_bin_YM(np.array(df_arms['dist']), np.array(df_arms['YM']), n_bins)
print('Arms YM (mean, s.e.m) in the central part (kPa): ',
      arm_eta_bin_ym[4:11].mean(), arm_eta_bin_ym[4:11].std() / len(arm_eta_bin_ym[4:11]))




# ---- Viscoelasticity (eta) data analysis ----

def lin_bin(x, y, n_bins):
    '''Bins the data in x and y, i.e., distance from CCA and eta, and compute the avg and s.e.m for each bin.'''
    boundaries = np.linspace(min(x), max(x), n_bins+1)
    x_m = []
    y_m = []
    x_er = []
    y_er = []
    for i in range(n_bins):
        sel = np.array((x >= boundaries[i]) & (x <= boundaries[i+1]))
        x_m.append(np.average([x[i] for i in range(len(sel)) if sel[i]]))
        y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]]))
        x_er.append(sem([x[i] for i in range(len(sel)) if sel[i]]))
        y_er.append(sem([y[i] for i in range(len(sel)) if sel[i]]))
    return np.array(x_m), np.array(y_m), np.array(x_er), np.array(y_er)


# Path for the results of all the viscoelasticity values
path = "20240424_ViscIndex_Indentation_parts_values_Contact150nm.txt"
df_eta = pd.read_csv(path, sep='\t')
# Separate the results of telomeres, centromeres and arms
df_telomere_eta = df_eta[df_eta['area'] == 'telomere']
df_centromere_eta = df_eta[df_eta['area'] == 'centromere']
df_arms_eta = df_eta[df_eta['area'] == 'arms']
# Discarding values for which is NOT  0<= eta <= 1
df_telomere_eta = df_telomere_eta[df_telomere_eta['eta'] <= 1]
df_centromere_eta = df_centromere_eta[df_centromere_eta['eta'] <= 1]
df_arms_eta = df_arms_eta[df_arms_eta['eta'] <= 1]
df_telomere_eta = df_telomere_eta[df_telomere_eta['eta'] > 0]
df_centromere_eta = df_centromere_eta[df_centromere_eta['eta'] > 0]
df_arms_eta = df_arms_eta[df_arms_eta['eta'] > 0]

# Analysis of the viscoelasticity results for Telomeres
tel_dist_bin, tel_eta_bin, tel_dist_err, tel_eta_err = \
    lin_bin(np.array(df_telomere_eta['dist']), np.array(df_telomere_eta['eta']), n_bins)
# Analysis of the viscoelasticity results for Centromeres
cent_dist_bin, cent_eta_bin, cent_dist_err, cent_eta_err =\
    lin_bin(np.array(df_centromere_eta['dist']), np.array(df_centromere_eta['eta']), n_bins)
# Analysis of the viscoelasticity results for Arms
arm_dist_bin, arm_eta_bin, arm_dist_err, arm_eta_err =\
    lin_bin(np.array(df_arms_eta['dist']), np.array(df_arms_eta['eta']), n_bins)


# ---- Figure avg profiles ----
# single column widht = 8,9 cm, 2-column width = 18,2 cm
cm = 1/2.54  # centimeters in inches
fig = plt.figure(figsize=(18*cm, 9.5*cm))
ax = fig.subplot_mosaic([['a', 'b']])
# Annotations parameters
for label, a in ax.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(1/72, 3/72, fig.dpi_scale_trans)
    a.text(0.0, 1.0, label, transform=a.transAxes + trans,
            fontsize='large', weight='bold', va='bottom', fontfamily='sans-serif')
# Palettes used in the figures
c_palette = sns.color_palette('rainbow')
palette = [c_palette[1], c_palette[3], c_palette[5]]
# Joining the Y.M., eta results of telomeres, centromeres and arms into two dataframe
df_eta = pd.concat([df_telomere_eta, df_arms_eta, df_centromere_eta])
df_ym = pd.concat([df_telomere, df_arms, df_centromere])

# Young Modulus
ax['a'].set_xlabel('Distance from CCA (nm)')
ax['a'].set_ylabel('Young Modulus (kPa)')
ax['a'].plot(tel_dist_bin_ym, tel_eta_bin_ym, color=c_palette[1], linestyle='solid', marker='o', label='Peripheral regions')
ax['a'].errorbar(tel_dist_bin_ym, tel_eta_bin_ym, xerr=tel_dist_err_ym, yerr=tel_eta_err_ym,
               fmt='o', c=c_palette[1])
ax['a'].plot(cent_dist_bin_ym, cent_eta_bin_ym, color=c_palette[3], linestyle='solid', marker='o', label='Centromeric regions')
ax['a'].errorbar(cent_dist_bin_ym, cent_eta_bin_ym, xerr=cent_dist_err_ym, yerr=cent_eta_err_ym,
               fmt='o', c=c_palette[3])
ax['a'].plot(arm_dist_bin_ym, arm_eta_bin_ym, color=c_palette[5], linestyle='solid', marker='o', label='Mid-arm regions')
ax['a'].errorbar(arm_dist_bin_ym, arm_eta_bin_ym, xerr=arm_dist_err_ym, yerr=arm_eta_err_ym,
               fmt='o', c=c_palette[5])
ax['a'].set_xlim(-1200, 1200)
# Viscoelasticity
ax['b'].set_xlabel('Distance from CCA (nm)')
ax['b'].set_ylabel('Viscoelasticity index')
ax['b'].set_ylim(0, 1)
ax['b'].plot(tel_dist_bin, tel_eta_bin, color=c_palette[1], linestyle='solid', marker='o', label='Peripheral regions')
ax['b'].errorbar(tel_dist_bin, tel_eta_bin, xerr=tel_dist_err, yerr=tel_eta_err,
               fmt='o', c=c_palette[1])
ax['b'].plot(cent_dist_bin, cent_eta_bin, color=c_palette[3], linestyle='solid', marker='o', label='Centromeric regions')
ax['b'].errorbar(cent_dist_bin, cent_eta_bin, xerr=cent_dist_err, yerr=cent_eta_err,
               fmt='o', c=c_palette[3])
ax['b'].plot(arm_dist_bin, arm_eta_bin, c=c_palette[5], linestyle='solid', marker='o', label='Mid-arm regions')
ax['b'].errorbar(arm_dist_bin, arm_eta_bin, xerr=arm_dist_err, yerr=arm_eta_err,
               fmt='o', c=c_palette[5])
ax['b'].legend(frameon=False)
ax['b'].set_xlim(-1200, 1200)
# storing all the data of the binned profiles in a dataframe in order to use them elsewhere
df_profiles = pd.DataFrame({
    "tel_dist_ym": tel_dist_bin_ym,
    "tel_ym": tel_eta_bin_ym,
    "tel_dist_ym_err": tel_dist_err_ym,
    "tel_ym_err": tel_eta_err_ym,
    "cent_dist_ym": cent_dist_bin_ym,
    "cent_ym": cent_eta_bin_ym,
    "cent_dist_ym_err": cent_dist_err_ym,
    "cent_ym_err": cent_eta_err_ym,
    "arm_dist_ym": arm_dist_bin_ym,
    "arm_ym": arm_eta_bin_ym,
    "arm_dist_ym_err": arm_dist_err_ym,
    "arm_ym_err": arm_eta_err_ym,
    "tel_dist_eta": tel_dist_bin,
    "tel_eta": tel_eta_bin,
    "tel_dist_eta_err": tel_dist_err,
    "tel_eta_err": tel_eta_err,
    "cent_dist_eta": cent_dist_bin,
    "cent_eta": cent_eta_bin,
    "cent_dist_eta_err": cent_dist_err,
    "cent_eta_err": cent_eta_err,
    "arm_dist_eta": arm_dist_bin,
    "arm_eta": arm_eta_bin,
    "arm_dist_eta_err": arm_dist_err,
    "arm_eta_err": arm_eta_err
})
# saving df_profiles to .csv
df_profiles.to_csv("Cent_Arm_Tel_BinnedProfiles_YM_eta.csv")
# saving the figure
plt.savefig("Figure_MechanicsAlongChromatids.png", dpi=500)

plt.tight_layout()
plt.show()

# ---- Supp Figure ----
# single column width = 8,9 cm, 2-column width = 18,2 cm
cm = 1/2.54  # centimeters in inches
fig = plt.figure(figsize=(18*cm, 12*cm))
ax = fig.subplot_mosaic([['a', 'b', 'c'],
                         ['d', 'e', 'f']])
# Annotations
for label, a in ax.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(1/72, 3/72, fig.dpi_scale_trans)
    a.text(0.0, 1.0, label, transform=a.transAxes + trans,
            fontsize='medium', weight='bold', va='bottom', fontfamily='sans-serif')
# Young Modulus
ax['a'].set_title('Peripheral regions')
ax['a'].set_xlabel('')
ax['a'].set_ylabel('YM (kPa)')
a = sns.scatterplot(x=df_telomere['dist'], y=df_telomere['YM'], alpha=0.15, ax=ax['a'], color=c_palette[1])
ax['a'].errorbar(tel_dist_bin_ym, tel_eta_bin_ym, xerr=tel_dist_err_ym, yerr=tel_eta_err_ym,
               fmt='o', c=c_palette[1])
ax['a'].plot(tel_dist_bin_ym, tel_eta_bin_ym, color=c_palette[1])
a.set(xlabel="")

ax['d'].set_ylim(0, 1)
ax['d'].set_xlabel('Distance from CCA (nm)')
ax['d'].set_ylabel('Viscoelasticity index')
d = sns.scatterplot(x=df_telomere_eta['dist'], y=df_telomere_eta['eta'], alpha=0.15, ax=ax['d'], color=c_palette[1])
ax['d'].errorbar(tel_dist_bin, tel_eta_bin, xerr=tel_dist_err, yerr=tel_eta_err,
                fmt='o', c=c_palette[1])
ax['d'].plot(tel_dist_bin, tel_eta_bin, c=c_palette[1])
# ---- CENTROMERES ----
ax['b'].set_title('Centromeric regions')
ax['b'].set_xlabel('')
ax['b'].set_ylabel('')
b = sns.scatterplot(x=df_centromere['dist'], y=df_centromere['YM'], alpha=0.15, ax=ax['b'], color=c_palette[3])
ax['b'].errorbar(cent_dist_bin_ym, cent_eta_bin_ym, xerr=cent_dist_err_ym, yerr=cent_eta_err_ym,
               fmt='o', c=c_palette[3])
ax['b'].plot(cent_dist_bin_ym, cent_eta_bin_ym, c=c_palette[3])
ax['b'].plot(cent_dist_bin_ym, cent_eta_bin_ym, c=c_palette[3])
b.set(xlabel="", ylabel="")

ax['e'].set_ylim(0, 1)
ax['e'].set_xlabel('Distance from CCA (nm)')
ax['e'].set_ylabel('')
e = sns.scatterplot(x=df_centromere_eta['dist'], y=df_centromere_eta['eta'], alpha=0.15, ax=ax['e'], color=c_palette[3])

ax['e'].errorbar(cent_dist_bin, cent_eta_bin, xerr=cent_dist_err, yerr=cent_eta_err,
               fmt='o', c=c_palette[3])
ax['e'].plot(cent_dist_bin, cent_eta_bin, c=c_palette[3])
e.set(ylabel="")
# ---- ARMS ----
ax['c'].set_title('Mid-arm regions')
ax['c'].set_xlabel('')
ax['c'].set_ylabel('')
c = sns.scatterplot(x=df_arms['dist'], y=df_arms['YM'], alpha=0.15, ax=ax['c'], color=c_palette[5])
ax['c'].errorbar(arm_dist_bin_ym, arm_eta_bin_ym, xerr=arm_dist_err_ym, yerr=arm_eta_err_ym,
               fmt='o', c=c_palette[5])
ax['c'].plot(arm_dist_bin_ym, arm_eta_bin_ym, c=c_palette[5])
ax['c'].plot(arm_dist_bin_ym, arm_eta_bin_ym, c=c_palette[5])
c.set(xlabel="", ylabel="")

ax['f'].set_ylim(0, 1)
ax['f'].set_xlabel('Distance from CCA (nm)')
ax['f'].set_ylabel('')
f = sns.scatterplot(x=df_arms_eta['dist'], y=df_arms_eta['eta'], alpha=0.15, ax=ax['f'], color=c_palette[5])
ax['f'].errorbar(arm_dist_bin, arm_eta_bin, xerr=arm_dist_err, yerr=arm_eta_err,
               fmt='o', c=c_palette[5])
ax['f'].plot(arm_dist_bin, arm_eta_bin, c=c_palette[5])
f.set(ylabel="")

# saving the figure
plt.savefig("SuppFig_MechanicsAlongChromatids_AllPoints.svg", dpi=800)

plt.tight_layout()
plt.show()
