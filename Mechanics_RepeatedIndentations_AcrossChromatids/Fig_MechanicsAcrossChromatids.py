import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.transforms as mtransforms
from matplotlib.gridspec import GridSpec
import seaborn as sns
from scipy.stats import sem
from PIL import Image


# ---- Plotting ----
# single column width = 8,9 cm, 2-column width = 18,2 cm
cm = 1/2.54  # centimeters in inches
# fig = plt.figure(figsize=(18.2*cm, 15*cm))

fig = plt.figure(figsize=(18*cm, 20*cm))

gs = GridSpec(5, 6, figure=fig)
# gs = fig.add_gridspec(5, 6, hspace=0.4, wspace=0.4)
# # axes = gs.subplots()
ax = fig.subplot_mosaic([['a', 'b'],
                         ['c', 'd'],
                         ['e', 'f']])
for val in ax.values():
    val.axis("off")
ax['a'] = fig.add_subplot(gs[0:3, :4])
ax['b'] = fig.add_subplot(gs[0:3, 4:6])
ax['c'] = fig.add_subplot(gs[3:4, :3])
ax['d'] = fig.add_subplot(gs[3:4, 3:])
ax['e'] = fig.add_subplot(gs[4:, :3])
ax['f'] = fig.add_subplot(gs[4:, 3:])


for label, a in ax.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(1/72, 3/72, fig.dpi_scale_trans)
    a.text(0.0, 1.0, label, transform=a.transAxes + trans,
            fontsize='large', weight='bold', va='bottom', fontfamily='sans-serif')

    
img_path = "Chromosome_indentationGrid_CCA.tiff"
# creating image object
img1 = Image.open(img_path)
# using convert method for img1
# img1 = img.convert("L")
ax['b'].imshow(img1, aspect="equal")
ax['b'].xaxis.set_major_locator(ticker.NullLocator())
ax['b'].yaxis.set_major_locator(ticker.NullLocator())

# ---- Chromosome cross-section cartoon ----
img_path = "Chromatids_CrossSection.png"
# creating image object
img = Image.open(img_path)
# using convert method for img1
# img1 = img.convert("L")
width, height = img.size[0], img.size[1]
print(width, height)
# basewidth, hsize = int(width * 5), int(height * 5)
# img = img.resize((basewidth, hsize), Image.Resampling.LANCZOS)
print(img.size)
ax['a'].imshow(img, aspect="equal")
ax['a'].axis('off')

# ---- Custom function for binning all the data points ----
def lin_bin_YM(x, y, n_bins, err):
    boundaries = np.linspace(min(x), max(x), n_bins+1)
    x_m = []
    y_m = []
    x_er = []
    y_er = []
    ym_err = []
    for i in range(n_bins):
        sel = np.array((x > boundaries[i]) & (x <= boundaries[i+1]))
        ym_err = [err[i] for i in range(len(sel)) if sel[i]]
        x_m.append(np.average([x[i] for i in range(len(sel)) if sel[i]]))
        y_m.append(np.average([y[i] for i in range(len(sel)) if sel[i]]))
        x_er.append(sem([x[i] for i in range(len(sel)) if sel[i]]))
        y_er.append(sem([y[i] for i in range(len(sel)) if sel[i]]))
    return np.array(x_m), np.array(y_m), np.array(x_er), np.array(y_er)

''' Processing the results of the Young Modulus analysis'''

df10 = pd.read_csv("20240424_YM_AcrossChromatids_AllData10.txt", sep='\t')
df20 = pd.read_csv("20240424_YM_AcrossChromatids_AllData20.txt", sep='\t')
df = pd.concat([df10, df20], axis=0, ignore_index=True)

df = df[df['ALL_ym'] > 0]
df = df[df['ALL_err'] > 0]

print(df["ALL_ym"].mean())
print(df["ALL_err"].mean())
print(df.columns)

print("Average Young Modulus (kPa) for inner part: {} ± {}"
      .format(df.query('ALL_dist <= 700 | ALL_dist > -700')["ALL_ym"].mean().round(2),
              df.query('ALL_dist <= 700 | ALL_dist > -700')["ALL_ym"].sem().round(2)))

df_var = df[(df["ALL_dist"] <= 700) & (df["ALL_dist"] >= -700)].groupby("chrm", as_index=False)["ALL_ym"].std()
print(df_var)
print("Average Standard dev: " + str(df_var["ALL_ym"].mean()) + " ± " + str(df_var["ALL_ym"].sem()))


# ---- Young Modulus obtained via fit ----
n_bins = 15
binned_data = lin_bin_YM(np.array(df['ALL_dist']), np.array(df['ALL_ym']), n_bins, np.array(df['ALL_err']))
dist_bin, ym_bin, xerr_bin, yerr_bin = binned_data[0], binned_data[1], binned_data[2], binned_data[3]

''' Processing the results of the Viscoelasticity analysis'''

df10_eta = pd.read_csv("20240424_ViscIndex_AcrossChromatids_AllData10.txt", sep='\t')
df20_eta = pd.read_csv("20240424_ViscIndex_AcrossChromatids_AllData20.txt", sep='\t')
df_eta = pd.concat([df10_eta, df20_eta], axis=0, ignore_index=True)
df_eta = df_eta[df_eta['ALL_eta'] >= 0]
df_eta = df_eta[df_eta['ALL_eta'] <= 1]
# print(df_eta)

# ---- Custom function for binning all the data points ----
def lin_bin_eta(x, y, n_bins):
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

# ---- ViscIndex ----
n_bins = 15
binned_data = lin_bin_eta(np.array(df_eta['ALL_dist']), np.array(df_eta['ALL_eta']), n_bins)
dist_bin, ym_bin, xerr_bin, yerr_bin = binned_data[0], binned_data[1], binned_data[2], binned_data[3]

# ---- Plotting Young Modulus and Viscoelasticity index ----
c_palette = sns.color_palette('rocket')
palette = 'rocket'
# # # Binning and calculating the average of the binned data
binned_data = lin_bin_YM(np.array(df['ALL_dist']), np.array(df['ALL_ym'] / 1000), n_bins, np.array(df['ALL_err']))
dist_bin, ym_bin, xerr_bin, yerr_bin = binned_data[0], binned_data[1], binned_data[2], binned_data[3]
sns.scatterplot(x=df['ALL_dist'], y=df['ALL_ym'] / 1000, s=55, alpha=0.2, ax=ax['c'], c='darkkhaki')
ax['c'].errorbar(dist_bin, ym_bin, xerr=xerr_bin, yerr=yerr_bin, fmt='.-', c='green', markeredgewidth=3, linewidth=2)
sns.scatterplot(x=df[df["chrm"] == 'PS220410']['ALL_dist'], y=df[df["chrm"] == 'PS220410']['ALL_ym'] / 1000,
                s=55, ax=ax['c'], facecolors='none', edgecolor='lightcoral', alpha=1.0)
sns.scatterplot(x=df[df["chrm"] == 'PS091154']['ALL_dist'], y=df[df["chrm"] == 'PS091154']['ALL_ym'] / 1000,
                s=55, ax=ax['c'], facecolors='none', edgecolor='royalblue', alpha=1.0)
# ax['c'].set_xlabel("Distance from CCA (nm)")
ax['c'].set_xlabel('')
ax['c'].set_xticks([])
ax['c'].set_ylabel("Young Modulus (MPa)")
ymax = df['ALL_ym'].max()/1000
ax['c'].set_yscale('log')
# ax['c'].set_ylim(0, 3)
# ax['c'].set_xlim(-1500, 1500)
ax['c'].vlines(x=[-700, 700], ymin=0, ymax=3, colors='gray', ls='--')


# # # Binning and calculating the average of the binned data
binned_data = lin_bin_eta(np.array(df_eta['ALL_dist']), np.array(df_eta['ALL_eta']), n_bins)
dist_bin, ym_bin, xerr_bin, yerr_bin = binned_data[0], binned_data[1], binned_data[2], binned_data[3]
sns.scatterplot(x=df_eta['ALL_dist'], y=df_eta['ALL_eta'], s=55, alpha=0.2, ax=ax['d'], c='darkkhaki')
ax['d'].errorbar(dist_bin, ym_bin, xerr=xerr_bin, yerr=yerr_bin, fmt='.-', c='green', markeredgewidth=3, linewidth=2)
# ax['d'].set_xlabel("Distance from CCA (nm)")
ax['d'].set_xlabel('')
ax['d'].set_xticks([])
ax['d'].set_ylabel("Viscoelasticity index")
ax['d'].set_ylim(0, 1)
ax['d'].set_xlim(-1500, 1500)
ax['d'].vlines(x=[-700, 700], ymin=0, ymax=ymax, colors='gray', ls='--')

# adding the profiles for Peripheral, mid-arm and centromeric regions
# Palettes used in the figures
c_palette = sns.color_palette('rainbow')
palette = [c_palette[1], c_palette[3], c_palette[5]]
df_profiles = pd.read_csv("Cent_Arm_Tel_BinnedProfiles_YM_eta.csv", sep=",")
ax['e'].set_xlabel('Distance from CCA (nm)')
ax['e'].set_ylabel('Young Modulus (kPa)')
ax['e'].plot(df_profiles["tel_dist_ym"], df_profiles["tel_ym"],
             color=c_palette[1], linestyle='solid', marker='o', label='Peripheral regions')
ax['e'].errorbar(df_profiles["tel_dist_ym"], df_profiles["tel_ym"],
                 xerr=df_profiles["tel_dist_ym_err"], yerr=df_profiles["tel_ym_err"],
                 fmt='o', c=c_palette[1])
ax['e'].plot(df_profiles["cent_dist_ym"], df_profiles["cent_ym"],
             color=c_palette[3], linestyle='solid', marker='o', label='Centromeric regions')
ax['e'].errorbar(df_profiles["cent_dist_ym"], df_profiles["cent_ym"],
                 xerr=df_profiles["cent_dist_ym_err"], yerr=df_profiles["cent_ym_err"],
                 fmt='o', c=c_palette[3])
ax['e'].plot(df_profiles["arm_dist_ym"], df_profiles["arm_ym"],
             color=c_palette[5], linestyle='solid', marker='o', label='Mid-arm regions')
ax['e'].errorbar(df_profiles["arm_dist_ym"], df_profiles["arm_ym"],
                 xerr=df_profiles["arm_dist_ym_err"], yerr=df_profiles["arm_ym_err"],
                 fmt='o', c=c_palette[5])
ax['e'].set_xlim(-1500, 1500)
# Viscoelasticity
ax['f'].set_xlabel('Distance from CCA (nm)')
ax['f'].set_ylabel('Viscoelasticity index')
ax['f'].set_ylim(0, 1)
ax['f'].plot(df_profiles["tel_dist_eta"], df_profiles["tel_eta"],
             color=c_palette[1], linestyle='solid', marker='o', label='Peripheral regions')
ax['f'].errorbar(df_profiles["tel_dist_eta"], df_profiles["tel_eta"],
                 xerr=df_profiles["tel_dist_eta_err"], yerr=df_profiles["tel_eta_err"],
                 fmt='o', c=c_palette[1])
ax['f'].plot(df_profiles["cent_dist_eta"], df_profiles["cent_eta"],
             color=c_palette[3], linestyle='solid', marker='o', label='Centromeric regions')
ax['f'].errorbar(df_profiles["cent_dist_eta"], df_profiles["cent_eta"],
                 xerr=df_profiles["cent_dist_eta_err"], yerr=df_profiles["cent_eta_err"],
                 fmt='o', c=c_palette[3])
ax['f'].plot(df_profiles["arm_dist_eta"], df_profiles["arm_eta"],
             c=c_palette[5], linestyle='solid', marker='o', label='Mid-arm regions')
ax['f'].errorbar(df_profiles["arm_dist_eta"], df_profiles["arm_eta"],
                 xerr=df_profiles["arm_dist_eta_err"], yerr=df_profiles["arm_eta_err"],
                 fmt='o', c=c_palette[5])
ax['f'].legend(frameon=False)
ax['f'].set_xlim(-1500, 1500)
# saving the figure
plt.savefig("Fig_MechanicsAcrossChromatids.svg", dpi=800)

plt.tight_layout()
plt.show()