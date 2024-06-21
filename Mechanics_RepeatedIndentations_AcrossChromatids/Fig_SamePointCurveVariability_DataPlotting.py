import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker
import pandas as pd
import numpy as np
import seaborn as sns



df1 = pd.read_csv("Point_variability_analyses/point3_analysis.txt", sep='\t')
df1['point_n'] = 'point 1'
df1['YM_Hertz'] = df1['YM_Hertz']*1000
df1['contact_pt_var'] = df1['contact_pt'][0] - df1['contact_pt']
df1['YM_Hertz_var'] = 1 - abs(df1['YM_Hertz'][0]*1000 / df1['YM_Hertz']*1000)
df1['viscoel_index_var'] = 1 - abs(df1['viscoel_index'][0] / df1['viscoel_index'])
df2 = pd.read_csv("Point_variability_analyses/point1_analysis.txt", sep='\t')
df2['point_n'] = 'point 2'
df2['YM_Hertz'] = df2['YM_Hertz']*1000
df2['contact_pt_var'] = df2['contact_pt'][0] - df2['contact_pt']
df2['YM_Hertz_var'] = 1 - abs(df2['YM_Hertz'][0]*1000 / df2['YM_Hertz']*1000)
df2['viscoel_index_var'] = 1 - abs(df2['viscoel_index'][0] / df2['viscoel_index'])
df3 = pd.read_csv("Point_variability_analyses/point2_analysis.txt", sep='\t')
df3['point_n'] = 'point 3'
df3['YM_Hertz'] = df3['YM_Hertz']*1000
df3['contact_pt_var'] = df3['contact_pt'][0] - df3['contact_pt']
df3['YM_Hertz_var'] = 1 - abs(df3['YM_Hertz'][0]*1000 / df3['YM_Hertz']*1000)
df3['viscoel_index_var'] = 1 - abs(df3['viscoel_index'][0] / df3['viscoel_index'])
df4 = pd.read_csv("Point_variability_analyses/point4_analysis.txt", sep='\t')
df4['point_n'] = 'point 4'
df4['YM_Hertz'] = df4['YM_Hertz']*1000
df4['contact_pt_var'] = df4['contact_pt'][0] - df4['contact_pt']
df4['YM_Hertz_var'] = 1 - abs(df4['YM_Hertz'][0]*1000 / df4['YM_Hertz']*1000)
df4['viscoel_index_var'] = 1 - abs(df4['viscoel_index'][0] / df4['viscoel_index'])
df = pd.concat([df4, df1, df2, df3])
df_all = [df1, df2, df3, df4]

# ---- Plotting ----
# single column widht = 8,9 cm, 2-column width = 18,2 cm
cm = 1/2.54  # centimeters in inches
fig = plt.figure(figsize=(18*cm, 15*cm))
ax = fig.subplot_mosaic([['a', 'b'],
                         ['a', 'c'],
                         ['a', 'd']], )
for label, a in ax.items():
    # label physical distance to the left and up:
    trans = mtransforms.ScaledTranslation(1/72, 3/72, fig.dpi_scale_trans)
    if label == 'E':
        continue
    a.text(0.0, 1.0, label, transform=a.transAxes + trans,
        fontsize='large', weight='bold', va='bottom', fontfamily='sans-serif')

# ax['E'].axis('off')

color = sns.color_palette("viridis")
df_curve = pd.read_csv('Representative_curve_Fit_area_ContacPt_data.txt', sep='\t')
print(df_curve)
# ax['a'].grid()
ax['a'].set_xlabel('Separation (nm)')
ax['a'].set_ylabel('Force (nN)')
ax['a'].set_xlim(0, 600)
sns.lineplot(x=df_curve['Rt_separation'], y=df_curve['Rt_force'] / 1000, linewidth=3,
             color=color[2], ax=ax['a'])
sns.lineplot(x=df_curve['Ex_separation'], y=df_curve['Ex_force'] / 1000, linewidth=3,
             color=color[4], ax=ax['a'])
ax['a'].fill_between(df_curve['Ex_separation'], df_curve['Ex_force'] / 1000, 0,
                      color='gray', alpha=0.3, hatch='\\')
ax['a'].fill_between(df_curve['Rt_separation'], df_curve['Rt_force'] / 1000, 0,
                      color='gray', alpha=0.3, hatch='x')
ax['a'].fill_between(df_curve['Rt_separation'], 0, df_curve['Rt_force'] / 1000,
                     where=(df_curve['Rt_force'] / 1000 <= 0), color='white', alpha=1)
sns.lineplot(x=df_curve['fit_x'].dropna() + df_curve['contact_pt(x,y)'][0], y=df_curve['fit_y'].dropna() / 1000,
             linewidth=2, linestyle='dashed', color='red', ax=ax['a'])
ax['a'].scatter(x=df_curve['contact_pt(x,y)'][0], y=df_curve['contact_pt(x,y)'][1] / 1000, color='red',
                 marker='X', s=45, zorder=2)
# Hide the right and top spines
ax['a'].spines[['right', 'top']].set_visible(False)

i = 0
color = sns.color_palette("RdYlBu_r")
for df in df_all:
    plot1 = sns.regplot(data=df, x='curve_n', y='contact_pt', ax=ax['b'],
                        color=color[i])
    a = (np.abs(df["contact_pt"] - df["contact_pt"][0])/ df["contact_pt"][0]) * 100
    print(np.mean(a))
    i += 1
# ax['b'].set_xlabel('Indentation n.')
ax['b'].set_ylabel('Contact point (nm)')
ax['b'].set_xlim(0, 30)
ax['b'].set_ylim(200, 600)
ax['b'].set_xlabel('')
ax['b'].set_xticks([])
ax['b'].spines[['right', 'top', 'bottom']].set_visible(False)
# ax['c'].set_title('Young Modulus vs curve number')
i = 0
for df in df_all:
    plot2 = sns.regplot(data=df, x='curve_n', y='YM_Hertz', ax=ax['c'],
                        color=color[i], label=None)
    i += 1
# ax['c'].set_xlabel('Indentation n.')

ax['c'].set_ylabel('Young Modulus (kPa)')
ax['c'].set_xlim(0, 30)
ax['c'].set_ylim(0)
ax['c'].set_xlabel('')
ax['c'].set_xticks([])
ax['c'].spines[['right', 'top', 'bottom']].set_visible(False)
# ax['d'].set_title('eta vs curve number')
i = 0
for df in df_all:
    curr_label = 'Chromosome ' + str(i+1)
    plot3 = sns.regplot(data=df, x='curve_n', y='viscoel_index', ax=ax['d'],
                        color=color[i], label=curr_label)
    i += 1
ax['d'].set_xlabel('Indentation n.')
ax['d'].set_ylabel('Viscoelasticity index')
ax['d'].set_xlim(0, 30)
ax['d'].set_ylim(0)
ax['d'].spines[['right', 'top']].set_visible(False)
# Define a custom formatting function
def format_y(value, pos):
    return f'{value:.1f}'  # Format value to two decimal places
# Apply the custom formatter to the y-axis
ax['d'].yaxis.set_major_formatter(ticker.FuncFormatter(format_y))
plt.legend(frameon=False, loc='lower right', ncol=2, fontsize=8)

plt.savefig("Fig_RepeatedIndentations.svg", dpi=800)

plt.tight_layout()
plt.show()

