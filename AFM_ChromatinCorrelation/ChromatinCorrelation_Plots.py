import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import seaborn as sns
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score


def ChromatinCorr(cuts_data, img_path, profile_type):
    """Reads the height profile .csv file and computes a rolling std for increasing windows sizes,
    returns the averaged std for the different windows sizes and applies an exponential fit to it,
    finally it prints the value of the scaling exponent."""
    # Reading and tidying up the .csv file containing the profiles
    df = pd.read_csv(cuts_data, sep=";", skiprows=3, header=None)
    df = df.iloc[:, :-1]
    header = []
    for i in range(1, len(df.columns)//2 + 1, 1):
        x_name = str(i) + "_x_(nm)"
        y_name = str(i) + "_y_(nm)"
        header.append(x_name)
        header.append(y_name)
    df.columns = header
    df = df*1e9

    vals = pd.DataFrame()
    # plt.ion()
    # Calculating the size of the starting window (3 points) and the
    # minimum window size increase "delta" (distance between two consecutive points)
    start_win_size = df.iloc[:, 0][2] - df.iloc[:, 0][0]
    delta_win_size = df.iloc[:, 0][1] - df.iloc[:, 0][0]
    print("Minimum window size increase (nm): ", delta_win_size)
    # Rolling sigma (std) on each column
    max_winsize = 60 # arbitrary, based on the length of the profiles
    min_size = 3
    step = 1
    for size in range(min_size, max_winsize, step):
        # Rolling std for all the columns in the df
        sigma_values = df.rolling(window=size, axis=0).std()
        # Calculating teh average for all the profiles
        mean_sigma = sigma_values.mean(skipna=True).values
        curr_df = pd.DataFrame({
            "sigma": mean_sigma[1::2],
            "profile": range(1, len(mean_sigma)//2 + 1, 1),
            "window": size,
            "window_size_nm": start_win_size
        })
        vals = pd.concat([vals, curr_df])
        start_win_size = start_win_size + step*delta_win_size
    # Fitting until window size < lin_lim, then no linear trend anymore
    lin_lim = 250
    fit_df = vals[vals["window_size_nm"] <= lin_lim]
    x = fit_df["window_size_nm"].unique()
    y = fit_df.groupby("window")["sigma"].mean()
    y_std = fit_df.groupby("window")["sigma"].std()

    def log_fit(x, theta, B):
        y = B * x ** (theta)
        return y

    # def lin_fit(x, theta, B):
    #     y = B + x * (theta)
    #     return y

    p, cov = curve_fit(log_fit, x, y, sigma=y_std, absolute_sigma=True)
    fit = log_fit(x, *p)
    err = np.sqrt(np.diag(cov))
    # p, cov = curve_fit(lin_fit, np.log(x), np.log(y))
    # fit = lin_fit(np.log(x), *p)

    print(f"Exponent: {p[0]}")
    exponent = p[0]
    print(f"error: {err[0]}")

    # r2 = r2_score(np.log(y), fit)
    # print(f"R-squared: {r2}")
    # img_path = "Chrms_PABuffer_OlympusTip.007_LongProfiles.tiff"


    # # Plotting each profile separately
    # fig, ax = plt.subplots()
    # for p in vals["profile"].unique():
    #     d = vals[vals["profile"] == p]
    #     ax.scatter(d["window_size_nm"], d["sigma"], label=p)
    # ax.set_xscale("log")
    # ax.set_yscale("log")
    # ax.set_xlabel("Window size (nm)")
    # ax.set_ylabel("$\sigma$")
    # plt.title(cuts_data)
    # plt.legend()
    # plt.show()

    # Plotting the average of all the profiles
    fig = plt.figure()
    ax = fig.subplots()
    # ax.errorbar(x=vals["window_size_nm"].unique(), y=vals.groupby("window")["sigma"].mean(),
    #                yerr=vals.groupby("window")["sigma"].sem(), marker='o', markersize=5,
    #                linestyle='', linewidth=2, label="Experimental $\\sigma$", c='dodgerblue',
    #                zorder=1, alpha=.85)
    # sns.lineplot(data=vals, x="window_size_nm", y="sigma", errorbar='se',
    #              err_style="bars", marker='o', markersize=8, ax=ax,
    #              c="dimgrey", label="Experimental $\\sigma$")
    # vals["window_size_nm"] = np.log(vals["window_size_nm"])
    # vals["sigma"] = np.log(vals["sigma"])
    sns.lineplot(data=vals, x="window_size_nm", y="sigma", errorbar='se',
                 err_style="bars", marker='o', markersize=8, ax=ax,
                 c="dimgrey", label="Experimental $\\sigma$")
    ax.plot(x, fit, c="firebrick", label="$\\sigma\\sim L^\\alpha$", linestyle='--',
            linewidth=2, zorder=2)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$L$ (nm)")
    ax.set_ylabel("$\sigma$ (nm)")
    ax.legend()
    ax.set_title(str(profile_type))
    plt.tight_layout()
    fig_name = profile_type + "_AvgRoughness_Fit.png"
    plt.savefig(fig_name, dpi=500)
    plt.show()

    return [vals["window_size_nm"].unique(), vals.groupby("window")["sigma"].mean(),
            vals.groupby("window")["sigma"].sem(), x, fit, exponent, img_path]


chrm_lp = "Chromosome/Chrms_PABuffer_OlympusTip.007_LongProfiles.csv"
chrm_lp_img = "Chromosome/Chrms_PABuffer_OlympusTip.007_LongProfiles.tiff"
chrm_cp = "Chromosome/Chrms_PABuffer_OlympusTip.007_CrossProfiles.csv"
chrm_cp_img = "Chromosome/Chrms_PABuffer_OlympusTip.007_CrossProfiles.tiff"


lp_var_x, lp_var_y, lp_var_yerr, lp_fit_x, lp_fit_y, lp_exponent, lp_img = \
    ChromatinCorr(chrm_lp, chrm_lp_img, "LongitudinalProfiles")

cp_var_x, cp_var_y, cp_var_yerr, cp_fit_x, cp_fit_y, cp_exponent, cp_img = \
    ChromatinCorr(chrm_cp, chrm_cp_img, "PerpendicularProfiles")


# Plotting Longitudinal and Perpendicular profiles together
# fig = plt.figure()
# ax = fig.subplots()
#
# # Longitudinal profiles
# ax.errorbar(x=lp_var_x, y=lp_var_y, yerr=lp_var_yerr, marker='o', markersize=3,
#             linestyle='', linewidth=1, capsize=2, capthick=1, label="rolling_std_lp",
#             c='navy', zorder=1)
# ax.plot(lp_fit_x, lp_fit_y, c="firebrick", label="fit_lp", linestyle='--', zorder=2)
# # Cross profiles
# ax.errorbar(x=cp_var_x, y=cp_var_y, yerr=cp_var_yerr, marker='o', markersize=3,
#             linestyle='', linewidth=1, capsize=2, capthick=1, label="rolling_std_cp",
#             c='lightskyblue', zorder=1)
# ax.plot(cp_fit_x, cp_fit_y, c="firebrick", label="fit_cp", linestyle='--', zorder=2)
#
# ax.set_xscale("log")
# ax.set_yscale("log")
# ax.set_xlabel("Window size (nm)")
# ax.set_ylabel("$\sigma$")
# ax.legend()
# plt.tight_layout()
# plt.show()
