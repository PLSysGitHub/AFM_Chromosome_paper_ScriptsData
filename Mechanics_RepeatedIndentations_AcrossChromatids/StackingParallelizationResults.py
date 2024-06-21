import pandas as pd
import glob  # this is for listing all the files in a given directory
import os

''' This script stacks all the results obtained from running the YM and/or
    ViscIndex analyses per each chromosome separately (using the 
    parallelization) into one single dataframe that can be then used for 
    further analyses and plotting'''


# Data are stored in "parallelization_results" folder
# path = "parallelization_results_20230628"
path = "parallelization_results_20240424"
all_results = [f for f in glob.glob(path + '**/*', recursive=True)]

YM_df_list_10 = [file for file in all_results if "_YM_AcrossChromatids_AllData10.txt" in file]
YM_df_list_20 = [file for file in all_results if "_YM_AcrossChromatids_AllData20.txt" in file]
print(YM_df_list_10, len(YM_df_list_10))
print(YM_df_list_20, len(YM_df_list_20))
#
df_YM_10 = pd.DataFrame()
for item in YM_df_list_10:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_YM_10 = pd.concat([df_YM_10, curr_df])
    df_YM_10.reset_index(drop=True, inplace=True)
print(df_YM_10, df_YM_10.shape)
df_YM_10.to_csv("20240424_YM_AcrossChromatids_AllData10.txt", sep='\t')


df_YM_20 = pd.DataFrame()
for item in YM_df_list_20:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_YM_20 = pd.concat([df_YM_20, curr_df])
    df_YM_20.reset_index(drop=True, inplace=True)
print(df_YM_20, df_YM_20.shape)
df_YM_20.to_csv("20240424_YM_AcrossChromatids_AllData20.txt", sep='\t')


ViscIndex_df_list_10 = [file for file in all_results if "_ViscIndex_AcrossChromatids_AllData10.txt" in file]
ViscIndex_df_list_20 = [file for file in all_results if "_ViscIndex_AcrossChromatids_AllData20.txt" in file]

df_ViscIndex_10 = pd.DataFrame()
for item in ViscIndex_df_list_10:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_ViscIndex_10 = pd.concat([df_ViscIndex_10, curr_df])
    df_ViscIndex_10.reset_index(drop=True, inplace=True)
print(df_ViscIndex_10, df_ViscIndex_10.shape)
df_ViscIndex_10.to_csv("20240424_ViscIndex_AcrossChromatids_AllData10.txt", sep='\t')

df_ViscIndex_20 = pd.DataFrame()
for item in ViscIndex_df_list_20:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_ViscIndex_20 = pd.concat([df_ViscIndex_20, curr_df])
    df_ViscIndex_20.reset_index(drop=True, inplace=True)
print(df_ViscIndex_20, df_ViscIndex_20.shape)
df_ViscIndex_20.to_csv("20240424_ViscIndex_AcrossChromatids_AllData20.txt", sep='\t')
