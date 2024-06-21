import pandas as pd
import glob  # this is for listing all the files in a given directory
import os

''' This script stacks all the results obtained from running the YM and/or
    ViscIndex analyses per each chromosome separately (using the 
    parallelization) into one single dataframe that can be then used for 
    further analyses and plotting'''


# Data are stored in "parallelization_results" folder
path = "20240424"
# path = "parallelization_results_20230628"
all_results = [f for f in glob.glob(path + '**/*', recursive=True)]

YM_df_list = [file for file in all_results if "YM_Indentation_parts_values_Contact150nm.txt" in file]
df_YM = pd.DataFrame()
for item in YM_df_list:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_YM = pd.concat([df_YM, curr_df])
    df_YM.reset_index(drop=True, inplace=True)
print(df_YM, df_YM.shape)
df_YM.to_csv("20240424_YM_Indentation_parts_values_Contact150nm.txt", sep='\t')

ViscIndex_df_list = [file for file in all_results if "ViscIndex_Indentation_parts_values_Contact150nm.txt" in file]
df_ViscIndex = pd.DataFrame()
for item in ViscIndex_df_list:
    curr_df = pd.read_csv(item, sep='\t', index_col=0)
    df_ViscIndex = pd.concat([df_ViscIndex, curr_df])
    df_ViscIndex.reset_index(drop=True, inplace=True)
print(df_ViscIndex, df_ViscIndex.shape)
df_ViscIndex.to_csv("20240424_ViscIndex_Indentation_parts_values_Contact150nm.txt", sep='\t')
