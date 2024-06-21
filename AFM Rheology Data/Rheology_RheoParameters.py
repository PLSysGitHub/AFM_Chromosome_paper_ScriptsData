import pandas as pd
pd.options.display.max_columns = 999
import os

'''Script for extracting the results of the AFM-rheological analysis on the chromosomes. 
   The script reads the rheology analysis files generated by the JPK-Bruker software for 
   each chromosome and stack all them together into a single .txt file called 'Rheology_results_allChrms.txt'.
   
   Each rheological analysis file is a .tsv that contains the following information:
   ['Filename', 'Position Index', 'X Position', 'Y Position', 'Baseline Offset [N]', 
   'Baseline Slope', 'Contact Point Offset [m]', 'Oscillation Frequency [Hz]', 
   'Force Amplitude [N]', 'Force Phase [deg]', 'Indentation Amplitude [m]', 
   'Indentation Phase [deg]', 'Indentation Mean [m]', 'Phase Offset [deg]', 
   'Phase Delta [deg]', 'Youngs Modulus [Pa]', 'Storage Modulus E' [Pa]', 
   'Loss Modulus E'' [Pa]', 'Shear Storage G' [Pa]', 'Shear Loss G'' [Pa]', 
   'Loss Tan', 'Drag Correction [Pa]'] '''

# Grouping all the directories of the file containing the results of the analyses
root = "Probed_chromosomes"
all_files = []
for path, subdirs, files in os.walk(root):
    for name in files:
        file = os.path.join(path, name)
        if ".tsv" in file:
            all_files.append(file)
# Storing the useful parameters of all the probed chromosomes within a single dataframe
df = pd.DataFrame()
for path in all_files:
    df_raw = pd.read_csv(path, sep='\t')
    # Changing the name of the chromosome in the Filename col
    df_raw["Filename"] = [item.split('_')[0] for item in df_raw["Filename"]]
    # Keeping only the meaningful information
    df_curr = df_raw[['Filename', 'Contact Point Offset [m]',
                      'Oscillation Frequency [Hz]', 'Force Amplitude [N]',
                      'Phase Delta [deg]', 'Indentation Amplitude [m]',
                      'Force Phase [deg]',
                      'Youngs Modulus [Pa]', "Storage Modulus E' [Pa]",
                      "Loss Modulus E'' [Pa]", "Shear Storage G' [Pa]",
                      "Shear Loss G'' [Pa]", 'Loss Tan']].copy()
    df_curr["Contact Point Offset [m]"] = df_curr["Contact Point Offset [m]"] * 1e9
    df_curr.rename(columns={"Contact Point Offset [m]": "Contact Point [nm]",
                            'Youngs Modulus [Pa]': 'Young Modulus [Pa]'}, inplace=True)
    df = pd.concat([df, df_curr], ignore_index=True)
# print(df.columns)

# Discarding indentations on the glass surface based on their Force Amplitude value
df = df[df["Force Amplitude [N]"] < .4e-9]
# print(df["Filename"].unique())


# Saving the DataFrame as a .txt file (uncomment to save and create the file)
# df.to_csv("Rheology_results_allChrms.txt", sep='\t')

