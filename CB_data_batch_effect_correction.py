#!/usr/bin/env python3

"""
Script to correct for batch effects in Cord Blood data based on FACS date.
We observed batch effects in the FACS data, likely due to antibody changes and instrument drift over time.
This script performs a simple median-based correction for each phenotype based on a smoothing window of 200 FACS dates.

@author: Antton Lamarka
"""

import numpy as np
import polars as pl
import matplotlib.pyplot as plt
import sys

########################################################################################################################
# Data loading
########################################################################################################################

path_to_metadata_file = "/media/antton/cbio3/data/BloodVariome/Cord Blood/CB sample info_antton20230227.xlsx"

# Read metadata file
metadata_df = pl.read_excel(path_to_metadata_file, sheet_name="CB sample info",
                            read_csv_options={"sep": ",",
                                              "null_values": ["-", "no count", "X", "skipped", "no label", "?",
                                                              "skipped"]})

metadata_df = metadata_df.drop(["Freeze date",
                                "Hospital",
                                "Birth date",
                                "Gender",
                                "Cell nr (10^6)",
                                "Blood (ml)",
                                "Tube scale",
                                "stock -150",
                                "Notes for sample collection",
                                "96V plate",
                                "Total events",
                                "Notes for analysis",
                                ""])


path_to_pheno_file = "/home/antton/Dropbox (BNilsson lab)/2021-CordBlood/CordBlood.final.2022-09-05.txt"

# Read phenotype file
pheno_df = pl.read_csv(path_to_pheno_file, sep='\t', skip_rows_after_header=3,
                       null_values=["NA", "N/A", "n/a", "na", "NaN", "nan", "NAN", "no count", "few", "-", "-"],
                       dtypes={})
# Rename column "Phenotype_Aligater_name_cleaned" to "ID"
pheno_df = pheno_df.rename({"Phenotype_AliGater_name_cleaned": "ID"})

# Add prefix "CD" to all elements "ID" in metadata_df
metadata_df = metadata_df.with_column(pl.col("ID").str.replace("", "CB").alias("ID"))

# Merge metadata and phenotype dataframes by ID
merged_df = metadata_df.join(pheno_df, on="ID")

merged_df = merged_df.unique(subset=["ID"])

########################################################################################################################
# Data processing
########################################################################################################################

# Remove rows with missing values in "ID" column
merged_df = merged_df.drop_nulls("ID")
print(merged_df.head())

# Convert "Birth date" from string to pl.Date, and from pl.Date to number. Add 365 to the number for each year after 2018
merged_df = merged_df.with_column(
    pl.col("FACS date ").cast(pl.Utf8()).str.slice(0, 6).str.replace("", "20")
    .str.strptime(pl.Date, fmt="%Y%m%d", strict=False)
    .dt.ordinal_day() +
    (pl.col("FACS date ").cast(pl.Utf8()).str.slice(0, 6).str.replace("", "20")
     .str.strptime(pl.Date, fmt="%Y%m%d", strict=False).
     dt.year() - 2018) * 365)
print(merged_df.select([pl.col("ID"), pl.col("FACS date ")]))

########################################################################################################################
# Plotting
########################################################################################################################

output_path = "/home/antton/Tiny_Projects/Investigating_missing_CB_data/CB_sanity_check/plots/unadjusted_plots"


def plot_all_phenotypes(df, pheno_list, output_folder):
    # Plot phenotype vs FACS date
    for column in pheno_list:  # Make
        # New figure
        plt.figure()
        print(column)
        plt.scatter(df["FACS date "], df[column], s=0.45)
        if "/" in column:
            column = column.replace("/", "_div_")
        if column == "unconventional_CD45RA":
            break
        plt.savefig(output_folder + column + ".png")


########################################################################################################################
# Correction
########################################################################################################################

def print_status(percent):
    """Prints a status bar to the console. Used to show progress of the script when we iterate through a loop."""
    sys.stdout.write("%3d%%\r" % percent)
    sys.stdout.flush()


print(pheno_df.columns[2:-1])
pandas_pheno_df = pheno_df[pheno_df.columns[2:-1]]
corrected_output_df = merged_df.to_pandas()
# Set ID as index
corrected_output_df.set_index("ID", inplace=True)
print(corrected_output_df)

for i, pheno in enumerate(pheno_df.columns[2:-1]):
    # print(pheno)
    full_pheno_pandas_df = pheno_df[pheno].to_pandas()
    FACS_date_list = merged_df["FACS date "].unique().to_list()[1:]
    for j, flow_session in enumerate(FACS_date_list):
        if j < 100:  # First 100 elements
            window_start = FACS_date_list[0]
        else:
            window_start = FACS_date_list[j - 100]
        if j > (len(FACS_date_list) - 101):  # Last 100 elements
            window_end = FACS_date_list[-1]
        else:
            window_end = FACS_date_list[j + 100]

        smoothing_window_pandas_df = merged_df.select(
            pl.col(pheno).where(pl.col("FACS date ").is_between(window_start, window_end))).to_pandas()
        smoothing_window_pandas_df.index = merged_df.select(
            pl.col("ID").where(pl.col("FACS date ").is_between(window_start, window_end))).to_pandas()

        day_pheno_pandas_df = merged_df.select(pl.col(pheno).where(pl.col("FACS date ") == flow_session)).to_pandas()
        day_pheno_pandas_df.index = merged_df.select(
            pl.col("ID").where(pl.col("FACS date ") == flow_session)).to_pandas()

        # print(smoothing_window_pandas_df)
        # print(day_pheno_pandas_df)

        if np.isnan(day_pheno_pandas_df[pheno]).all():
            continue

        median_diff = smoothing_window_pandas_df[pheno].median() - day_pheno_pandas_df[pheno].median()

        for sample in day_pheno_pandas_df.index.to_list():
            sample = sample[0]
            # print(sample)
            # Extract the phenotype value
            og_value = merged_df.select(pl.col(pheno).where((pl.col("ID") == sample) & (pl.col("FACS date ") == flow_session)))[0, 0]
            # print(og_value)
            # Compute the corrected value
            if og_value:
                corrected_value = og_value + median_diff
            else:
                corrected_value = np.nan  # None
            # print(corrected_value)
            # Store the corrected value
            corrected_output_df.loc[sample, pheno] = corrected_value

    # After correcting, prune all values that are outside 4 standard deviations from the mean
    #upper_bound = corrected_output_df.loc[:, pheno].mean() + 4 * corrected_output_df.loc[:, pheno].std()
    #lower_bound = corrected_output_df.loc[:, pheno].mean() - 4 * corrected_output_df.loc[:, pheno].std()

    # Prune all values that are outside 4 standard deviations from the mean in phase 2
    # Those values should be set to nan
    #corrected_output_df.loc[(corrected_output_df.loc[:, pheno] > upper_bound), pheno] = np.nan
    #corrected_output_df.loc[(corrected_output_df.loc[:, pheno] < lower_bound), pheno] = np.nan

    # Print status every 10 phenotypes
    percentage = 100 * (i / len(pheno_df.columns[2:-1]))
    print_status(percentage)

output_path = "/home/antton/Tiny_Projects/Investigating_missing_CB_data/CB_sanity_check/plots/corrected_plots_pruned/"

plot_all_phenotypes(corrected_output_df, pheno_df.columns[2:-1], output_path)
corrected_output_df = corrected_output_df.drop("FACS date ", axis=1)

########################################################################################################################
# Saving
########################################################################################################################
pheno_file_headers_df = pl.read_csv(path_to_pheno_file, sep='\t')
# Drop every row after the 3rd row
pheno_file_headers_df = pheno_file_headers_df[0:3]
pheno_file_headers_df = pheno_file_headers_df.rename({"Phenotype_AliGater_name_cleaned": "ID"})

corrected_output_df.to_csv(
    "/home/antton/Tiny_Projects/Investigating_missing_CB_data/CB_sanity_check/batch_corrected_step1.txt", sep='\t', na_rep="nan")

corrected_output_df = pl.from_pandas(corrected_output_df)
out_df = pl.concat([pheno_file_headers_df, corrected_output_df])

out_df.write_csv(
    "/home/antton/Tiny_Projects/Investigating_missing_CB_data/CB_sanity_check/batch_corrected_final_CordBlood_2023-03-08.txt", sep='\t')
