import argparse
import collections
import numpy as np
import pandas as pd


import logging
from typing import List

logger = logging.getLogger(__name__)

# parse command line arguments
parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--input",
    nargs="+",
    help="List of SNP array files.",
    required=True,
)

parser.add_argument(
    "-o",
    "--output",
    help="Path to output file.",
    required=True,
)


def read_csv_files(file_paths: List[str]) -> list[pd.DataFrame]:
    """Reads multiple SNP array CSV files into a list of dataframes.

    Args:
        file_paths (List[str]): List of paths to SNP array CSV files.

    Returns:
        list[pd.DataFrame]: List of dataframes containing SNP array data.
    """

    dfs = []
    for file_path in file_paths:
        df = pd.read_csv(file_path, sep="\t")
        dfs.append(df)
    return dfs


def check_input_dfs(dfs: list[pd.DataFrame]) -> None:
    """Check that the "Probeset ID" column is identical for each dataframe provided.
    To avoid merging samples form different SNP arrays.
    """

    probeset_ids = []
    for df in dfs:
        probeset_ids.append(df["Probeset ID"].values)

    if all([np.array_equal(probeset_ids[0], x) for x in probeset_ids[1:]]) == False:
        raise ValueError(
            "Probeset IDs are not identical for all input files. Check you have not mixed SNP arrays data."
        )

    column_names = [col_name for df in dfs for col_name in df.columns.values]
    column_names = [
        col_name
        for col_name in column_names
        if col_name
        not in {
            "Probeset ID",
            "Chr",
            "Position",
        }  # Filter out the column names that are expected to be duplicated
    ]

    # Check that there are no other duplicated columns
    if len(column_names) != len(set(column_names)):
        raise ValueError(
            f"Duplicate columns found in input files, {[item for item, count in collections.Counter(column_names).items() if count > 1]}. Check you have not provided duplicate files."
        )


def merge_array_files(dfs_to_merge: list[pd.DataFrame]) -> pd.DataFrame:
    """Merges SNP array files into one file.  Due to limitations with the software used the
    SNP array CSV files produced can only contain a limited number of SNPs.  Large runs are
    therefore split over multiple files which need to be merged into one file before passing
    to BASHer.

    Args:
        dfs_to_merge (list[Dataframes]): List of dataframes to merge from imported from SNP array files.

    Returns:
        pd.DataFrame: Merged dataframe.
    """
    check_input_dfs(dfs_to_merge)
    result = dfs_to_merge[0]
    for df in dfs_to_merge[1:]:
        result = pd.merge(
            result,
            df.loc[:, ~df.columns.isin(["Chr", "Position"])],
            how="inner",
            on="Probeset ID",
            sort=False,
            copy=True,
            validate="one_to_one",
        )
    return result


def write_merged_file(merged_df: pd.DataFrame, merged_df_name: str) -> None:
    """
    Writes merged SNP array file to disk.
    Args:
        merged_df (pd.DataFrame): Merged dataframe.
        merged_df_name (str): Merged dataframe name.
    """

    merged_df.to_csv(
        merged_df_name,
        index=False,
        header=True,
        sep="\t",
    )


def main(input_files: List[str]) -> pd.DataFrame:
    input_dfs = read_csv_files(input_files)
    merged_df = merge_array_files(input_dfs)
    return merged_df


# run the script
if __name__ == "__main__":
    args = parser.parse_args()
    output_df = main(args.input)
    write_merged_file(output_df, args.output)
