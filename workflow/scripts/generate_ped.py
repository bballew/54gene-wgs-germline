#!/usr/bin/env python

import argparse as ap
import re
import shutil
import sys

import pandas as pd


def get_args() -> list:
    """Get command line arguments"""
    parser = ap.ArgumentParser(
        description="Creates a plink-formatted ped from a tsv with subject and sex"
    )
    parser.add_argument("infile", type=str, help="linker filename")
    parser.add_argument(
        "outprefix", type=str, help="output filename prefix (.ped extension will be added)"
    )
    results = parser.parse_args()
    return (results.infile, results.outprefix)


def check_if_ped(infile: str) -> bool:
    """Look for .ped file extension

    This script can either be provided a tsv, which it will
    convert to a plink-style ped, or it can be provided a ped
    which it will simply copy to the expected output name/
    location.
    """
    if re.search(".ped$", infile, re.IGNORECASE):
        return True
    else:
        return False


def read_in_linker(infile: str) -> pd.DataFrame:
    """Read tsv linker file into a dataframe"""
    df = pd.read_table(infile, sep="\t")
    return df


def check_column_number(df: pd.DataFrame):
    """Check for two and only two columns in input tsv"""
    if len(df.columns) != 2:
        raise


def check_column_headers(df: pd.DataFrame):
    """Check for required headers in input tsv"""
    if ["Sample", "Sex"] != list(df.columns)[0:2]:
        raise


def add_ped_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add dummy columns to conform to plink ped format

    Should end up with a dataframe like so:
    Sample Sample 0 0 [F|M] -9
    """
    df.insert(1, "a", df.loc[:, "Sample"])
    df.insert(2, "b", 0)
    df.insert(2, "c", 0)
    df["d"] = "-9"
    return df


def check_sex_values(df: pd.DataFrame) -> pd.DataFrame:
    """Check for allowed values for sex information

    Input tsv may only have (case-insensitive) f, female,
    m, or male to encode sex.  Missing data can be encoded
    by one of the standard NA values that pandas can deal
    with automatically.  The error message emitted when this
    error is raised lists out the allowed values.
    """
    df["Sex"] = df["Sex"].str.lower()
    if not all(df[~df["Sex"].isna()].isin(["f", "female", "m", "male"])["Sex"]):
        raise


def encode_sex(df: pd.DataFrame) -> pd.DataFrame:
    """Convert sex information to ped encoding

    Missing sex information is converted to 0; male and
    female are converted to 1 and 2 respectively.
    """
    df["Sex"] = df["Sex"].str.lower()
    df["Sex"] = df["Sex"].replace(["f", "female", "m", "male"], [2, 2, 1, 1])
    df = df.fillna(0)
    df["Sex"] = df["Sex"].astype(int)
    return df


if __name__ == "__main__":
    infile, outprefix = get_args()
    if check_if_ped(infile):
        shutil.copyfile(infile, outprefix + ".ped")
        sys.exit(0)
    df = read_in_linker(infile)
    try:
        check_column_number(df)
    except Exception:
        sys.exit("Error: {} is required to have exactly two columns.".format(infile))
    try:
        check_column_headers(df)
    except Exception:
        sys.exit(
            "Error: Correct headers not detected in {} (requires 'Sample' and 'Sex', in that order).".format(
                infile
            )
        )
    try:
        check_sex_values(df)
    except Exception:
        sys.exit(
            "Error: Sex reported in {} must be represented by f, female, m, or male (case insensitive).  Missing data can be represented by any of the following strings: '', '#N/A', '#N/A N/A', '#NA', '-1.#IND', '-1.#QNAN', '-NaN', '-nan', '1.#IND', '1.#QNAN', '<NA>', 'N/A', 'NA', 'NULL', 'NaN', 'n/a', 'nan', or 'null'.".format(
                infile
            )
        )

    df = add_ped_columns(df)
    df = encode_sex(df)
    df.to_csv(outprefix + ".ped", sep="\t", header=False, index=False)
