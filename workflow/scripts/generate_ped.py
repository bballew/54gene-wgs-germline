#!/usr/bin/env python

import argparse as ap
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
    """a"""
    if infile.endswith(".ped"):
        return True
    else:
        return False


def read_in_linker(infile: str) -> pd.DataFrame:
    """a"""
    df = pd.read_table(infile, sep="\t")
    return df


def check_column_number(df: pd.DataFrame):
    if len(df.columns) != 2:
        raise


def check_column_headers(df: pd.DataFrame):
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


def encode_sex(df: pd.DataFrame) -> pd.DataFrame:
    """a"""
    df["Sex"] = df["Sex"].str.lower()
    df["Sex"] = df["Sex"].replace(["f", "female", "m", "male"], [2, 2, 1, 1])
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

    df = add_ped_columns(df)
    df = encode_sex(df)
    df.to_csv(outprefix + ".ped", sep="\t", header=False, index=False)
