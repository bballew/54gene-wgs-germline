#!/usr/bin/env python3

import argparse
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd


def get_args():
    """Handle command line arguments"""
    parser = argparse.ArgumentParser(
        description="Generates a list of samples to exclude based on previously-run QC metrics."
    )
    parser.add_argument("bstats", type=str, help="bcftools stats file name")
    parser.add_argument("outfile", type=str, help="Output file prefix")
    parser.add_argument(
        "-v",
        "--verify",
        type=str,
        help="concatenated verifyBamID *.selfSM output files",
        default=False,
    )
    parser.add_argument(
        "-r", "--ratio", type=float, help="maximum allowed het/hom_alt ratio", default=2.5
    )
    parser.add_argument(
        "-d", "--depth", type=float, help="minimum allowed average depth", default=20.0
    )
    parser.add_argument(
        "-c", "--contam", type=float, help="maximum allowed contamination", default=0.03
    )
    results = parser.parse_args()
    return (
        results.bstats,
        results.outfile,
        results.verify,
        results.ratio,
        results.depth,
        results.contam,
    )


def read_in_bcftools_PSC(bstats: str) -> pd.DataFrame:
    """Read into a dataframe just the PSC portion of bcftools stats output"""
    skip: list[int] = [i for i, line in enumerate(open(bstats)) if not line.startswith("PSC")]
    return pd.read_csv(
        bstats,
        sep="\t",
        skiprows=skip,
        header=None,
        names=[
            "PSC",
            "id",
            "sample",
            "nRefHom",
            "nNonRefHom",
            "nHets",
            "nTransitions",
            "nTransversions",
            "nIndels",
            "average depth",
            "nSingletons",
            "nHapRef",
            "nHapAlt",
            "nMissing",
        ],
    )


def exclude_high_het_hom(df: pd.DataFrame, r: float) -> pd.DataFrame:
    """Exclude any samples with a het/hom ratio over 2.5"""
    df.loc[:, "het_hom_ratio"] = df["nHets"] / df["nNonRefHom"]
    df_het = df.loc[(df["het_hom_ratio"] > r)].copy()
    if df_het.empty:
        return pd.DataFrame({"sample": [], "exclude_reason": []})
    else:
        df_het.loc[:, "exclude_reason"] = "high_het_hom"
        return df_het[["sample", "exclude_reason"]]


def exclude_low_depth(df: pd.DataFrame, d: float) -> pd.DataFrame:
    """Setting a 20x avg depth cutoff

    This is reported by the lab as well, but this should catch samples
    that are borderline, that were accidentally included despite low
    coverage, or that are from outside collaborators."""
    df_depth = df.loc[df["average depth"] < d].copy()
    if df_depth.empty:
        return pd.DataFrame({"sample": [], "exclude_reason": []})
    else:
        df_depth.loc[:, "exclude_reason"] = "low_depth"
        return df_depth[["sample", "exclude_reason"]]


def read_in_verifybamid(verify: str) -> pd.DataFrame:
    """Requires single *.selfSM files to be concatenated without headers"""
    return pd.read_csv(
        verify,
        sep="\t",
        header=None,
        names=[
            "#SEQ_ID",
            "RG",
            "CHIP_ID",
            "#SNPS",
            "#READS",
            "AVG_DP",
            "FREEMIX",
            "FREELK1",
            "FREELK0",
            "FREE_RH",
            "FREE_RA",
            "CHIPMIX",
            "CHIPLK1",
            "CHIPLK0",
            "CHIP_RH",
            "CHIP_RA",
            "DPREF",
            "RDPHET",
            "RDPALT",
        ],
    )


def exclude_contam(df_v: pd.DataFrame, c: float) -> pd.DataFrame:
    """Setting a 3% contamination threshold"""
    df_contam = df_v.loc[df_v["FREEMIX"] > c].copy()
    if df_contam.empty:
        return pd.DataFrame({"sample": [], "exclude_reason": []})
    else:
        df_contam.loc[:, "exclude_reason"] = "contamination"
        t = df_contam["#SEQ_ID"].str.split(":", expand=True)
        df_contam["sample"] = t[1]
        return df_contam[["sample", "exclude_reason"]]


def combine_exclusions(df_list: list) -> pd.DataFrame:
    """"""
    all_df = pd.concat(df_list).groupby("sample")["exclude_reason"].apply(",".join).reset_index()
    return all_df


if __name__ == "__main__":
    bstats, outfile, verify, r, d, c = get_args()
    df = read_in_bcftools_PSC(bstats)
    exclude1 = exclude_high_het_hom(df, r)
    exclude2 = exclude_low_depth(df, d)
    if verify:
        df_v = read_in_verifybamid(verify)
        exclude3 = exclude_contam(df_v, c)
        exclude_df = combine_exclusions([exclude1, exclude2, exclude3])
    else:
        exclude_df = combine_exclusions([exclude1, exclude2])

    exclude_df.to_csv(outfile + "_with_annotation.tsv", sep="\t", header=None, index=None)
    exclude_df[["sample"]].to_csv(outfile + ".tsv", sep="\t", header=None, index=None)
