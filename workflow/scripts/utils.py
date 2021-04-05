"""Utilities and helper functions for the WGS pipeline"""

import glob
import os
import sys

sampleDict = {}


def read_in_manifest(s, full):
    """Parse whitespace delimited manifest file.

    Note the distinct requirements for the two different run types.

    We can handle an arbitrary number of r1/r2 fastq pairs per sample,
    as long as the readgroups are distinct but the sample names are the same.
    """
    global sampleDict
    with open(s) as f:
        for line in f:
            if full:
                (rg, sm, read1, read2) = line.split()
                sampleDict[rg] = sm, read1, read2
            else:
                (sm, gvcf) = line.split()
                sampleDict[sm] = gvcf
    return sampleDict


def create_samples_set(d):
    """Generate a list of unique sample IDs from the manifest."""
    sm_set = set()
    for (key, value) in d.items():
        sm_set.add(value[0])
    return list(sm_set)


def _get_dict(**kwargs):
    """Return test dict for testing or global dict for workflow

    Functions called as snakemake input can't pass more than the wildcards
    parameter.  To enable testing, use **kwargs to optionally pass in a dict.
    Otherwise, for use in snakemake rules, use the globally-defined sample
    dict.
    """
    if kwargs:
        for d in kwargs.values():
            myDict = d.copy()
    else:
        myDict = sampleDict.copy()
    return myDict


def get_read1_fastq(wildcards, **kwargs):
    """Return r1 fastq for a given readgroup."""
    myDict = _get_dict(**kwargs)
    (sm, read1, read2) = myDict[wildcards.rg]
    return read1


def get_read2_fastq(wildcards, **kwargs):
    """Return r2 fastq for a given readgroup."""
    myDict = _get_dict(**kwargs)
    (sm, read1, read2) = myDict[wildcards.rg]
    return read2


def get_gvcf(wildcards, **kwargs):
    """"""
    myDict = _get_dict(**kwargs)
    gvcf = myDict[wildcards.sample]
    return gvcf


def get_gvcf_index(wildcards, **kwargs):
    """"""
    myDict = _get_dict(**kwargs)
    gvcf = myDict[wildcards.sample]
    return gvcf + ".tbi"


def get_sm(wildcards, **kwargs):
    """Return sample name for a given readgroup."""
    myDict = _get_dict(**kwargs)
    (sm, read1, read2) = myDict[wildcards.rg]
    return sm


def get_inputs_with_matching_SM(wildcards, **kwargs):
    """
    Lists full path and file name for each readgroup bam that
    has the same sample name (SM).  Used to generate input list
    for rule mark_duplicates.
    """
    myDict = _get_dict(**kwargs)
    l1 = []
    for (keys, vals) in myDict.items():
        if wildcards.sample in vals:
            l1.append("results/mapped/" + keys + ".bam")
    return l1


def list_markdup_inputs(wildcards, **kwargs):
    """List path and filename for each readgroup bam with the same sample name.
    Used to create input string for Picard MarkDuplicates (INPUT=/path/to/sample.bam
    INPUT=/path/to/sample2.bam ...).  This allows merging of bams when there are
    multiple fastq pairs per sample, resulting in multiple pre-dedup bams per sample.
    """
    myDict = _get_dict(**kwargs)
    l1 = []
    for (keys, vals) in myDict.items():
        if wildcards.sample in vals:
            l1.append("results/mapped/" + keys + ".bam")
    s1 = " INPUT=".join(l1)
    return s1


def get_batch_limit_number(jobs, n):
    """Assign weight to limit concurrent running rules
    Required to avoid hitting throughput limits of fsx.  The limit for resources
    is set to the max number of jobs in the config and wrapper script, and if n=10,
    the rule weight is set to max jobs/10.

    Example:
        max jobs = 1000
        default batch resource = 1
        only want to run max 10 at a time

        1000 / 10 = 100 is assigned as the "batch" resource for alignment jobs,
        so only max 10 will fit under the limit at one time.

    """
    limit = round(int(jobs) / n)
    if limit == 0:
        limit = 1
    return limit


def get_chrom_list(bed):
    """Retrieve list of chromosomes from bed file in config

    The user provides a bed file (e.g. each chromosome's start and end for
    WGS, or a list of targeted regions for WES).  This function returns a
    list of unique chromosomes included in the bed for parallelization.
    """
    with open(bed) as file:
        chromList = list(set([line.split()[0] for line in file]))
    return sorted(chromList, key=_karyotypic_sort)


def _karyotypic_sort(c):
    """"""
    c = c.replace("chr", "")
    if c == "X":
        return 23
    elif c == "Y":
        return 24
    elif c in ("M", "MT"):
        return 25
    elif c.isdigit():
        return int(c)
    else:
        sys.exit("Non-canonical contigs detected in bed file.")


def get_DBImport_path1(wildcards):
    """Define input files for rule HC_genotype_gvcfs."""
    return glob.glob(
        "results/HaplotypeCaller/DBImport/"
        + wildcards.chrom
        + "/"
        + wildcards.chrom
        + "*/genomicsdb_meta_dir/genomicsdb_meta*.json"
    )


def get_DBImport_path2(wildcards):
    """Define input files for rule HC_genotype_gvcfs."""
    path = "".join(glob.glob("results/HaplotypeCaller/DBImport/" + wildcards.chrom + "/*/__*/"))
    myList = []
    if os.path.exists(path):
        myList = [
            "AD.tdb",
            "AD_var.tdb",
            "ALT.tdb",
            "ALT_var.tdb",
            "BaseQRankSum.tdb",
            "__book_keeping.tdb.gz",
            "__coords.tdb",
            "DP_FORMAT.tdb",
            "DP.tdb",
            "END.tdb",
            "ExcessHet.tdb",
            "FILTER.tdb",
            "FILTER_var.tdb",
            "GQ.tdb",
            "GT.tdb",
            "GT_var.tdb",
            "ID.tdb",
            "ID_var.tdb",
            "InbreedingCoeff.tdb",
            "MIN_DP.tdb",
            "MLEAC.tdb",
            "MLEAC_var.tdb",
            "MLEAF.tdb",
            "MLEAF_var.tdb",
            "MQRankSum.tdb",
            "PGT.tdb",
            "PGT_var.tdb",
            "PID.tdb",
            "PID_var.tdb",
            "PL.tdb",
            "PL_var.tdb",
            "QUAL.tdb",
            "RAW_MQandDP.tdb",
            "ReadPosRankSum.tdb",
            "REF.tdb",
            "REF_var.tdb",
            "SB.tdb",
            "__tiledb_fragment.tdb",
        ]
        myList = [path + file for file in myList]
    return myList


def allow_blanks(c):
    return c if c is not None else ""
