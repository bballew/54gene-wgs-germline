from unittest import mock

import pytest
from scripts import utils

sampleDict_full = {
    "sampleA_rg1": ["sampleA", "A1_r1.fq", "A1_r2.fq"],
    "sampleA_rg2": ["sampleA", "A2_r1.fq", "A2_r2.fq"],
    "sampleB_rg1": ["sampleB", "B_r1.fq", "B_r2.fq"],
}


sampleDict_jointgeno = {
    "sampleA": "sampleA.g.vcf.gz",
    "sampleB": "sampleB.g.vcf.gz",
}


class wildcards:
    rg = ""
    sample = ""


def test_read_in_manifest_full():
    m = mock.mock_open(
        read_data="rg1\tsample1\tsample1_r1.fq\tsample1_r2.fq\nrg2 sample2 sample2_r1.fq sample2_r2.fq"
    )
    with mock.patch("scripts.utils.open", m):
        test_out = utils.read_in_manifest("file", "full", "fastq_qc_only")
    exp_dict = {
        "rg1": ("sample1", "sample1_r1.fq", "sample1_r2.fq"),
        "rg2": ("sample2", "sample2_r1.fq", "sample2_r2.fq"),
    }
    assert test_out == exp_dict


# test with non-unique readgroups?

# def test_read_in_manifest_jointgeno():
#     m = mock.mock_open(read_data="sample1\tsample1.g.vcf.gz\nsample2 sample2.g.vcf.gz")
#     with mock.patch("utils.open", m):
#         test_out = utils.read_in_manifest("file", False)
#     exp_dict = {"sample1": "sample1.g.vcf.gz", "sample2": "sample2.g.vcf.gz"}
#     assert test_out == exp_dict


def test_create_samples_set():
    test_in = {
        "rg1": ("sample1", "x", "y"),
        "rg2": ("sample1", "x", "y"),
        "rg3": ("sample2", "x", "y"),
    }
    test_out = utils.create_samples_set(test_in)
    exp_out = ["sample1", "sample2"]
    assert set(test_out) == set(exp_out)


@pytest.mark.parametrize(
    "test_in, exp_out",
    [("sampleA_rg1", "A1_r1.fq"), ("sampleA_rg2", "A2_r1.fq"), ("sampleB_rg1", "B_r1.fq")],
)
def test_get_read1_fastq(test_in, exp_out):
    wildcards.rg = test_in
    assert utils.get_read1_fastq(wildcards, d=sampleDict_full) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out",
    [("sampleA_rg1", "A1_r2.fq"), ("sampleA_rg2", "A2_r2.fq"), ("sampleB_rg1", "B_r2.fq")],
)
def test_get_read2_fastq(test_in, exp_out):
    wildcards.rg = test_in
    assert utils.get_read2_fastq(wildcards, d=sampleDict_full) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out",
    [("sampleA_rg1", "sampleA"), ("sampleA_rg2", "sampleA"), ("sampleB_rg1", "sampleB")],
)
def test_get_sm(test_in, exp_out):
    wildcards.rg = test_in
    assert utils.get_sm(wildcards, d=sampleDict_full) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out", [("sampleA", "sampleA.g.vcf.gz"), ("sampleB", "sampleB.g.vcf.gz")],
)
def test_get_gvcf(test_in, exp_out):
    wildcards.sample = test_in
    assert utils.get_gvcf(wildcards, d=sampleDict_jointgeno) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out", [("sampleA", "sampleA.g.vcf.gz.tbi"), ("sampleB", "sampleB.g.vcf.gz.tbi")],
)
def test_get_gvcf_index(test_in, exp_out):
    wildcards.sample = test_in
    assert utils.get_gvcf_index(wildcards, d=sampleDict_jointgeno) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out",
    [
        ("sampleA", ["results/mapped/sampleA_rg1.bam", "results/mapped/sampleA_rg2.bam"]),
        ("sampleB", ["results/mapped/sampleB_rg1.bam"]),
    ],
)
def test_get_inputs_with_matching_SM(test_in, exp_out):
    wildcards.sample = test_in
    assert utils.get_inputs_with_matching_SM(wildcards, d=sampleDict_full) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out",
    [
        ("sampleA", "results/mapped/sampleA_rg1.bam INPUT=results/mapped/sampleA_rg2.bam"),
        ("sampleB", "results/mapped/sampleB_rg1.bam"),
    ],
)
def test_list_markdup_inputs(test_in, exp_out):
    wildcards.sample = test_in
    assert utils.list_markdup_inputs(wildcards, d=sampleDict_full) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out", [([100, 10], 10), ([1, 20], 1), ([10, 6], 2)],
)
def test_get_batch_limit_number(test_in, exp_out):
    assert utils.get_batch_limit_number(test_in[0], test_in[1]) == exp_out


@pytest.mark.parametrize(
    "test_in, exp_out",
    [("chr1", 1), ("chr10", 10), ("20", 20), ("chrM", 25), ("chrY", 24), ("X", 23)],
)
def test_karyotypic_sort(test_in, exp_out):
    assert utils._karyotypic_sort(test_in) == exp_out


def test_karyotypic_sort_exit(capsys):
    with pytest.raises(SystemExit):
        utils._karyotypic_sort("chrUn")
        test_out, err = capsys.readouterr()
        assert test_out == "Non-canonical contigs detected in bed file."


# def test_get_DBImport_path1(wildcards):
# def test_get_DBImport_path2(wildcards):


@pytest.mark.parametrize(
    "test_in, exp_out", [("-opt 1 -and another", "-opt 1 -and another"), (None, "")],
)
def test_allow_blanks(test_in, exp_out):
    assert utils.allow_blanks(test_in) == exp_out
