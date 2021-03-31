from unittest import mock

import pytest
import scripts.utils as utils


def test_read_in_manifest_full():
    m = mock.mock_open(
        read_data="rg1\tsample1\tsample1_r1.fq\tsample1_r2.fq\nrg2 sample2 sample2_r1.fq sample2_r2.fq"
    )
    with mock.patch("scripts.utils.open", m):
        test_out = utils.read_in_manifest("file", True)
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


# def test_get_read1_fastq(wildcards):
# def test_get_read2_fastq(wildcards):
# def test_get_gvcf(wildcards):
# def test_get_gvcf_index(wildcards):
# def test_get_sm(wildcards):
# def test_get_inputs_with_matching_SM(wildcards):
# def test_list_markdup_inputs(wildcards):


@pytest.mark.parametrize(
    "test_in, exp_out", [([100, 10], 10), ([1, 20], 1), ([10, 6], 2)],
)
def test_get_batch_limit_number(test_in, exp_out):
    assert utils.get_batch_limit_number(test_in[0], test_in[1]) == exp_out


def test_get_chrom_list():
    m = mock.mock_open(read_data="chr1 123 124\nchr20 456 457\nchr3 123 124\nchrY 123 154")
    with mock.patch("scripts.utils.open", m):
        test_out = utils.get_chrom_list("bed")
    assert test_out == ["chr1", "chr3", "chr20", "chrY"]


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
