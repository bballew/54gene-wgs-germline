import create_exclude_list as ex
import pandas as pd
import pytest

bcf_columns = [
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
]
verifybamid_cols = [
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
]
out_cols = ["sample", "exclude_reason"]


def test_exclude_high_het_hom_pass():
    # both samples pass
    test_df = pd.DataFrame(
        [
            (
                "PSC",
                0,
                "sample-0049",
                28290884,
                1575665,
                3330169,
                3248597,
                1657237,
                1207087,
                26.7,
                95340,
                0,
                0,
                804565,
            ),
            (
                "PSC",
                0,
                "sample-0012",
                26838824,
                1643045,
                3239289,
                4172939,
                2106139,
                1456745,
                27.9,
                122761,
                0,
                0,
                600877,
            ),
        ],
        columns=bcf_columns,
    )
    test_out = ex.exclude_high_het_hom(test_df, 2.5)
    expected_out = None
    assert test_out == expected_out


def test_exclude_high_het_hom_fail():
    # one sample fails
    test_df = pd.DataFrame(
        [
            (
                "PSC",
                0,
                "sample-0049",
                28290884,
                974913,
                5304165,
                3248597,
                1657237,
                1207087,
                26.7,
                95340,
                0,
                0,
                804565,
            ),
            (
                "PSC",
                0,
                "sample-0012",
                26838824,
                1575665,
                3330169,
                4172939,
                2106139,
                1456745,
                27.9,
                122761,
                0,
                0,
                600877,
            ),
        ],
        columns=bcf_columns,
    )
    test_out = ex.exclude_high_het_hom(test_df, 2.5)
    expected_out = pd.DataFrame([("sample-0049", "high_het_hom")], columns=out_cols)
    assert pd.testing.assert_frame_equal(test_out, expected_out) is None


def test_exclude_depth_pass():
    test_df = pd.DataFrame(
        [
            (
                "PSC",
                0,
                "sample-0049",
                28290884,
                1575665,
                3330169,
                3248597,
                1657237,
                1207087,
                26.7,
                95340,
                0,
                0,
                804565,
            ),
            (
                "PSC",
                0,
                "sample-0012",
                26838824,
                1643045,
                3239289,
                4172939,
                2106139,
                1456745,
                27.9,
                122761,
                0,
                0,
                600877,
            ),
        ],
        columns=bcf_columns,
    )
    test_out = ex.exclude_low_depth(test_df, 20)
    expected_out = None
    assert test_out == expected_out


def test_exclude_depth_fail():
    test_df = pd.DataFrame(
        [
            (
                "PSC",
                0,
                "sample-0049",
                28290884,
                1575665,
                3330169,
                3248597,
                1657237,
                1207087,
                16.7,
                95340,
                0,
                0,
                804565,
            ),
            (
                "PSC",
                0,
                "sample-0012",
                26838824,
                1643045,
                3239289,
                4172939,
                2106139,
                1456745,
                27.9,
                122761,
                0,
                0,
                600877,
            ),
        ],
        columns=bcf_columns,
    )
    test_out = ex.exclude_low_depth(test_df, 20)
    expected_out = pd.DataFrame([("sample-0049", "low_depth")], columns=out_cols)
    assert pd.testing.assert_frame_equal(test_out, expected_out) is None


def test_exclude_contam_pass():
    test_df = pd.DataFrame(
        [
            (
                "file1:sample-0295",
                "ALL",
                "sample-0295",
                1931410,
                36724056,
                19.01,
                0.00967,
                7749750.48,
                7767800.39,
                "NA",
                "NA",
                0.01478,
                7054259.63,
                7094649.09,
                "NA",
                "NA",
                19.301,
                1.0170,
                0.9929,
            ),
            (
                "file2:sample-0214",
                "ALL",
                "sample-0214",
                1931410,
                35895529,
                18.59,
                0.00915,
                7641448.19,
                7657447.68,
                "NA",
                "NA",
                0.01520,
                6931430.04,
                6971883.56,
                "NA",
                "NA",
                18.960,
                1.0210,
                0.9958,
            ),
        ],
        columns=verifybamid_cols,
    )
    test_out = ex.exclude_contam(test_df, 0.03)
    expected_out = None
    assert test_out == expected_out


def test_exclude_contam_fail():
    test_df = pd.DataFrame(
        [
            (
                "file1:sample-0295",
                "ALL",
                "sample-0295",
                1931410,
                36724056,
                19.01,
                0.10967,
                7749750.48,
                7767800.39,
                "NA",
                "NA",
                0.01478,
                7054259.63,
                7094649.09,
                "NA",
                "NA",
                19.301,
                1.0170,
                0.9929,
            ),
            (
                "file2:sample-0214",
                "ALL",
                "sample-0214",
                1931410,
                35895529,
                18.59,
                0.00915,
                7641448.19,
                7657447.68,
                "NA",
                "NA",
                0.01520,
                6931430.04,
                6971883.56,
                "NA",
                "NA",
                18.960,
                1.0210,
                0.9958,
            ),
        ],
        columns=verifybamid_cols,
    )
    test_out = ex.exclude_contam(test_df, 0.03)
    expected_out = pd.DataFrame([("sample-0295", "contamination")], columns=out_cols)
    assert pd.testing.assert_frame_equal(test_out, expected_out) is None
