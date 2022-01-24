import pytest
import pandas as pd
from scripts import generate_ped

@pytest.fixture
def passing_linker():
    sample_data = {
        "Sample": ["Sample1", "Sample2"],
        "Sex": ["M", "F"]
    }
    return pd.DataFrame(sample_data)

@pytest.fixture
def dropped_column(passing_linker):
    passing_linker.drop('Sex', axis=1, inplace = True)
    return passing_linker

@pytest.fixture
def incorrect_headers(passing_linker):
    passing_linker.rename({"Sample" : "Sample_name"}, axis=1, inplace = True)
    return passing_linker

def ped_format():
    return pd.DataFrame({
        "Sample": ["Sample1", "Sample2"],
        "a": ["Sample1", "Sample2"],
        "c": [0, 0],
        "b": [0, 0],
        "Sex": ["M", "F"],
        "d": "-9"
    })

@pytest.fixture
def exp_ped_format():
    return pd.DataFrame({
        "Sample": ["Sample1", "Sample2"],
        "a": ["Sample1", "Sample2"],
        "c": [0, 0],
        "b": [0, 0],
        "Sex": ["M", "F"],
        "d": "-9"
    })

def encoded_sex(exp_ped_format):
    exp_ped_format.replace({"M": 1, "F": 2}, inplace=True)
    return exp_ped_format 

# testing functions
# check_if_ped
@pytest.mark.parametrize(
    "test_in, exp_out", [("testfile.ped", True), ("testfile.tsv", False), ("testfile.PED", True)],
)
def test_check_if_ped(test_in, exp_out):
    test_out = generate_ped.check_if_ped(test_in)
    assert exp_out == test_out

# read_in_linker

# check_column_number
def test_check_column_number(passing_linker):
    assert generate_ped.check_column_number(passing_linker) is None

def test_check_column_number_fail(dropped_column):
    with pytest.raises(Exception):
        generate_ped.check_column_number(dropped_column)

# check_column_headers
def test_check_column_headers(passing_linker):
    assert generate_ped.check_column_headers(passing_linker) is None

def test_check_column_headers_fail(incorrect_headers):
    with pytest.raises(Exception):
        generate_ped.check_column_headers(incorrect_headers)

# # add_ped_columns
def test_add_ped_columns(passing_linker):
     test_out = generate_ped.add_ped_columns(passing_linker)
     pd.testing.assert_frame_equal(test_out, ped_format())

# encode_sex
def test_encode_sex(exp_ped_format):
     test_out = generate_ped.encode_sex(exp_ped_format)
     pd.testing.assert_frame_equal(test_out, encoded_sex(exp_ped_format))

