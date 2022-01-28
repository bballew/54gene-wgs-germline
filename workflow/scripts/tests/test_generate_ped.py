import pytest
import pandas as pd
from scripts import generate_ped

# create mock linker dataframe to generate PED for
@pytest.fixture
def passing_linker():
    return pd.DataFrame({
        "Sample": ["Sample1", "Sample2"],
        "Sex": ["M", "F"]
    })

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

""" Can't call multiple fixtures directly in an elegant way for 
test_add_ped_columns so here I create a fixture that just calls 
the ped_format() function to return the dataframe and use as 
input for test_encoded_sex.

Reference: https://github.com/pytest-dev/pytest/issues/3950 """

@pytest.fixture(name="ped_format_indirect")
def ped_format_indirect():
    return ped_format() 

# testing functions

# check_if_ped
@pytest.mark.parametrize(
    "test_in, exp_out", [("testfile.ped", True), ("testfile.tsv", False)],
)
def test_check_if_ped(test_in, exp_out):
    test_out = generate_ped.check_if_ped(test_in)
    assert exp_out == test_out

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

# add_ped_columns
def test_add_ped_columns(passing_linker):
    ped_out = generate_ped.add_ped_columns(passing_linker)
    exp_ped_out = ped_format()
    pd.testing.assert_frame_equal(ped_out, exp_ped_out)

# encode_sex
def test_encoded_sex(ped_format_indirect):
    encoded_out = generate_ped.encode_sex(ped_format_indirect)
    ped_format_example = ped_format()
    exp_encode_out = ped_format_example.replace({'Sex': {'M':1, 'F':2}})
    pd.testing.assert_frame_equal(encoded_out, exp_encode_out)
