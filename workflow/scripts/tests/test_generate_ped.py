import pytest
from scripts import generate_ped


@pytest.mark.parametrize(
    "test_in, exp_out", [("testfile.ped", True), ("testfile.tsv", False), ("testfile.PED", True)],
)
def test_check_if_ped(test_in, exp_out):
    test_out = generate_ped.check_if_ped(test_in)
    assert exp_out == test_out


# read_in_linker
# check_column_number
# check_column_headers
# add_ped_columns
# encode_sex
