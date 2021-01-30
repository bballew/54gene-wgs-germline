import os
import shutil
import subprocess as sp
import sys
from pathlib import Path, PurePosixPath
from tempfile import TemporaryDirectory

sys.path.insert(0, os.path.dirname(__file__))

import common  # noqa: E402


def test_fastqc():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/fastqc/data")
        expected_path = PurePosixPath(".tests/unit/fastqc/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print(
            "results/fastqc/NA12878_readgroup1_r1_fastqc.html results/fastqc/NA12878_readgroup1_r1_fastqc.zip results/fastqc/NA12878_readgroup1_r2_fastqc.html results/fastqc/NA12878_readgroup1_r2_fastqc.zip",
            file=sys.stderr,
        )

        # Run the test job.
        sp.check_output(
            [
                "python",
                "-m",
                "snakemake",
                "results/fastqc/NA12878_readgroup1_r1_fastqc.html results/fastqc/NA12878_readgroup1_r1_fastqc.zip results/fastqc/NA12878_readgroup1_r2_fastqc.html results/fastqc/NA12878_readgroup1_r2_fastqc.zip",
                "-F",
                "-j1",
                "--keep-target-files",
                "--directory",
                workdir,
            ]
        )

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file),
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
