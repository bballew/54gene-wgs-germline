Tests
=====

Unit tests
----------

Unit tests for the python modules used in this workflow can be found in ``workflow/scripts/tests`` and run using Pytest which is included in the conda run-time environment for this pipeline.

To run all available unit tests::

    conda activate 54gene-wgs-germline
    pytest -s workflow/scripts/tests/*.py


Pipeline/Integration tests
--------------------------
To test the core pipeline, we provide a small test dataset and instructions on how to use this dataset available in a repository `here <https://gitlab.com/data-analysis5/dna-sequencing/54gene-wgs-test-data>`_.

To summarize, this test dataset contains a small region of chromosome 21 from the NA12878 platinum reference genome. The above repository contains all necessary inputs (configuration file, manifest, intervals, sex_linker files) required to run the pipeline in all three run-modes. The README provides instructions on how to use these files to execute a test using the 54gene-wgs-germline pipeline.

Snakemake unit tests
--------------------
*In development (TBD)*

CI/CD
-----
The aforementioned python unit tests and integration tests (in all three run-modes) are run as part of the `Gitlab Continuous Integration (CI) <https://docs.gitlab.com/ee/ci/>`_ pipeline for this codebase. You can find the status of the CI pipeline on the main repository page.

*Note*: The test suite and CI pipeline are still a work in progress.
