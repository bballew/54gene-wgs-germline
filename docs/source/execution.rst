Execution
=========

Deploying the pipeline
----------------------

With the ``config.yaml`` configured to your run-mode of choice with paths to the necessary manifest and input files, the workflow can be executed on any infrastructure using the ``snakemake`` command, supplied with further Snakemake command-line arguments (e.g. specifying a profile with ``--profile`` or ``--cluster`` to submit jobs to an HPC) depending on your environment.

Test your configuration by performing a dry-run::

    snakemake --use-conda -n

Execute the workflow locally via::

    snakemake --use-conda --cores $N

Execute the workflow on a cluster using something like::

    snakemake --use-conda --cluster sbatch --jobs 100


The pipeline will automatically create a subdirectory for logs in ``logs/`` and temporary workspace at the path specified for ``tempDir`` in the ``config.yaml``.

Wrapper scripts
---------------

We have provided two convenience scripts in the 54gene-wgs-germline repository to execute the workflow in a cluster environment: ``run.sh`` and ``wrapper.sh``.  You may customize these scripts for your needs, or run using a profile (e.g. `this <https://github.com/Snakemake-Profiles/slurm>`_ profile for a slurm job scheduler).

The ``wrapper.sh`` script embeds the ``snakemake`` command and other command-line flags to control submission of jobs to an HPC using the ``cluster_mode`` string pulled from the ``config.yaml``. This script also directs all stdout from Snakemake to a log file named ``WGS_${DATE}.out`` which will include the latest git tag and version of the pipeline, if cloned from our repository. For additional logging information, see :ref:`logging`.

This wrapper script can be edited to your needs and run using ``bash run.sh``.

.. _logging:

Logging
-------

All job-specific logs will be directed to a ``logs/`` subdirectory in the home analysis directory of the pipeline. This directory is automatically created for you upon execution of the pipeline. These log files will follow the naming structure of ``snakejob.<name_of_rule>.<job_number>``.

If you choose to use the ``wrapper.sh`` script provided and modified for your environment, a ``WGS_${DATE}.out`` log file containing all stdout from snakemake will also be available in the home analysis directory.
