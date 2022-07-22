Installation
======================

This workflow was designed to use `conda <https://docs.conda.io/en/latest/>`_ for dependency management and utilizes  `Snakemake <https://snakemake.readthedocs.io/en/stable/>`_ as the workflow management system, ensuring reproducible and scalable analyses.

This installation guide is designed for Unix/Linux environments; however, this pipeline has been minimally tested on OSX as well.

Obtain a copy of this workflow
--------------------------------------

Clone this repository to your local system, into the place where you want to perform the data analysis::

    git clone git@gitlab.com:data-analysis5/54gene-wgs-germline.git

Install the run-time environment
----------------------------------------------

If needed, follow this `guide <https://docs.conda.io/en/latest/miniconda.html#installing>`_ to install Miniconda.

Once installed, create the run-time conda environment with minimal dependencies defined using the following command::
    
    conda env create -f environment.yaml

Activate the environment as follows::

    conda activate 54gene-wgs-germline

