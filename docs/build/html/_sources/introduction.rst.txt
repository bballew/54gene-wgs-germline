Introduction
===============

This workflow was designed by the Genomics & Data Science team at 54gene and is used to analyze germline whole-genome sequencing data in the form of either FASTQs or gVCFs. This pipeline emit a joint-called multi-sample VCF. It is currently optimized to be run on HPC infrastructure and was developed on AWS' ParallelCluster.

Features:

- Read filtering and trimming
- Read alignment
- Variant calling and filtering
- Joint-genotyping
- Check sex and relatedness
- Generate a multiQC report 

To install the latest release, type::
    git clone <link here>

See the :doc:`installation` and :doc:`usage` for details.
