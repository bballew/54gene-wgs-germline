Usage and Running the pipeline
==============================

Below are descriptions and usage options for the various config parameters specified in ``config.yaml``.

+---------------+-----------+----------------+---------------------------+
| Parameter     |  Required |  Default value |         Description       |
+===============+===========+================+===========================+
| sampleFile    |     Y     |      NA        |  Manifest file with IDs   |
+---------------+-----------+----------------+---------------------------+
| intervalsFile |     Y     |      NA        | File with interval names  |
|               |           |                |  and file paths           |
+---------------+-----------+----------------+---------------------------+ 
| jobs          |     Y     |     1000       | Jobs to submit to cluster |
+---------------+-----------+----------------+---------------------------+
| sexLinker     |     Y     |      NA        | File with reported sex of |
|               |           |                | sexes of each sample ID   |     
+---------------+-----------+----------------+---------------------------+
| tempDir       |     Y     |     temp/      | Location of temp directory|
+---------------+-----------+----------------+---------------------------+
| runType       |     Y     |    full: yes   | Specify run mode to use   |
+---------------+-----------+----------------+---------------------------+
| global_vars   |     N     |  As specified  | Set global java options   |
+---------------+-----------+----------------+---------------------------+
| cluster_mode  |     N     |  As specified  | Used to submit jobs       |
+---------------+-----------+----------------+---------------------------+
| default_queue |     Y     |     "big"      | Name of your default      |
|               |           |                | cluster partition/queue   |
+---------------+-----------+----------------+---------------------------+
| compute_queue |     Y     |     "big"      | Name of queue/partition   |
|               |           |                | with compute-heavy nodes  |
+---------------+-----------+----------------+---------------------------+
| memory_queue  |     Y     |     "big"      | Name of queue/parition    |
|               |           |                | with nodes for memory-    |
|               |           |                | intensive jobs            |
+---------------+-----------+----------------+---------------------------+
| bed           |     N     | Homo_sapiens_  | Split genome by chromosome|
|               |           | assembly38.bed | when calling variants     |
+---------------+-----------+----------------+---------------------------+
| max_concurent |     Y     |      40        |Max concurrent running jobs|
+---------------+-----------+----------------+---------------------------+
| max_het_ratio |     Y     |      2.5       | Max het/hom ratio to allow|
|               |           |                |                           |
+---------------+-----------+----------------+---------------------------+
| min_avg_depth |     Y     |      20        | Minimum depth required for|
|               |           |                | sample                    |
+---------------+-----------+----------------+---------------------------+
| max_contam    |     Y     |      0.03      | Max % of contam allowed   |
+---------------+-----------+----------------+---------------------------+



Resource Allocation
-------------------
This pipeline was originally developed to be run on an Unix-like HPC system for scalable and efficient analyses. As a result, there are several configuration options available for allocating resources and running jobs within the workflow.


Within the ``config/config.yaml`` there are several options for allocating memory, threads, and setting Java-specific variables (for GATK suite of tools) for the various tools used in the workflow.

- To modify the number of maximum jobs to run, you can specify a number for ``jobs``
- Set a maximum number of jobs to run concurrently if you have bandwith constraints using the ``max_concurrent`` variable
- Specify which paritions/queues to submit your jobs to, depending on your HPC using the ``default_queue``, ``compute_queue`` and ``memory_queue`` variables (some jobs such as alignment, require more memory and can use nodes with more memory available when specified with these variables)
- Specify the number of threads and memroy in MB for each tool, where available using the ``threads`` and ``memory`` variables
- Specify the space to allocate for Java lass metadata using the ``global_vars`` variable 
  