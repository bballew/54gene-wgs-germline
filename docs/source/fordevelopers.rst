For developers
==============

Contribute back
---------------

For developers, if you would like to contribute back to this repository, consider forking the original repo and creating a pull request:

1. `Fork <https://help.github.com/en/articles/fork-a-repo>`_ the original repo to a personal or lab account.
2. `Clone <https://help.github.com/en/articles/cloning-a-repository>`_ the fork to your local system, to a different place than where you ran your analysis.
3. Copy the modified files from your analysis to the clone of your fork, e.g., ``cp -r workflow path/to/fork``. (Make sure to not accidentally copy config file contents or sample sheets. Instead, manually update the example config files if necessary)
4. Commit and push your changes to your fork.
5. Create a `pull request <https://help.github.com/en/articles/creating-a-pull-request>`_ against the original repository.

Obtain updates from upstream
----------------------------
Whenever you want to synchronize your workflow copy with new developments from upstream, do the following.

1. Once, register the upstream repository in your local copy: ``git remote add -f upstream git@github.com:snakemake-workflows/54gene-wgs-germline.git`` or ``git remote add -f upstream https://github.com/snakemake-workflows/54gene-wgs-germline.git`` if you do not have setup ssh keys.
2. Update the upstream version: ``git fetch upstream``
3. Create a diff with the current version: ``git diff HEAD upstream/master workflow > upstream-changes.diff``
4. Investigate the changes: ``vim upstream-changes.diff``
5. Apply the modified diff via: ``git apply upstream-changes.diff``
6. Carefully check whether you need to update the config files: ``git diff HEAD upstream/master config``. If so, do it manually, and only where necessary, since you would otherwise likely overwrite your settings and samples.
