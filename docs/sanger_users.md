
On the Sanger Farm the user group must be specified to submit jobs. Since Nextflow does this for you, you need to set a default user group in your `~/.bashrc` like this so Farm understands which group you're running your jobs under

```bash
export LSB_DEFAULT_USERGROUP=your-unix-group
```

For the first time running the pipeline, remember to `source ~/.bashrc`. You don't have to do this again next time you log in.

After the pipeline successfully runs and the expected outputs have been checked, please delete the `work` and the `.nextflow` directory, as well as the `.nextflow.log*` files. This is because Farm has a limit on the number of files that your group can have, and this pipeline creates a huge number of them.

```bash
rm -r ${workdir}/work ${workdir}/.nextflow
rm ${workdir}/.nextflow.log*
```