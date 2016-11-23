# iclip-seq-snakemake
Snakemake based pipeline for analysing iCLIP-Seq data
## To Run

In order to make full use of the pipeline, it should be run on cluster.
It submits independent qsub job for all steps which are independent.
If there are 3 replicates the pipeline will run paralley for all 3, using 
3 qsub jobs.

### Step 1

Install dependencies:

`conda create -n clipseq -f environment.yml`


### Step 2a

Edit `config.py`. All the variables in the file
are self explanatory.

The `SAMPLES` variable refers to the filename of fastqs leaving out the extension
which is assumed to be `.fq`.

`SRC_DIR`: Path to the `scripts/` directory(absolute)

`RAWDATA_DIR`: Where should the fasts be read fromi.(`ROOT_DIR` is redundant)

`STAR_INDEX`: Directory location of STAR index

`_BED`: Location to BED files, required for annotation

The `_OTHER` fields are relevant only if you have data from some 'other' species. For example, if you have human
and mouse iCLIP data, these fields are used to do a lot of liftover operations. If you have single specie data, you
can leave the `_OTHER` fields as they are.

## Step 2b
Edit `jobscript.sh` so that it has the correct [`PATH` variable](https://github.com/saketkc/iclip-seq-snakemake/blob/master/jobscript.sh#L5). Make sure it includes the conda environment as the first.
In my case it is `home/cmb-panasas2/skchoudh/software_frozen/anaconda2/envs/clipseq/bin`. The reason to do this hard coding of `PATH` is because some hpc nodes fail to run the job because they can't detect the environment often.

## Step 3 

On hpc-cmb login node:

`bash submitall.sh`

Error log will be printed on STDOUT indicating if jobs are running or not.

