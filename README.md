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


### Step 2

Edit `config.py`. All the variables in the file
are self explanatory. 

The `SAMPLES` variable refers to the filename of fastqs leaving out the extension
whihch is assumed to be `.fq`.

`SRC_DIR`: Path to the `scripts/` directory(absolute)

`RAWDATA_DIR`: Where should the fasts be read fromi.(`ROOT_DIR` is redundant)

`STAR_INDEX`: Directory location of STAR index

`_BED`: Location to BED files, required for annotation


## Step 3 

On hpc-cmb login node:

`bash submitall.sh`

Error log will be printed on STDOUT indicating if jobs are running or not.

