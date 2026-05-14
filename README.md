# msk-ancestry

This is a snakemake workflow for ancestry inference from MSK-IMPACT data. It takes BAM file(s) to genotype ancestry informative markers selected for MSK-IMPACT panels, uses 1000 genomes populations as reference and runs supervised ADMIXTURE to estimate the ancestral proportions of European (EUR), African (AFR), Native American (NAM), East Asian (EAS) and South Asian (SAS) for the samples. Additionally, it also genotypes Ashkenazi Jewish (ASJ) ancestry informative markers to infer ASJ ancestry.
The output file contains ancestry proportions for EUR, AFR, NAM, EAS, SAS; the ASJ ancestry inference and final ancestry labels assigned to each sample. 

### Requirements:
- python3
  - pandas
  - numpy
- conda
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [snakemake lsf profile](https://github.com/Snakemake-Profiles/lsf) (if running on LSF cluster)
OR
- [snakemake executor plugin: slurm](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html) (if running on SLURM cluster)
### Example commands:

#### Local execution 
`snakemake --use-conda`

#### LSF cluster:
`snakemake --use-conda --profile lsf`

#### SLURM cluster:
`snakemake --use-conda \
  --executor slurm \
  --jobs 50 \
  --default-resources \
    slurm_partition=<SLURM_PARTITION_NAME> \
    slurm_account=<SLURM_ACCOUNT_NAME> \
    runtime=100
`

Notes on SLURM execution:
- `--jobs` controls the maximum number of SLURM jobs submitted concurrently
- `runtime` specifies the requested walltime (in minutes) for each job.
- `slurm_partition` specifies the SLURM partition to submit jobs to. Since most jobs in this pipeline complete within minutes, a short or fast partition is typically appropriate. 
- `slurm_account` is your SLURM billing account and can be found using: `sacctmgr show user $USER format=User,DefaultAccount,Account`

This expects the Snakemake file in the current working directory.
You can modify the ***config.yaml*** file to change the input file or output directory. By default, it uses the ***metadata/samples.tsv*** file as input.

### NOTE:
- Currently the config file points to reference FASTA file in IRIS. This file will have to copied and the config file will have to be modified if you want to run this on a machine that cannot access these files.
- The markers used in this workflow are chosen for MSK-IMPACT data. You would have to create new set of marker files to run this on data from different sequencing panels. See instructions [here](https://github.com/mskcc/msk-ancestry/blob/main/MarkerSelection.md).
- The reference FASTA file and marker files are for GRCh37 reference genome. You would have to create new set of marker files and provide a different reference FASTA file to run this on BAM files aligned to a different reference genome version.
- The sample name can only have alphanumeric characters, underscore (\_) or hyphen (-). If it has any other characters, the workflow will fail.  

## Workflow diagram:
<img src="https://github.com/mskcc/msk-ancestry/blob/main/misc/workflow_diagram.png" width="40%" height="40%">
