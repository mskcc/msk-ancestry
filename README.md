# msk-ancestry

This is a snakemake workflow for ancestry inference from MSK-IMPACT data. It uses 1000 genomes populations as reference and runs supervised ADMIXTURE to estimate the ancestral proportions of European (EUR), African (AFR), Native American (NAM), East Asian (EAS) and South Asian (SAS) for the samples. Additionally, it also genotypes Ashkenazi Jewish (ASJ) ancestry informative markers to infer ASJ ancestry.
The output file contains ancestry proportions for EUR, AFR, NAM, EAS, SAS; the ASJ ancestry inference and final ancestry labels assigned to each sample. 
Currently it is only tested on the juno cluster. 

### Requirements:
- python3
  - pandas
  - numpy
- conda
- [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
- [snakemake lsf profile](https://github.com/Snakemake-Profiles/lsf) (if running on LSF cluster)

### Example command:
`snakemake --use-conda`
If running on LSF cluster: `snakemake --use-conda --profile lsf`

This expects the Snakemake file in the current working directory.
You can modify the ***config.yaml*** file to change the input file or output directory. By default, it uses the ***metadata/samples.tsv*** file as input.

### NOTE:
- Currently the config file points to reference FASTA file in juno. This file will have to copied and the config file will have to be modified if you want to run this on a machine that cannot access these files.
- The markers used in this workflow are chosen for MSK-IMPACT data. You would have to create new set of marker files to run this on data from different sequencing panels.
- The sample name can only have alphanumeric characters, underscore (\_) or hyphen (-). If it has any other characters, the workflow will fail.  

## Workflow diagram:
<img src="https://github.com/mskcc/msk-ancestry/blob/main/misc/workflow_diagram.png" width="40%" height="40%">
