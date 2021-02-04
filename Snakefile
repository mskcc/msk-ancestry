import pandas as pd
import os

configfile: "config.yaml"
snakemakedir =  workflow.basedir
cwd = os.getcwd()

input = pd.read_table(config["input"], dtype="str").drop_duplicates()
bams = pd.read_table(config["input"], dtype=str).drop_duplicates().set_index("SAMPLE")["BAM"].to_dict()
samples = bams.keys()

outdir = os.path.abspath(config["outdir"])
## create temp directory for the run ##
tmpdir = os.path.join(outdir,"tmp_snakemake_out")

reference = config["reference"]
kgbed = config["kgbed"]
markers_txt = config["markers_txt"]
markers_vcf = config["markers_vcf"]


rule final:
    input:
        expand(os.path.join(tmpdir,{sample}.genotypes.vcf), sample=samples)
# rule final:
#     input:
#         files: expand(tmpdir + "/{sample}.admixture_results.txt", sample=samples)
#     output:
#         outdir + "/admixture_results.tsv",
#         outdir + "/admixture_results.pdf"
#     conda:
#         "envs/plotting.yaml"
#     shell:
#         "scripts/merge_and_plot_results.R"


include: "rules/genotype_markers.smk"
