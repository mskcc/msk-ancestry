import pandas as pd
import os

configfile: "config.yaml"
snakemakedir =  workflow.basedir
cwd = os.getcwd()

input = pd.read_table(config["input"], dtype="str").drop_duplicates()
bams = pd.read_table(config["input"], dtype=str).drop_duplicates().set_index("SAMPLE")["BAM"].to_dict()
samples = list(bams.keys())
outdir = os.path.abspath(config["outdir"])
## create temp directory for the run ##
tmpdir = os.path.join(outdir,"tmp_snakemake_out")

reference = config["reference"]
kgbed = config["kgbed"]
markers_txt = config["markers_txt"]
markers_vcf = config["markers_vcf"]
kgpops = config["kgpops"]
numpops = len(pd.read_table(kgpops, dtype="str")["SuperPop"].drop_duplicates())

wildcard_constraints:
    sample="[a-zA-Z0-9_\-]+"

localrules: merge_and_plot, create_pop_file, sample_admixture_output, cleanup

rule cleanup:
    input:
        os.path.join(outdir,"admixture_results.txt")
    params:
        tmpdir = tmpdir
    shell:
        """
        if [ -d "{params.tmpdir}" ]; then
            rm -rf {params.tmpdir}
        fi
        """

rule merge_and_plot:
    input:
        expand(os.path.join(outdir, "individual_results", "{sample}.admixture_results.txt"), sample=samples)
    output:
        txt = os.path.join(outdir,"admixture_results.txt"),
        pdf = os.path.join(outdir,"admixture_results.pdf")
    conda:
        "envs/rplot.yaml"
    script:
        "scripts/plot_admixture.R"


include: "rules/genotype_markers.smk"
include: "rules/merge_with_1KGP.smk"
include: "rules/run_admixture.smk"
