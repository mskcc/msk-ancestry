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

localrules: merge_and_plot, create_pop_file, sample_admixture_output, cleanup, merge_nummarkers, merge_numnonrefasj, ash_nummarkers

rule cleanup:
    input:
        os.path.join(outdir,"ancestry_results.txt"),
        os.path.join(outdir,"QC", "num_markers_genotyped.txt"),
        expand(os.path.join(outdir, "individual_admixture_results", "{sample}.nummarkers.txt"), sample=samples),
        expand(os.path.join(outdir, "individual_ASJ_results", "{sample}.nummarkers.txt"), sample=samples)
    params:
        tmpdir = tmpdir
    shell:
        """
        if [ -d "{params.tmpdir}" ]; then
            rm -rf {params.tmpdir}
        fi
        """

rule merge_nummarkers:
    input:
        admix=expand(os.path.join(outdir, "individual_admixture_results", "{sample}.nummarkers.txt"), sample=samples),
        asj=expand(os.path.join(outdir, "individual_ASJ_results", "{sample}.nummarkers.txt"), sample=samples)
    output:
        os.path.join(outdir,"QC", "num_markers_genotyped.txt")
    shell:
        """
        cat {input.admix} > {output}.1
        cat {input.asj} > {output}.2
        echo -e "SAMPLE\tNUM_ADMIXTURE_MARKERS\tNUM_ASJ_MARKERS" > {output}
        paste {output}.1 {output}.2 | cut -f 1,2,4 >> {output}
        rm {output}.1 {output}.2
        """

rule merge_numnonrefasj:
    input:
        expand(os.path.join(outdir, "individual_ASJ_results", "{sample}.nonrefasjmarkers.txt"), sample=samples)
    output:
        os.path.join(outdir,"individual_ASJ_results", "num_nonref_ASJ_markers.txt")
    shell:
        """
        cat {input} > {output}
        """
rule merge_and_plot:
    input:
        admixture = expand(os.path.join(outdir, "individual_admixture_results", "{sample}.admixture_results.txt"), sample=samples),
        asj = os.path.join(outdir,"individual_ASJ_results", "num_nonref_ASJ_markers.txt")
    output:
        txt = os.path.join(outdir,"ancestry_results.txt"),
        pdf = os.path.join(outdir,"admixture_results.pdf")
    params:
        asjcutoff = config["asj_marker_cutoff"],
        admixcutoff = config["admix_frac_cutoff"]
    conda:
        "envs/rplot.yaml"
    script:
        "scripts/plot_admixture.R"


include: "rules/genotype_markers.smk"
include: "rules/get_asj_counts.smk"
include: "rules/merge_with_1KGP.smk"
include: "rules/run_admixture.smk"
