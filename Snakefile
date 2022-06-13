import pandas as pd
import os
import numpy as np

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

localrules: create_pop_file, sample_admixture_output, cleanup, asj_nummarkers, create_final_ancestry_files

rule cleanup:
    input:
        os.path.join(outdir,"ancestry_results.txt"),
        os.path.join(outdir,"QC", "num_markers_genotyped.txt"),
        expand(os.path.join(outdir,"final_ancestry_results", "{sample}.ancestry_results.txt"), sample=samples)
    params:
        tmpdir = tmpdir
    shell:
        """
        if [ -d "{params.tmpdir}" ]; then
            rm -rf {params.tmpdir}
        fi
        """

rule create_qc_files:
    input:
        admix=expand(os.path.join(tmpdir, "{sample}", "{sample}.admix.nummarkers.txt"), sample=samples),
        asj=expand(os.path.join(tmpdir, "{sample}", "{sample}.asj.nummarkers.txt"), sample=samples)
    output:
        out=os.path.join(outdir,"QC", "num_markers_genotyped.txt")
    run:
        dfs=[]
        for filename in input.admix:
            dfs.append(pd.read_table(filename, sep="\t", header=None))
        adm=pd.concat(dfs, ignore_index=True)
        adm.columns=["Sample","num_admixture_markers"]
        dfs=[]
        for filename in input.asj:
            dfs.append(pd.read_table(filename, sep="\t", header=None))
        asj=pd.concat(dfs, ignore_index=True)
        asj.columns=["Sample","num_ASJ_markers"]
        final=adm.merge(asj, on="Sample")
        final.to_csv(output.out, index=False, sep="\t")

rule create_final_ancestry_files:
    input:
        admixture = os.path.join(tmpdir, "{sample}", "{sample}.admixture_results.txt"),
        asj = os.path.join(tmpdir, "{sample}", "{sample}.nonrefasjmarkers.txt")
    output:
        out = os.path.join(outdir,"final_ancestry_results", "{sample}.ancestry_results.txt")
    params:
        asjcutoff = config["asj_marker_cutoff"],
        admixcutoff = config["admix_frac_cutoff"]
    run:
        adm=pd.read_table(input.admixture, sep="\t")
        asj=pd.read_table(input.asj, sep="\t")
        final=adm.merge(asj, on="Sample")
        final["maxfrac"] = final.drop(columns=["Sample","num_nonref_asj_markers"]).max(axis=1)
        final["admixture_label"] = final.drop(columns=["Sample","num_nonref_asj_markers"]).idxmax(axis=1)
        final["admixture_label"] = np.where(final["maxfrac"] < params.admixcutoff, "ADM", final["admixture_label"])
        final["ASJ"] = np.where(final['num_nonref_asj_markers'] >= params.asjcutoff, "yes", "no")
        final["ancestry_label"] = np.where((final["admixture_label"]=="EUR") & (final["ASJ"]=="yes"), "ASJ", final['admixture_label'])
        final=final.drop(columns=["maxfrac"])
        final.to_csv(output.out, index=False, sep="\t")

rule merge_ancestry_files:
    input:
        filename = expand(os.path.join(outdir,"final_ancestry_results", "{sample}.ancestry_results.txt"), sample=samples)
    output:
        out = os.path.join(outdir,"ancestry_results.txt")
    run:
        dfs=[]
        for filename in input.filename:
            dfs.append(pd.read_table(filename, sep="\t"))
            final=pd.concat(dfs, ignore_index=True)
            final.to_csv(output.out, index=False, sep="\t")

include: "rules/genotype_markers.smk"
include: "rules/get_asj_counts.smk"
include: "rules/merge_with_1KGP.smk"
include: "rules/run_admixture.smk"
