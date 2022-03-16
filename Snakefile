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

localrules: merge_and_plot, create_pop_file, sample_admixture_output, cleanup, create_qc_files, merge_numnonrefasj, asj_nummarkers, create_final_ancestry_files, merge_ancestry_files

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

rule create_qc_files:
    input:
        admix=expand(os.path.join(outdir, "individual_admixture_results", "{sample}.nummarkers.txt"), sample=samples),
        asj=expand(os.path.join(outdir, "individual_ASJ_results", "{sample}.nummarkers.txt"), sample=samples)
    output:
        os.path.join(outdir,"QC", "num_markers_genotyped.txt")
    shell:
        """
        cat {input.admix} > {output}.1
        cat {input.asj} > {output}.2
        echo -e "Sample\tnum_admixture_markers\tnum_ASJ_markers" > {output}
        paste {output}.1 {output}.2 | cut -f 1,2,4 >> {output}
        rm {output}.1 {output}.2
        """

rule create_final_ancestry_files:
    input:
        admixture = os.path.join(outdir, "individual_admixture_results", "{sample}.admixture_results.txt"),
        asj = os.path.join(outdir, "individual_ASJ_results", "{sample}.nonrefasjmarkers.txt")
    output:
        out = os.path.join(outdir,"final_ancestry_results", "{sample}.ancestry_results.txt")
    params:
        asjcutoff = config["asj_marker_cutoff"],
        admixcutoff = config["admix_frac_cutoff"]
    run:
        print(output)
        adm=pd.read_table(input.admixture, sep="\t")
        asj=pd.read_table(input.asj, sep="\t")
        final=adm.merge(asj, on="Sample")
        final["maxfrac"] = final.drop(columns=["Sample","num_nonref_asj_markers"]).max(axis=1)
        final["admixture_label"] = final.drop(columns=["Sample","num_nonref_asj_markers"]).idxmax(axis=1)
        final["admixture_label"] = np.where(final["maxfrac"] < params.admixcutoff, "ADM", final["admixture_label"])
        final["ASJ"] = np.where(final['num_nonref_asj_markers'] >= params.admixcutoff, "yes", "no")
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
