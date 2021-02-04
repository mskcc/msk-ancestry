rule vcf_to_plink:
    input:
        os.path.join(tmpdir, "{sample}", "{sample}.genotypes.vcf")
    output:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.bed"),
        ids = os.path.join(tmpdir, "{sample}", "{sample}.rsids.txt")
    conda:
        "../envs/plink.yaml"
    shell:
        """
        bedprefix=`echo {output.bed} | sed 's/.bed$//'`

        plink \
        --vcf {input} \
        --const-fid \
        --out $bedprefix \
        --make-bed \
        --geno

        cut -f 2 $bedprefix.bim > {output.ids}
        """
rule merge_with_1kgp:
    input:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.bed"),
        ids = os.path.join(tmpdir, "{sample}", "{sample}.rsids.txt")
    output:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.tmp.bed"),
        ids = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.tmp.prune.in")
    conda:
        "../envs/plink.yaml"
    params:
        ignore_fams = config["ignore_fams"],
        kgbed = kgbed
    shell:
        """
        inputbedprefix=`echo {input.bed} | sed 's/.bed$//'`
        outputbedprefix=`echo {output.bed} | sed 's/.bed$//'`

        plink \
        --bfile {params.kgbed} \
        --bmerge $inputbedprefix \
        --extract {input.ids} \
        --make-bed \
        --out $outputbedprefix \
        --indep-pairwise 500kb 50 0.2 \
        --remove {params.ignore_fams}
        """

rule keep_ld_pruned:
    input:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.tmp.bed"),
        ids = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.tmp.prune.in")
    output:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.bed")
    conda:
        "../envs/plink.yaml"
    shell:
        """
        inputbedprefix=`echo {input.bed} | sed 's/.bed$//'`
        outputbedprefix=`echo {output.bed} | sed 's/.bed$//'`

        plink \
        --bfile $inputbedprefix \
        --extract {input.ids} \
        --out $outputbedprefix \
        --make-bed
        """
