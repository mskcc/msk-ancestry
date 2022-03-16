def get_mem_mb(wildcards, attempt):
    return attempt * 2000

rule genotype:
    input:
        pileup = os.path.join(tmpdir, "{sample}", "{sample}.pileup.txt"),
        reference = reference,
        markers_txt = markers_txt
    output:
        os.path.join(tmpdir, "{sample}", "{sample}.genotypes.vcf")
    conda:
        "../envs/genotyping.yaml"
    params:
        snakemakedir = snakemakedir
    shell:
    	"{params.snakemakedir}/scripts/pileup_to_vcf.py \
    	-P {input.pileup} \
    	-M {input.markers_txt} \
    	-F {input.reference} \
    	-O {output} \
    	-N {wildcards.sample}"

rule pileup:
    input:
        bam = lambda wildcards: bams[wildcards.sample],
        markers_vcf = markers_vcf,
        reference = reference
    output:
        pileup = os.path.join(tmpdir, "{sample}", "{sample}.pileup.txt")
    conda:
        "../envs/gatk.yaml"
    resources:
        mem_mb=get_mem_mb
    shell:
        "gatk Pileup -R {input.reference} \
        -I {input.bam} \
        -L {input.markers_vcf} \
        -O {output.pileup} \
        -VS SILENT \
        -RF NotDuplicateReadFilter \
        -RF CigarContainsNoNOperator \
        -RF AmbiguousBaseReadFilter \
        -RF GoodCigarReadFilter \
        -RF MatchingBasesAndQualsReadFilter \
        --show-verbose"
