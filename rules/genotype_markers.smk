def get_mem_mb(wildcards, attempt):
    return attempt *  4000

rule genotype:
    input:
        pileup = os.path.join(tmpdir, "{sample}", "{sample}.pileup.txt"),
        reference = reference,
        markers_txt = markers_txt
    conda:
        "../envs/python310.yaml"
    output:
        os.path.join(tmpdir, "{sample}", "{sample}.genotypes.vcf")
    params:
        snakemakedir = snakemakedir
    resources:
        mem_mb=2000
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
        mem_mb=lambda wc, attempt: 4000 * attempt
    retries: 3
    threads: 1
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
