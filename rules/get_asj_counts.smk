rule asj_nummarkers:
    input:
        os.path.join(tmpdir, "{sample}", "{sample}.genotypes.vcf")
    output:
        asj = os.path.join(tmpdir, "{sample}", "{sample}.asj.genotypes.txt"),
        gt = os.path.join(tmpdir, "{sample}", "{sample}.asj.tsv"),
        nonref = os.path.join(tmpdir, "{sample}", "{sample}.nonrefasjmarkers.txt"),
        nummarkers = os.path.join(tmpdir, "{sample}", "{sample}.asj.nummarkers.txt")
    params:
        include = config["asj_markers"]
    shell:
        r"""
        for i in `cat {params.include}`; do grep -w $i {input}; done > {output.asj}

        grep -v "#" {output.asj} | cut -f 3,10 | grep -v '\./\.' > {output.gt}
 
        echo -e Sample"\t"num_nonref_asj_markers > {output.nonref}
        echo -e {wildcards.sample}"\t"`grep -vc "0\/0" {output.gt}` >> {output.nonref}

        echo -e {wildcards.sample}"\t"`wc -l {output.gt} | cut -d ' ' -f 1` > {output.nummarkers}
        """
