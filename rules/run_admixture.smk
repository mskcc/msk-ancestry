rule create_pop_file:
    input:
        os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.fam")
    output:
        os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.pop")
    params:
        kgpops = kgpops
    shell:
        """
        while read -r fam sample remain; do
            pop=`grep "^{wildcards.sample}" {params.kgpops}` | cut -f 5`
            if [[ $pop == "" ]]; then
                echo "-"
            else
                echo $pop
            fi
        done < {input} | sed 's/AMR/NAM/' > {output}
        """

rule run_supervised_admixture:
    input:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.bed"),
        pop = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.pop")
    output:
        os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref."+str(numpops)+".Q")
    params:
        numpops = numpops
    conda:
        "../envs/admixture.yaml"
    shell:
        "admixture --supervised {input.bed} {params.numpops}"
