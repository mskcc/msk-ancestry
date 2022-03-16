rule create_pop_file:
    input:
        fam = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.fam")
    output:
        pops = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.pop")
    params:
        kgpops = kgpops
    run:
        pops = pd.read_table(kgpops, dtype="str").set_index("#IID")["SuperPop"].to_dict()
        samples=pd.read_table(input.fam,dtype="str",sep=" ",header=None).iloc[:,1].to_list()
        o=open(output.pops, "w")
        for s in samples:
            if s in pops:
                o.write(pops[s]+"\n")
            else:
                o.write("-\n")
        o.close()

rule run_supervised_admixture:
    input:
        bed = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.bed"),
        pops = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.pop")
    output:
        os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref."+str(numpops)+".Q")
    params:
        numpops = numpops,
        tmpdir = tmpdir,
        cwd = os.getcwd()
    conda:
        "../envs/admixture.yaml"
    shell:
        """
        cd {params.tmpdir}/{wildcards.sample}
        admixture --supervised {input.bed} {params.numpops}
        cd {params.cwd}
        """

rule sample_admixture_output:
    input:
        admix = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref."+str(numpops)+".Q"),
        fam = os.path.join(tmpdir, "{sample}", "{sample}.w1kgpref.fam")
    output:
        os.path.join(outdir, "individual_admixture_results", "{sample}.admixture_results.txt")
    params:
        kgpops = kgpops
    run:
        admix = pd.read_table(input.admix, dtype="float",header=None,sep=" ")
        fam = pd.read_table(input.fam,dtype="str",sep=" ",header=None).iloc[:,1]
        pops = pd.read_table(kgpops, dtype="str").set_index("#IID")["SuperPop"].to_dict()
        df = pd.concat([fam,admix],axis=1)
        header=["Sample"]
        # for the populations in admixture results
        for i in range(1,len(df.columns)):
            # samples that are predominantly this population
            samples = df[df.iloc[:, i] > 0.99].iloc[:,0]
            for s in samples:
                if s in pops:
                    header.append(pops[s])
                    break
        df.columns=header
        final = df[df["Sample"] == wildcards.sample]
        final.to_csv(output[0], index=False, sep="\t")
