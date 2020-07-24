# vim: set ft=python:


def get_sample_reads(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    r1 = list(df[df["coassembly"] == sample]["r1"])[0].split(",")
    r2 = list(df[df["coassembly"] == sample]["r2"])[0].split(",")
    dict = {"r1": r1, "r2": r2}
    return dict


rule rename_fasta:
    input:
        MAG=join(DATA_DIR, binning_dir, "singlerun/{sample}/refined_bins_50_10/{file}.fa"),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    params:
        name="{sample}_metawrap_refined_{file}",
    shell:
        """
        scripts/rename_multifasta_prefix.py -f {input.MAG} -p {params.name} > {output}
        """


rule prokka:
    input:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
        prefix="{sample}_metawrap_refined_{file}",
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """


rule rename_fasta_coas:
    input:
        MAG=join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10/{file}.fa"),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    params:
        name="{sample}_metawrap_refined_{file}",
    shell:
        """
        scripts/rename_multifasta_prefix.py -f {input.MAG} -p {params.name} > {output}
        """


rule prokka_coas:
    input:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    singularity:
        "docker://staphb/prokka:1.14.5"
    params:
        out_prokka=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/"),
        prefix="{sample}_metawrap_refined_{file}",
    shell:
        """
        prokka {input} --kingdom Bacteria --outdir {params.out_prokka} \
        --prefix {params.prefix} --force --locustag {params.prefix} --cpus {threads}
        """


rule cmseq:
    input:
        r1=join(DATA_DIR, "00_preprocessing/processed/singlerun/{sample}_1.fastq"),
        r2=join(DATA_DIR, "00_preprocessing/processed/singlerun/{sample}_2.fastq"),
        MAG=join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
        prokka=join(
            DATA_DIR,
            binning_analyses,
            "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.cmseq.csv"),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/Containers:cmseq"
    params:
        name=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}"),
        dir=join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
    shell:
        """
        scripts/cmseq.sh -t {threads} -i {input.r1} -n {input.r2} -r {input.MAG} -g {input.prokka} -o {params.name}
        """


rule cmseq_coas:
    input:
        MAG=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.fa",
        ),
        prokka=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.gff",
        ),
    output:
        cmseq=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}.cmseq.csv",
        ),
        done=join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_{file}_done.txt",
        ),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/Containers:cmseq"
    params:
        r1=lambda wildcards: get_sample_reads(wildcards.sample)["r1"],
        name=join(
            DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/{sample}_metawrap_refined_{file}"
        ),
    shell:
        """
        rm -f {output.cmseq}
        rm -f {output.done}

        for i in {params.r1}; do run=$(basename ${{i}} _1.fastq); scripts/cmseq.sh\
        -t {threads} -i ${{i}} -n ${{i%%_1.fastq}}_2.fastq -r {input.MAG} -g {input.prokka}\
        -o {params.name}_${{run}}; awk 'NR==2' {params.name}_${{run}}.cmseq.csv >> {output.cmseq};done

        touch {output.done}
        """


rule append_cmseq:
    input:
        lambda wildcards: expand(rules.cmseq.output, file=aggregate_bins(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/sum_cmseq_{sample}.csv"),
    params:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/refined_bins_50_10_renamed/"),
    run:
        if os.path.exists(str(output)):
            os.remove(str(output))
        else:
            print("Aggregating single run strain heterogeneity results")
        path = str(params)
        file1 = open(str(output), "w")
        for filename in glob.glob(os.path.join(path, "*.csv")):
            with open(filename, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    lines_sub = lines[1:2][0]
                    L = join(filename, "\t", lines_sub)
                    file1.writelines(L)
        file1.close()



rule append_cmseq_coas:
    input:
        lambda wildcards: expand(rules.cmseq_coas.output, file=aggregate_bins_coas(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/sum_cmseq_{sample}.csv"),
    params:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/refined_bins_50_10_renamed/"),
    run:
        if os.path.exists(str(output)):
            os.remove(str(output))
        else:
            print("Aggregating coassembly strain heterogeneity results")
        path = str(params)
        file1 = open(str(output), "w")
        for filename in glob.glob(os.path.join(path, "*SRR*.csv")):
            with open(filename, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    lines_sub = lines[1:2][0]
                    L = join(filename, "\t", lines_sub)
                    file1.writelines(L)
        file1.close()



rule aggregate_cmseq:
    input:
        sr=expand(join(DATA_DIR, binning_analyses, "singlerun/cmseq/{sample}/sum_cmseq_{sample}.csv"), sample=RUN),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/cmseq/summ_cmseq_all.csv"),
    shell:
        """
        cat {input}>{output}         
        """


rule aggregate_cmseq_coas:
    input:
        sr=expand(join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/{sample}/sum_cmseq_{sample}.csv"), sample=COAS),
    output:
        coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/summ_cmseq_all.csv"),
    shell:
        """
        cat {input}>{output}
        """


rule plot_cmseq:
    input:
        coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/cmseq/summ_cmseq_all.csv"),
        sr=join(DATA_DIR, binning_analyses, "singlerun/cmseq/summ_cmseq_all.csv"),
    output:
        join(DATA_DIR, "figures/cmseq_plot.png"),
    singularity:
        "shub://sskashaf/Containers:r"
    shell:
        """
        Rscript scripts/plotting/plot_cmseq.R {input.sr} {input.coas}
        """
