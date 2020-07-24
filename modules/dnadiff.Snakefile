# vim: set ft=python:


rule mash_dist:
    input:
        bins=join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
        db=config["db"]["mash_RefSeq"],
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/mashdist/{sample}_metawrap_refined_{i}.tab"),
    threads: workflow.cores
    singularity:
       # "shub://sskashaf/Containers:isolatescompare"
        "docker://quay.io/biocontainers/mash:2.2.1--h3d38be6_0"
    shell:
        "mash dist -p {threads} {input.db} {input.bins} > {output}"


rule best_mash:
    input:
        mashdist=join(
            DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/mashdist/{sample}_metawrap_refined_{i}.tab"
        ),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/best_mash/{sample}_metawrap_refined_{i}.tab"),
    threads: workflow.cores
#    singularity:
 #       "shub://sskashaf/Containers:isolatescompare"
    shell:
        """
        sort -gk3 {input.mashdist}|sed -n 1p >{output}
        """


rule dnadiff:
    input:
        bestmash=join(
            DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/best_mash/{sample}_metawrap_refined_{i}.tab"
        ),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/{sample}_metawrap_refined_{i}.report"),
    params:
        outdir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff"),
        bins="{sample}_metawrap_refined_{i}",
    singularity:
     #   "shub://sskashaf/Containers:isolatescompare"
         "docker://quay.io/biocontainers/mummer:3.23--pl526_7"
    shell:
        """
        while read col1 col2 rem
          do
            echo 'dnadiff ${{col1}} ${{col2}} -p ${{col1%%.fasta}}_${{col2%%.fa}}_'
            dnadiff ${{col1}} ${{col2}} -p {params.outdir}/{params.bins}
          done < {input.bestmash}
        """


rule parse_dnadiff:
    input:
        dnadiff=join(
            DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/{sample}_metawrap_refined_{i}.report"
        ),
    output:
        join(
            DATA_DIR,
            binning_analyses,
            "singlerun_coassembly/MAG_RefSeq/dnadiff/{sample}_metawrap_refined_{i}_dnadiff_parsed.tsv",
        ),
    run:
        outfile = str(output)
        f = open(input.dnadiff)
        data = f.read()
        first_line = data.split("\n", 1)[0]
        a = first_line.split(" ")
        ref = a[0]
        quer = a[1]
        with open(outfile, "w") as outf:
            path_dna = input.dnadiff
            base = os.path.basename(path_dna)
            base = base.split(".report")[0]
            with open(path_dna) as f:
                for line in f:
                    if "TotalBases" in line:
                        cols = line.split()
                        lenref = int(cols[1])
                        lenquer = int(cols[2])
                    if "AlignedBases" in line:
                        cols = line.split()
                        aliref = cols[1].split("(")[-1].split("%")[0]
                        alique = cols[2].split("(")[-1].split("%")[0]
                    if "AvgIdentity" in line:
                        cols = line.split()
                        ident = float(cols[1])
            line = "%s\t%s\t%i\t%.2f\t%i\t%.2f\t%.2f" % (ref, quer, lenref, float(aliref), lenquer, float(alique), float(ident))
            outf.writelines(line, "\n")



rule aggregate_dnadiff:
    input:
        lambda wildcards: expand(rules.parse_dnadiff.output, i=aggregate_bins(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/singlerun/dnadiff_{sample}_summary.tsv"),
    shell:
        """
        cat {input}>{output}
        """


rule aggregate_dnadiff_coas:
    input:
        lambda wildcards: expand(rules.parse_dnadiff.output, i=aggregate_bins_coas(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/coassembly/dnadiff_{sample}_summary.tsv"),
    shell:
        """
        cat {input}>{output}
        """


rule summarize_dnadiff:
    input:
        expand(
            join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/singlerun/dnadiff_{sample}_summary.tsv"),
            sample=RUN,
        ),
        expand(
            join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/coassembly/dnadiff_{sample}_summary.tsv"),
            sample=COAS,
        ),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/dnadiff_summary.tsv"),
    shell:
        """
        cat {input}>{output}
        """


rule plot_dnadiff:
    input:
        summ=join(DATA_DIR, binning_analyses, "singlerun_coassembly/MAG_RefSeq/dnadiff/dnadiff_summary.tsv"),
        checkm_sr=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
        checkm_coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
    output:
        join(DATA_DIR, "figures/dnadiff.png"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:r"
    shell:
        """
        Rscript scripts/plotting/dnadiff_plot.R {input.checkm_sr} {input.checkm_coas} {input.summ}
        """
