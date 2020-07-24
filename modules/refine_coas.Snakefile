# vim: set ft=python:


rule bin_refinement_coas:
    input:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/done.txt"),
    output:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/done.txt"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:metawrap"
    params:
        metabat=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/metabat2_bins"),
        maxbin=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/maxbin2_bins"),
        concoct=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/concoct_bins"),
        outdir=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement"),
    threads: workflow.cores
    shell:
        """
        checkm data setRoot ~/checkm_database/
        rm -rf {params.outdir}
        rm -rf {output}
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.metabat} -B {params.maxbin} -C {params.concoct} -c 50 -x 10
         touch {output}
        """


checkpoint refine_bins_init_coas:
    input:
        rules.bin_refinement_coas.output,
    output:
        directory(join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10")),
    params:
        refined_folder=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins"),
        folder=join(DATA_DIR, binning_dir, "coassembly/{sample}/refined_bins_50_10/"),
    shell:
        """
        rm -rf {params.folder}
        mkdir -p {params.folder}
        cp {params.refined_folder}/*.fa {params.folder}
        """


def aggregate_bins_coas(wildcards):
    refine_out = checkpoints.refine_bins_init_coas.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(refine_out, "{bin}.fa")).bin


rule copy_bins_coas:
    input:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/{i}.fa"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/{sample}_metawrap_refined_{i}.fa"),
        join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    params:
        dir1=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/"),
        dir2=join(DATA_DIR, "all_bins/"),
        out1=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas/{sample}_metawrap_refined_{i}.fa"),
        out2=join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    shell:
        """
        mkdir -p {params.dir1}
        mkdir -p {params.dir2}
        cp {input} {params.out1}
        cp {input} {params.out2}
        """


rule copy_refined_coas:
    input:
        lambda wildcards: expand(rules.copy_bins_coas.output, i=aggregate_bins_coas(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/copied.txt"),
    shell:
        """
        touch {output}
        """


rule checkm_coas:
    input:
        expand(join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap_bin_refinement/copied.txt"), sample=COAS),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.tsv"),
    params:
        ext="fa",
        indir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas"),
        outdir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:metawrap"
    shell:
        """
        checkm data setRoot ~/checkm_database/
        rm -rf {params.outdir}
        checkm lineage_wf -t {threads} -x {params.ext} --tab_table -f {output} {params.indir} {params.outdir}
        """


rule parse_checkm_coas:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.tsv"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
    params:
        checkm=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics_tmp.csv"),
        checkm2=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics_tmp2.csv"),
    shell:
        """
        sed -i '1d' {input}
        cut -f1,12,13,14 {input} | tr '\t' ','>{params.checkm}
        sed 's/,/.fa,/' {params.checkm}>{params.checkm2}
        echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.checkm2} > {output}
        rm {params.checkm}
        rm {params.checkm2}
        """
