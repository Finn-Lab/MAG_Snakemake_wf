# vim: set ft=python:


checkpoint dRep:
    input:
        checkm=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
    output:
        drep=directory(join(DATA_DIR, binning_analyses, "singlerun/dRep/dereplicated_genomes")),
        out=join(DATA_DIR, binning_analyses, "singlerun/dRep/data_tables/Sdb.csv"),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/Containers:drep"
    params:
        indir=join(DATA_DIR, binning_analyses, "singlerun/all_metawrap_bins/"),
        outdir=join(DATA_DIR, binning_analyses, "singlerun/dRep/"),
    shell:
        """
         rm -rf {params.outdir}
         dRep dereplicate -p {threads} {params.outdir} -g {params.indir}/*.fa \
        -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo {input.checkm} -comp 50 -con 5
        """


rule combine_checkm:
    input:
        checkm=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
        checkm_coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
    output:
        checkm=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm_metrics_combined.csv"),
    params:
        checkm=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm_metrics_combined_tmp.csv"),
    shell:
        """
        awk 'FNR!=1' {input.checkm} {input.checkm_coas}>{params.checkm}
        echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.checkm} > {output}
        rm {params.checkm}
        """


checkpoint dRep_coas:
    input:
        checkm=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm_metrics_combined.csv"),
        out=join(DATA_DIR, binning_analyses, "singlerun/dRep/data_tables/Sdb.csv"),
    output:
        drep=directory(join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/dereplicated_genomes")),
        out=join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/data_tables/Sdb.csv"),
    threads: workflow.cores
    singularity:
        "shub://sskashaf/Containers:drep"
    params:
        indir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/bins_dRep_singlerun_coas/"),
        singlerun_dRep=join(DATA_DIR, binning_analyses, "singlerun/dRep/dereplicated_genomes"),
        all_coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/all_metawrap_bins_coas"),
        outdir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/"),
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.indir}
        rsync {params.singlerun_dRep}/* {params.indir}
        rsync {params.all_coas}/* {params.indir}
        dRep dereplicate -p {threads} {params.outdir} -g {params.indir}/*.fa \
        -pa 0.9 -sa 0.95 -nc 0.30 -cm larger --genomeInfo {input.checkm} -comp 50 -con 5
        """


rule GTDB_TK:
    input:
        done=join(DATA_DIR, binning_analyses, "singlerun_coassembly/dRep/data_tables/Sdb.csv"),
        gtdbrelease=config["db"]["gtdb_db"],
    output:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/GTDB/gtdbtk.bac120.summary.tsv"),
    threads: workflow.cores
    params:
        indir=directory(join(DATA_DIR, binning_analyses, "singlerun_coassembly/dereplicated_genomes/")),
        outdir=join(DATA_DIR, binning_analyses, "singlerun_coassembly/GTDB/"),
        ext="fa"
    singularity:
        "shub://sskashaf/Containers:gtdbtk"
    shell:
        """
        real=$(realpath {input.gtdbrelease})
        rm -rf {params.outdir} 
        export GTDBTK_DATA_PATH=${{real}}
        gtdbtk classify_wf --cpus {threads} --genome_dir {params.indir} --out_dir {params.outdir} -x {params.ext}
        """


rule plot_GTDB:
    input:
        join(DATA_DIR, binning_analyses, "singlerun_coassembly/GTDB/gtdbtk.bac120.summary.tsv"),
    output:
        join(DATA_DIR, "figures/gtdb_bacteria.png"),
    singularity:
        "shub://sskashaf/Containers:r"
    shell:
        "Rscript scripts/plotting/plot_gtdb.R {input}"
