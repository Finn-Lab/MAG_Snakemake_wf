# vim::w set ft=python:


rule sra_download:
    output:
        join(DATA_DIR, "raw/{run}_1.fastq.gz"),
        join(DATA_DIR, "raw/{run}_2.fastq.gz"),
    resources:
        time=5,
        mem=40,
    threads: workflow.cores
    singularity: "docker://quay.io/biocontainers/parallel-fastq-dump:0.6.3--py36_1"
#"shub://sskashaf/MAG_wf_containers:fastqdump"
    params:
        outdir=join(DATA_DIR, "raw/"),
        run="{run}",
    shell:
        """
        echo {threads}
        prefetch {params.run} && vdb-validate {params.run} && parallel-fastq-dump --threads {threads} \
        --outdir {params.outdir} --skip-technical --split-3 --sra-id {params.run} --gzip
        """
