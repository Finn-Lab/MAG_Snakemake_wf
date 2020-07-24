# vim: set ft=python


checkpoint spades:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq"),
    output:
        join(DATA_DIR, assembly_dir, "singlerun/{run}/scaffolds.fasta"),
        join(DATA_DIR, assembly_dir, "singlerun/{run}/contigs.fasta"),
    threads: workflow.cores
    params:
        outdir=join(DATA_DIR, assembly_dir, "singlerun/{run}/"),
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem=lambda wildcards, attempt: 100 * attempt,
    singularity:
        "shub://sskashaf/MAG_wf_containers:assembly"
    shell:
        """
        rm -rf {params.outdir}
        dmesg -T
        spades.py --meta -1 {input.fwd} -2 {input.rev} \
        -o {params.outdir} -t {threads} -m {resources.mem}
        """


checkpoint spades_coas:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq.gz"),
    output:
        join(DATA_DIR, assembly_dir, "coassembly/{run}/scaffolds.fasta"),
        join(DATA_DIR, assembly_dir, "coassembly/{run}/contigs.fasta"),
    resources:
        time=lambda wildcards, attempt: 10 * attempt,
        mem=lambda wildcards, attempt: 100 * attempt,
    threads: workflow.cores
    params:
        outdir=join(DATA_DIR, assembly_dir, "coassembly/{run}/"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:assembly"
    shell:
        """
        rm -rf {params.outdir}
        dmesg -T
        spades.py --meta -1 {input.fwd} -2 {input.rev} \
        -o {params.outdir} -t {threads} -m {resources.mem}
        """
