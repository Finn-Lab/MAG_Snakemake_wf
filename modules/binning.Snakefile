# vim: set ft=python:


def get_sample_reads_1(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    reads1 = df.loc[df["coassembly"] == wildcards.sample, "r1"]
    r1 = reads1.str.split(",", expand=True)
    r1 = r1.values.tolist()
    return r1


def get_sample_reads_2(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    reads2 = df.loc[df["coassembly"] == wildcards.sample, "r2"]
    r2 = reads2.str.split(",", expand=True)
    r2 = r2.values.tolist()
    return r2


def metawrap_cmmd(wildcards):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    reads2 = df.loc[df["coassembly"] == wildcards.sample, "r2"]
    r2 = reads2.str.cat(sep=" ")
    r2 = r2.replace(",", " ")
    reads1 = df.loc[df["coassembly"] == wildcards.sample, "r1"]
    r1 = reads1.str.cat(sep=" ")
    r1 = r1.replace(",", " ")
    cmmd = r1, " ", r2
    return cmmd


checkpoint metawrap_binning:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{sample}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{sample}_2.fastq"),
        fasta=join(DATA_DIR, assembly_dir, "singlerun/{sample}/scaffolds.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/done.txt"),
    singularity:
        "shub://sskashaf/Containers:metawrap"
    params:
        outdir=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap"),
        mincontiglength=2500,
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        rm -rf {params.outdir}
        metawrap binning -t {threads} -m {resources.mem}\
        -a {input.fasta} --maxbin2 --metabat2 --concoct \
        -l {params.mincontiglength} -o {params.outdir} {input.fwd} {input.rev}
        touch {output.outfile}
        """


checkpoint metawrap_binning_coas:
    input:
        unpack(get_sample_reads_1),
        unpack(get_sample_reads_2),
        join(DATA_DIR, assembly_dir, "coassembly/{sample}/scaffolds.fasta"),
    output:
        outfile=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/done.txt"),
    params:
        assembly=join(DATA_DIR, assembly_dir, "coassembly/{sample}/scaffolds.fasta"),
        outdir=join(DATA_DIR, binning_dir, "coassembly/{sample}/metawrap/"),
        reads=lambda wildcards: metawrap_cmmd(wildcards),
        mincontiglength=2500,
    singularity:
        "shub://sskashaf/Containers:metawrap"
    threads: workflow.cores
    resources:
        time=lambda wildcards, attempt: 20 * attempt,
        mem=lambda wildcards, attempt: 80 * attempt,
    shell:
        """
        rm -rf {params.outdir}
        metawrap binning -t {threads} -m {resources.mem} \
        -a {params.assembly} --metabat2 --maxbin2 --concoct \
        -l {params.mincontiglength} -o {params.outdir} {params.reads}
        touch {output}
        """
