# vim: set ft=python:


def get_sample_reads(sample):
    sample_reads = []
    sample_file = "coassembly_runs.txt"
    df = pd.read_csv(sample_file, sep="\t")
    r1 = list(df[df["coassembly"] == sample]["r1"])[0].split(",")
    r2 = list(df[df["coassembly"] == sample]["r2"])[0].split(",")
    dict = {"r1": r1, "r2": r2}
    return dict


rule coas:
    input:
        unpack(lambda wildcards: get_sample_reads(wildcards.sample)),
    output:
        r1=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{sample}_1.fastq"),
        r2=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{sample}_2.fastq"),
    shell:
        """
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}
        """


rule sort:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq"),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq"),
    singularity:
        "modules/envs/singularity/assembly"
    shell:
        """
        repair.sh in={input.fwd} in2={input.rev} out={output.fwd} out2={output.rev} repair
        rm -f {input.fwd}
        rm -f  {input.rev}
        """


rule zip_file:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq"),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_2.fastq.gz"),
    params:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_1.fastq.gz"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_sorted_2.fastq.gz"),
    shell:
        """
        rm -f {params.fwd}
        rm -f {params.rev}
        rm -f {output.fwd}
        rm -f {output.rev}
        gzip {input.fwd}
        gzip {input.rev}
        mv {params.fwd} {output.fwd}
        mv {params.rev} {output.rev} 
        """




rule readcounts_coas:
    input:
        raw=expand(join(DATA_DIR, preprocessing_dir, "processed/coassembly/{run}_{read}.fastq"), run=COAS, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "readcounts_coas.tsv"),
    run:
        outfile = str(output)
        if os.path.exists(outfile):
            os.remove(outfile)
        with open(outfile, "w") as outf:
            outf.writelines("Run\tReadcount\n")
            for run in input:
                readcount = int(linecount(run))
                line = "\t".join([run, str(readcount)])
                outf.writelines(line, "\n")
