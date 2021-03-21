# This file is part of MAG Snakemake workflow.
#
# MAG Snakemake workflow is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# MAG Snakemake workflow is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with MAG Snakemake workflow.  If not, see <https://www.gnu.org/licenses/>.

# vim::w set ft=python:


rule raw_fastqc_fwd:
    input:
        fwd=join(DATA_DIR, "raw/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, "raw/{run}_2.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_1_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.fwd} --outdir {params.outdir}
        """


rule raw_fastqc_rev:
    input:
        fwd=join(DATA_DIR, "raw/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, "raw/{run}_2.fastq.gz"),
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_2_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.rev} --outdir {params.outdir}
        """



rule raw_multiqc:
    input:
        expand(join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/{run}_{read}_fastqc.html"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/raw_multiqc_report.html"),
    params:
        indir=join(DATA_DIR, preprocessing_dir, "raw_qc/fastqc/"),
        outdir=join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/"),
        outfile=join(DATA_DIR, preprocessing_dir, "raw_qc/multiqc/multiqc_report.html"),
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.3--py35_2"
    shell:
        """
        rm -f {params.outfile}
        multiqc {params.indir} -o {params.outdir}
        mv {params.outfile} {output}
        """


rule kneaddata_download_database:
    output:
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.1.bt2"),
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.2.bt2"),
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.3.bt2"),
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.4.bt2"),
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.rev.1.bt2"),
        join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.rev.2.bt2"),
    params:
        outdir="data/databases/human_genome_index/",
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:metagenome_preprocessing"
    shell:
        """
        kneaddata_database --download human_genome bowtie2 {params.outdir}
        """

rule kneaddata_bowtie:
    input:
        fwd=join(DATA_DIR, "raw/{run}_1.fastq.gz"),
        rev=join(DATA_DIR, "raw/{run}_2.fastq.gz"),
        indx1=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.1.bt2"),
        indx2=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.2.bt2"),
        indx3=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.3.bt2"),
        indx4=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.4.bt2"),
        indx5=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.rev.1.bt2"),
        indx6=join(DATA_DIR, "databases/human_genome_index/hg37dec_v0.1.rev.2.bt2"),
    output:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq"),
    params:
        tmp_fwd=join(DATA_DIR, "raw/{run}_tmp_1.fastq"),
        tmp_rev=join(DATA_DIR, "raw/{run}_tmp_2.fastq"),
        tmp_fwd2=join(DATA_DIR, "raw/{run}_tmp2_1.fastq.gz"),
        tmp_rev2=join(DATA_DIR, "raw/{run}_tmp2_2.fastq.gz"),
        fwd=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/{run}_tmp2_1_kneaddata_paired_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/{run}_tmp2_1_kneaddata_paired_2.fastq"),
        outdir=directory(join(DATA_DIR, preprocessing_dir, "kneaddata_bowtie/")),
        indx=join(DATA_DIR, "databases/human_genome_index/"),
    singularity:
        "shub://sskashaf/MAG_wf_containers_2021:metagenome_preprocessing"
    threads: workflow.cores
    shell:
        """
        #this version of kneaddata requires read identifiers (/1, /2). Given we did not download the file from the sra using the option --readids, we need to remove spaces in the headers
        # and then add /1 and /2 using reformat.sh. Kneaddata by default does not sort the reads so Next we sort the reads. 
        # Alternatively, one could sort using kneaddata with the option --reorder
        rm -f {params.tmp_fwd} {params.tmp_rev} {params.tmp_fwd2} {params.tmp_rev2} {params.fwd} {params.rev} {output.fwd} {output.rev}
        mkdir -p {params.outdir}
        seqtk seq -C {input.fwd} > {params.tmp_fwd}
        seqtk seq -C {input.rev} > {params.tmp_rev}
        reformat.sh in={params.tmp_fwd} in2={params.tmp_rev} out1={params.tmp_fwd2} out2={params.tmp_rev2} addslash spaceslash=f
        kneaddata --remove-intermediate-output --threads {threads} \
        --input {params.tmp_fwd2} --input {params.tmp_rev2}\
        --output {params.outdir} \
        --reference-db {params.indx} \
        --trimmomatic-options "ILLUMINACLIP:/data/adapters/TruSeq3-PE.fa:2:30:10: SLIDINGWINDOW:4:20 MINLEN:50" --trimmomatic /data/\
        --bowtie2-options "--very-sensitive --dovetail" 
        repair.sh in={params.fwd} in2={params.rev} out={output.fwd} out2={output.rev} repair
        rm {params.tmp_fwd} {params.tmp_rev} {params.tmp_fwd2} {params.tmp_rev2}
        """


rule postpreprocessing_fastqc_fwd:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq"),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_1_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.fwd} --outdir {params.outdir}
        """

rule postpreprocessing_fastqc_rev:
    input:
        fwd=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_1.fastq"),
        rev=join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_2.fastq"),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_2_fastqc.html"),
    params:
        outdir=directory(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/")),
    threads: workflow.cores
    singularity:
        "docker://quay.io/biocontainers/fastqc:0.11.7--4"
    shell:
        """
        fastqc {input.rev} --outdir {params.outdir}
        """



rule postpreprocessing_multiqc:
    input:
        expand(join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/{run}_{read}_fastqc.html"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/post_multiqc_report.html"),
    params:
        indir=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/fastqc/"),
        outdir=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/"),
        outfile=join(DATA_DIR, preprocessing_dir, "postprocessing_qc/multiqc/multiqc_report.html"),
    singularity:
        "docker://quay.io/biocontainers/multiqc:1.3--py35_2"
    shell:
        """
        multiqc --force {params.indir} -o {params.outdir}
        mv {params.outfile} {output}
        """


def linecount(fastq):
    out = subprocess.Popen("zcat -f "+fastq+" | wc -l", stdout=subprocess.PIPE, shell=True).communicate()[0]
    count=int(out.strip())/4
    return count

rule readcount_fq:
    input:
        raw=expand(join(DATA_DIR, preprocessing_dir, "processed/singlerun/{run}_{read}.fastq"), run=RUN, read=["1", "2"]),
    output:
        join(DATA_DIR, preprocessing_dir, "readcounts.tsv"),
    run:
        outfile = str(output)
        if os.path.exists(outfile):
            os.remove(outfile)
        with open(outfile, "w") as outf:
            outf.writelines("Run\tReadcount\n")
            for run in input:
                readcount = int(linecount(run))
                line = run + "\t" + str(readcount)
                outf.writelines(line + "\n")

