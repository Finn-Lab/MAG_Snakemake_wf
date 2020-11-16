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
