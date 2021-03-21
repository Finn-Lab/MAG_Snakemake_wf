# MAG Snakemake Workflow


## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Workflow setup](#Workflow-setup)
- [Running pipeline](#running-pipeline)
- [CPU time](#CPU-time)
- [Results](#results)
- [License](./LICENSE)
- [Issues](https://github.com/Finn-Lab/MAG_Snakemake_wf/issues)
- [Citation](#citation)

# Overview

This pipeline at its current form generates prokaryotic MAGs for a subset of our previous gut analyses. To run this pipeline on your own dataset, you need two files, which detail the co-assembly and single runs that you think are most suitable for your analyses. The file runs.txt has the SRA accession for the single run samples with a different accession in each line. The file coassembly_runs.txt specifies the co-assembly samples. This file is in tabular format with three columns, the first column specifies the name of the resulting co-assembly sample. The r1 and r2 columns specify the path of the forward and reverse reads constituting each co-assembly sample, with each read path separated by a comma. The forward and reverse reads must have the extensions _1.fastq and _2.fastq respectively. These files currently include a small subset of gut dataset previously examined by Almeida et al. There are 44 single runs and 3 co-assembly samples based on the metadata of age and geography. If the files are present locally, they should be placed in the subdirectory data/raw with respect to the Snakefile.  If the runs are not present locally the sra_download module will attempt to download the runs from the SRA using their SRA accession to this directory. 


# System Requirements

## Hardware Requirements

HPC with at least 350 gigabytes of memory

The CPU times below are generated using the cluster config file included in the repo

## Software Requirements

MAG Snakemake pipeline (https://github.com/Finn-Lab/MAG_Snakemake_wf)

Singularity 3.5.0 (https://github.com/hpcng/singularity)

Snakemake (version 5.18) (https://github.com/snakemake/snakemake) 

Running the MAG Snakemake pipeline will automatically download the sequencing data from the SRA. It will also download the relevant singularity containers so the relevant software needed for our pipeline can be used. Alternatively, the tools can be manually downloaded from:

ncbi-genome-download (version 0.3.0) (https://github.com/kblin/ncbi-genome-download)

mash (version 2.2.1) (https://github.com/marbl/Mash)

parallel-fastq-dump (version 0.6.6) & fastq-dump (version 2.8.0) (https://github.com/rvalieris/parallel-fastq-dump)

fastqc=0.11.7 (https://github.com/s-andrews/FastQC)

multiqc=1.3 (https://github.com/ewels/MultiQC)

kneaddata=0.7.4 with Trimmomatic=0.39 & Bowtie=2.4.2 (https://github.com/biobakery/kneaddata)

metaSPAdes (version 3.14.0) (https://github.com/ablab/spades)

metaWRAP=1.2.2 (version ) (https://github.com/bxlab/metaWRAP)

CheckM (version 1.0.12) (https://github.com/Ecogenomics/CheckM)

Bowtie (version 2.4.1) (https://github.com/BenLangmead/bowtie2)

Prokka (version 1.14.5) (https://github.com/tseemann/prokka)

CMSeq (version 1.0) (https://bitbucket.org/CibioCM/cmseq)

mummer (version 3.23) (https://github.com/mummer4/mummer)

dRep (version 2.3.2) (https://github.com/MrOlm/drep)

GTDB_Tk (version 1.2.0) (https://github.com/Ecogenomics/GTDBTk)

bwa (version 0.7.17) (https://github.com/lh3/bwa)

samtools (version 1.9) (https://github.com/samtools/samtools) 

## Other

RefSeq complete bacterial genomes (downloaded May 2020) (https://www.ncbi.nlm.nih.gov/refseq/)

GTDB database (release 89) (https://data.ace.uq.edu.au/public/gtdb/data/releases/) 


# Workflow setup


Download the GTDB database using:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release95/95.0/auxillary_files/gtdbtk_r95_data.tar.gz
```

Download all RefSeq bacterial genomes using:
```
ncbi-genome-download bacteria --formats fasta --section refseq --assembly-levels complete
```
Next generate a Mash sketch of the database with default k-mer and sketch size from the main directory using:

```
mash sketch -o refseq.msh /path/to/RefSeq/*fasta
```

Download the code for the pipeline from (https://github.com/Finn-Lab/MAG_Snakemake_wf) in a location that has at least 1.5 TB of disk space. Change directory to this folder. Move the GTDB database and the RefSeq Mash sketch to the subfolder /data/databases using:

```
cd /path/to/MAG_Snakemake_wf/
mkdir -p data/databases
mv /path/to/refseq.msh data/databases
mv /path/to/gtdbtk_r89_data.tar.gz data/databases
tar -xvzf data/databases/gtdbtk_r95_data.tar.gz
```

Install snakemake into an environment using:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake=5.18
```

Then activate the environment before using snakemake: 
```
conda activate snakemake
```


# Running pipeline 

### Submitting jobs

To run pipeline on the small gut dataset specified in runs.txt and coassembly_runs.txt, submit jobs with SLURM scheduler:
```
snakemake --use-singularity --restart-times 3 -k -j 50 --cluster-config clusterconfig.yaml --cluster "sbatch -n {cluster.nCPU} --mem {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```

Submit jobs with LSF scheduler:
```
snakemake --use-singularity --restart-times 3 -k --jobs 50 --cluster-config clusterconfig.yaml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -e {cluster.error} -o {cluster.output}"
```

# CPU-time

The CPU time for the demo dataset and the cluster configuration file provided is as follows:

Data Download: 16 hours

Preprocessing: 48 hours 

Assembly/Co-assembly:  2600 hours

Binning: 100 hours 

Quality Assessment & Bin Refinement; Estimate completeness and contamination of MAGs: 25 hours

Quality Assessment & Bin Refinement; Estimate strain heterogeneity of MAGs: 700 hours

Quality Assessment & Bin Refinement; Compare MAGs to RefSeq genomes: 1 hour

Quality Assessment & Bin Refinement; Bin refinement: 260 hours

Dereplicate MAGs: 4 hours

Taxonomic Classification: 3 hours

Evaluate Bottlenecks: 120 hours


# Results

To generate results from the associated paper, please select a figure to reproduce in the associated Snakefile.

# Citation




