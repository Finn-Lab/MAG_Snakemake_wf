# MAG Snakemake Workflow


## Contents

- [Overview](#overview)
- [System Requirements](#system-requirements)
- [Workflow setup](#Workflow-setup)
- [Running pipeline](#running-pipeline)
- [CPU time](#CPU-time)
- [License](./LICENSE)
- [Issues](https://github.com/Finn-Lab/MAG_Snakemake_wf/issues)
- [Citation](#citation)

# Overview

This pipeline can be used for recovery and quality assessment of prokaryotic MAGs from short-read, host-associated metagenomic datasets. The data analyzed is specificed using two files, which detail the co-assembly samples and single runs that are to be analyzed. The file runs.txt has the SRA accession for the single runs with a different accession in each line. The file coassembly_runs.txt specifies the co-assembly samples. This file is in tabular format with three columns, the first column specifies the name of the resulting co-assembly sample. The r1 and r2 columns specify the path of the forward and reverse reads constituting each co-assembly sample, with each read path separated by a comma. In this pipeline, we used a small subset of gut dataset previously analyzed by Almeida et al. There are 40 single runs and 2 co-assembly samples. The co-assembly samples are named based on the metadata that was chosen to perform co-assembly. If the runs are not publicly available, they should be placed in the subdirectory data/raw with respect to the Snakefile. The forward and reverse reads must have the extensions _1.fastq and _2.fastq respectively. Note that the version of kneaddata used in this pipeline still requires the read suffixes /1 and /2, so the headers must be formatted accordingly. If the runs are publicly available, specify them in the runs.txt file and the sra_download module will download the runs from the SRA using their SRA accession to the directory data/raw. 

# System Requirements

## Hardware Requirements

HPC with at least 500 gigabytes of memory

The CPU times below are generated using the cluster config file included in the repo

## Software Requirements

MAG Snakemake pipeline (https://github.com/Finn-Lab/MAG_Snakemake_wf)

Singularity 3.5.0 (https://github.com/hpcng/singularity)

Snakemake (version 5.18) (https://github.com/snakemake/snakemake) 

Running the MAG Snakemake pipeline will automatically download the sequencing data from the SRA. It will also download the relevant singularity containers so the relevant software needed for our pipeline can be used. Alternatively, the tools can be manually downloaded from:

ncbi-genome-download (version 0.3.0) (https://github.com/kblin/ncbi-genome-download)

mash (version 2.2.1) (https://github.com/marbl/Mash)

parallel-fastq-dump (version 0.6.6) & fastq-dump (version 2.8.0) (https://github.com/rvalieris/parallel-fastq-dump)

fastqc (version 0.11.7) (https://github.com/s-andrews/FastQC)

multiqc (version 1.3) (https://github.com/ewels/MultiQC)

kneaddata (version 0.7.4) with Trimmomatic (version 0.39) & Bowtie (version 2.4.2) (https://github.com/biobakery/kneaddata)

metaSPAdes (version 3.14.0) (https://github.com/ablab/spades)

metaWRAP (version 1.2.2) (https://github.com/bxlab/metaWRAP)

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

GTDB database (release 95) (https://data.ace.uq.edu.au/public/gtdb/data/releases/) 


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
mv /path/to/gtdbtk_r95_data.tar.gz data/databases
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




# Citation
For a walk-through of this pipeline, please visit: 
Saheb Kashaf, S., Almeida, A., Segre, J.A. et al. Recovering prokaryotic genomes from host-associated, short-read shotgun metagenomic sequencing data. Nat Protoc (2021). https://doi.org/10.1038/s41596-021-00508-2

To generate results from the paper associated with this pipeline, please select a figure to reproduce in the Snakefile.



