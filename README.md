## Getting Started

This protocol assumes that you have two files, which detail the co-assembly and single runs that you think are most suitable for your analyses. The file coassembly_runs.txt specifies the co-assembly samples. This file is in tabular format with three columns, the first column specifies the name of the resulting co-assembly sample. The r1 and r2 columns specify the path of the forward and reverse reads constituting each co-assembly sample, with each read path separated by a comma. The forward and reverse reads must have the extensions _1.fastq and _2.fastq respectively. If the runs are not present locally the sra_download module will attempt to download the runs from the SRA using their SRA accession. The file runs.txt has the SRA accession for the single run samples with a different accession in each line. These files currently include a small subset of gut dataset previously examined by Almeida et al. There are 44 single runs and 3 co-assembly samples based on the metadata of age and geography. If the files are present locally, they can be placed in the path data/raw.



### Setup workflow

Download the GTDB database using:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz 
tar -xvzf gtdbtk_r89_data.tar.gz
```

Download all RefSeq bacterial genomes using:
```
ncbi-genome-download bacteria --formats fasta --section refseq
```
Next generate a Mash sketch of the database with default k-mer and sketch size from the main directory using:

```
mash sketch -o refseq.msh /path/to/RefSeq/*fasta
```
Move the GTDB database and the Mash sketch to the subfolder /data/databases/


### Installing

Install snakemake into an environment using:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake=5.18
```

Then activate the environment before using snakemake: 
```
conda activate snakemake
```

Install Singularity version 3.5.0


### Submitting jobs


Submit jobs with SLURM scheduler:
```
snakemake --use-singularity --restart-times 3 -k -j 50 --cluster-config clusterconfig.yaml --cluster "sbatch -n {cluster.nCPU} --mem {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```

Submit jobs with LSF scheduler:
```
snakemake --use-singularity --restart-times 3 -k --jobs 50 --cluster-config clusterconfig.yaml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```


