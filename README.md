## Getting Started

These instructions will help you setup this workflow for your own analyses

### Setup workflow

Download the GTDB database using:
```
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release89/89.0/gtdbtk_r89_data.tar.gz 
tar xvzf gtdbtk_r89_data.tar.gz
```

Download the RefSeq database using:
```
ncbi-genome-download bacteria --formats fasta --section refseq
```
Next generate a Mash sketch of the database with default k-mer and sketch size using:

```
mash sketch -o refseq.msh /path/to/RefSeq/*fasta
```

### Installing

Install snakemake into an environment using:

```
conda create -c conda-forge -c bioconda -n snakemake snakemake=5.18
```

Then activate the environment before using snakemake: 
```
conda activate snakemake
```

### Submitting jobs

Submit jobs with SLURM scheduler:
```
snakemake --use-singularity --restart-times 3 -k -j 50 --cluster-config clusterconfig.yaml --cluster "sbatch -n {cluster.nCPU} --mem {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```

Submit jobs with LSF scheduler:
```
snakemake --use-singularity --restart-times 3 -k --jobs 50 --cluster-config clusterconfig.yaml --cluster "bsub -n {cluster.nCPU} -M {cluster.mem} -e {cluster.error} -o {cluster.output} -t {cluster.time}"
```


