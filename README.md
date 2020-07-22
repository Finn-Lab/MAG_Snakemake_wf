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
