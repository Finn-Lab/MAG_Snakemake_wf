# -*- coding: utf-8

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

#'''
#   This is a basic framework for recovery and basic quality control of MAGs
#   To visualize the pipeline: snakemake --dag | dot -Tpng > dag.png
#'''

__maintainer__ = "Sara Kashaf"
__email__ = "sskashaf@ebi.ac.uk"


import os
from os.path import join
import sys
import glob
import pandas as pd
import csv

configfile: "config.yaml"

# Directory structure
DATA_DIR = "data" #config['data']
preprocessing_dir = "00_preprocessing"
assembly_dir = "01_assembly"
binning_dir = "02_binning"
binning_analyses = "03_binning_analyses"
if not os.path.exists("logs"):
    os.makedirs("logs")
os.system("chmod +x scripts/*")
os.system("chmod +x scripts/plotting/*")

# LOAD METADATA
df_run = pd.read_csv("runs.txt")
RUN = df_run["Run"]

df_coas = pd.read_csv("coassembly_runs.txt", sep="\t")
COAS = df_coas["coassembly"]


all_outfiles = [
     # Figure 2
     join(DATA_DIR,preprocessing_dir, "raw_qc/multiqc/raw_multiqc_report.html"),
     join(DATA_DIR,preprocessing_dir, "postprocessing_qc/multiqc/post_multiqc_report.html"),
     # Figure 3a
     join(DATA_DIR,"figures/cmseq_plot.png"),
     join(DATA_DIR,"figures/checkm_contam.png"),
     join(DATA_DIR,"figures/checkm_completeness.png"),
     # Figure 3b
     join(DATA_DIR, "figures/dnadiff.png"),
     # Figure 3c
     join(DATA_DIR,"figures/gtdb_bacteria.png"),
     # Figure 4
     join(DATA_DIR,"figures/perassemb_perref.png")
]

rule all:
   input: all_outfiles


include: "modules/sra_download.Snakefile"
include: "modules/preprocessing.Snakefile"
include: "modules/coas.Snakefile"
include: "modules/assembly.Snakefile"
include: "modules/binning.Snakefile"
include: "modules/refine.Snakefile"
include: "modules/dRep_GTDB.Snakefile"
include: "modules/framework.Snakefile"
include: "modules/refine_coas.Snakefile"
include: "modules/cmseq.Snakefile"
include: "modules/dnadiff.Snakefile"
