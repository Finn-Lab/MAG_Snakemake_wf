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

# vim: set ft=python:


rule bin_refinement:
    input:
        join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/done.txt"),
    output:
        join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/done.txt"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:metawrap"
    params:
        metabat=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/metabat2_bins"),
        maxbin=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/maxbin2_bins"),
        concoct=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap/concoct_bins"),
        outdir=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement"),
    threads: workflow.cores
    shell:
        """
        checkm data setRoot ~/checkm_database/ 
        rm -rf {params.outdir}
        metawrap bin_refinement -o {params.outdir} -t {threads} -A {params.metabat} -B {params.maxbin} -C {params.concoct} -c 50 -x 10
        touch {output}
        """


checkpoint refine_bins_init:
    input:
        rules.bin_refinement.output,
    output:
        directory(join(DATA_DIR, binning_dir, "singlerun/{sample}/refined_bins_50_10")),
    params:
        refined_folder=join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement/metawrap_50_10_bins"),
        folder=directory(join(DATA_DIR, binning_dir, "singlerun/{sample}/refined_bins_50_10")),
    shell:
        """
        rm -rf {params.folder}
        mkdir -p {params.folder}
        cp {params.refined_folder}/*.fa {params.folder}
        """


def aggregate_bins(wildcards):
    refine_out = checkpoints.refine_bins_init.get(**wildcards).output[0]
    return glob_wildcards(os.path.join(refine_out, "{bin}.fa")).bin



rule copy_bins:
    input:
        join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement/metawrap_50_10_bins/{i}.fa"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/all_metawrap_bins/{sample}_metawrap_refined_{i}.fa"),
        join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    params:
        dir1=join(DATA_DIR, binning_analyses, "singlerun/all_metawrap_bins/"),
        dir2=join(DATA_DIR, "all_bins/"),
        out1=join(DATA_DIR, binning_analyses, "singlerun/all_metawrap_bins/{sample}_metawrap_refined_{i}.fa"),
        out2=join(DATA_DIR, "all_bins/{sample}_metawrap_refined_{i}.fa"),
    shell:
        """
        mkdir -p {params.dir1}
        mkdir -p {params.dir2}
        cp {input} {params.out1}
        cp {input} {params.out2}
        """


rule copy_refined:
    input:
        lambda wildcards: expand(rules.copy_bins.output, i=aggregate_bins(wildcards), sample=wildcards.sample),
    output:
        join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement/copied.txt"),
    shell:
        """
        touch {output}
        """


rule checkm:
    input:
        expand(join(DATA_DIR, binning_dir, "singlerun/{sample}/metawrap_bin_refinement/copied.txt"), sample=RUN), 
    output:
        join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.tsv"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:metawrap"
    params:
        ext="fa",
        indir=join(DATA_DIR, binning_analyses, "singlerun/all_metawrap_bins"),
        outdir=join(DATA_DIR, binning_analyses, "singlerun/checkm/"),
    shell:
        """
        checkm data setRoot ~/checkm_database/
        rm -rf {params.outdir}
        checkm lineage_wf -t {threads} -x {params.ext} --tab_table -f {output} {params.indir} {params.outdir}
        """


rule parse_checkm:
    input:
        join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.tsv"),
    output:
        join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
    params:
        checkm=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics_tmp.csv"),
        checkm2=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics_tmp2.csv"),
    shell:
        """
        sed -i '1d' {input}
        cut -f1,12,13,14 {input} | tr '\t' ','>{params.checkm}
        sed 's/,/.fa,/' {params.checkm}>{params.checkm2}
        echo -e "genome,completeness,contamination,strain_heterogeneity" | cat - {params.checkm2} > {output}
        rm {params.checkm}
        rm {params.checkm2}
        """


rule plot_checkm:
    input:
        sr=join(DATA_DIR, binning_analyses, "singlerun/checkm/checkm_metrics.csv"),
        coas=join(DATA_DIR, binning_analyses, "singlerun_coassembly/checkm/checkm_metrics.csv"),
    output:
        barplot=join(DATA_DIR, "figures/checkm_completeness.png"),
        barplot2=join(DATA_DIR, "figures/checkm_contam.png"),
    singularity:
        "shub://sskashaf/MAG_wf_containers:r"
    shell:
        """
        Rscript scripts/plotting/plot_checkm_mags.R {input.sr} {input.coas}
        """
