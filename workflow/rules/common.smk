from snakemake.utils import validate
import pandas as pd
import numpy as np
import os
import re
import csv

##### load config and sample sheets #####


configfile: "config/config.yaml"


# validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], dtype=str).set_index(
    ["sample_id", "datatype", "unit"], drop=False
)
samples.index = samples.index.set_levels(
    [i.astype(str) for i in samples.index.levels]
)  # enforce str in index
# validate(samples, schema="../schemas/samples.schema.yaml")


def get_concat_bcftools_vcfs_input(w):
    with checkpoints.generate_freebayes_regions.get(**w).output[0].open() as f:
        regions = [
            r
            for r in csv.DictReader(
                f,
                fieldnames=["chr", "start", "end"],
                dialect="excel-tab",
            )
        ]
    # Handle bcftools region format by adjusting start values
    return [
        f"results/bcftools/calls/vcfs/{r['chr']}/variants.{int(r['start'])+1}-{r['end']}.vcf"
        for r in regions
    ]


def concat_vcfs_input(w):
    with checkpoints.generate_freebayes_regions.get(**w).output[0].open() as f:
        regions = [
            r
            for r in csv.DictReader(
                f,
                fieldnames=["chr", "start", "end"],
                dialect="excel-tab",
            )
        ]
    return [
        f"results/freebayes/variants/vcfs/{r['chr']}/variants.{r['start']}-{r['end']}.vcf"
        for r in regions
    ]


def concat_gvcfs_input(w):
    with checkpoints.generate_freebayes_regions.get(**w).output[0].open() as f:
        regions = [
            r
            for r in csv.DictReader(
                f,
                fieldnames=["chr", "start", "end"],
                dialect="excel-tab",
            )
        ]
    return [
        f"results/freebayes/variants/vcfs/{r['chr']}/variants.{r['start']}-{r['end']}.gvcf"
        for r in regions
    ]


def get_base_multiqc_input(w):
    return (
        read_mapping.get_multiqc_input_from_samtools_stats(w)
        + read_mapping.get_multiqc_input_from_fastqc(w)
        + read_mapping.get_multiqc_input_from_deeptools_plotcoverage(w)
        + expand(
            "results/qc/fastp/{s.sample_id}-{s.unit}_fastp.json",
            s=samples.query('datatype=="illumina"').itertuples(),
        )
    )


def get_bcftools_stats_input(caller):
    input = dict(
        bcftools_vcf=rules.fill_bcftools_vcf_tags.output,
        freebayes_vcf=rules.compress_vcf.output,
        freebayes_gvcf=rules.compress_gvcf.output,
        delly=rules.delly_bcftools_merge.output,
        delly_cnv=rules.bcftools_view_cnv.output,
    )
    return input[caller]
