from snakemake.utils import validate
import pandas as pd
import numpy as np
import os
import re

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


def concat_vcfs_input(w):
    with checkpoints.samtools_faidx.get(**w).output[0].open() as f:
        chroms = [l.split("\t")[0] for l in f.readlines()]
    return expand(
        "results/freebayes/variants/vcfs/{chrom}/variants.{i}.vcf",
        chrom=chroms,
        i=range(1, config["freebayes"]["chunks"] + 1),
    )


def concat_gvcfs_input(w):
    with checkpoints.samtools_faidx.get(**w).output[0].open() as f:
        chroms = [l.split("\t")[0] for l in f.readlines()]
    return expand(
        "results/freebayes/variants/vcfs/{chrom}/variants.{i}.gvcf",
        chrom=chroms,
        i=range(1, config["freebayes"]["chunks"] + 1),
    )
