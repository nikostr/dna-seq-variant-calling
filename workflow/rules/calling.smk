rule samtools_index:
    input:
        config["genome"],
    output:
        config["genome"] + ".fai",
    log:
        "results/logs/samtools_faidx/index.log",
    wrapper:
        "v3.7.0/bio/samtools/faidx"


rule freebayes:
    input:
        alns=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
        idxs=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f],
        ref=config["genome"],
        ref_idx=rules.samtools_index.output,
    output:
        vcf = "results/freebayes/calls.vcf",
    log:
        "results/logs/freebayes/calls.log",
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
        "v3.7.0/bio/freebayes"
