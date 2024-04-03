rule freebayes:
    input:
        alns=lambda w: set(read_mapping.get_collect_bams_input(w)),
        idxs=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w)],
        ref=config["genome"],
    output:
        vcf = "results/calls/calls.vcf",
    log:
        "results/logs/freebayes/calls.log",
    threads: 16
    resources:
        mem_mb=1024,
    wrapper:
        "v3.7.0/bio/freebayes"
