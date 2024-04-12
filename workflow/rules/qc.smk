rule bcftools_stats:
    input:
        rules.concat_vcfs.output,
    output:
        "results/qc/bcftools_stats/calls.stats.log",
    log:
        "results/logs/bcftools_stats/bcftools_stats.log",
    params:
        extra="-s-"
    wrapper:
        "master/bio/bcftools/stats"

