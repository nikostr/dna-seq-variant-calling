rule bcftools_stats:
    input:
        rules.compress_vcf.output,
    output:
        "results/qc/bcftools_stats/freebayes_vcf.stats.log",
    log:
        "results/logs/bcftools_stats/freebayes_vcf.log",
    params:
        extra=lambda w: f"-s- --depth {config['bcftools_stats']['min_depth']},{config['bcftools_stats']['max_depth']},{config['bcftools_stats']['depth_bin']}",
    wrapper:
        "v5.10.0/bio/bcftools/stats"
