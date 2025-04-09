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


# use rule some_task from other_workflow as other_some_task with:
#    output:
#        "results/some-result.txt"
use rule multiqc from read_mapping as read_mapping_multiqc with:
    input:
        get_base_multiqc_input,
        rules.bcftools_stats.output,
