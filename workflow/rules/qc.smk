rule bcftools_stats:
    input:
        lambda w: get_bcftools_stats_input(w.caller),
    output:
        "results/qc/bcftools_stats/{caller}.stats.log",
    log:
        "results/logs/bcftools_stats/{caller}.log",
    params:
        extra=lambda w: f"-s- --depth {config['bcftools_stats']['min_depth']},{config['bcftools_stats']['max_depth']},{config['bcftools_stats']['depth_bin']}",
    wrapper:
        "v5.10.0/bio/bcftools/stats"


rule multiqc:
    input:
        get_base_multiqc_input,
        rules.bcftools_stats.output,
    output:
        "results/qc/multiqc_{caller}/multiqc.html",
        directory("results/qc/multiqc_{caller}/multiqc_data"),
    params:
        extra="--data-dir",  # Optional: extra parameters for multiqc.
        use_input_files_only=True,
    log:
        "results/logs/multiqc_{caller}.log",
    wrapper:
        "v5.10.0/bio/multiqc"
