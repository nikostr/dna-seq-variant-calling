rule samtools_faidx:
    input:
        config["genome"],
    output:
        f"{config["genome"]}.fai",
    log:
        "results/logs/samtools_faidx/index.log",
    wrapper:
        "v3.7.0/bio/samtools/faidx"


rule generate_bam_list:
    input:
        bams=lambda w: set(
            [f for f in read_mapping.get_collect_bams_input(w) if "illumina" in f]
        ),
    output:
        temp("results/freebayes/bams.txt"),
    localrule: True
    log:
        "results/logs/freebayes/generate_bam_list.log",
    run:
        with open(output[0], "w") as f:
            f.writelines([f"{l}\n" for l in input.bams])


checkpoint generate_freebayes_regions:
    input:
        ref_idx=rules.samtools_faidx.output,
        index=lambda w: [
            f"{f}.bai"
            for f in read_mapping.get_collect_bams_input(w)
            if "illumina" in f
        ],
        samples=rules.generate_bam_list.output,
    output:
        "results/freebayes/regions.bed",
    log:
        "results/logs/freebayes/generate_freebayes_regions.log",
    params:
        target_data_size=config["freebayes"]["target_data_size"],
    conda:
        "../envs/freebayes-env.yaml"
    localrule: True
    shell:
        "split_ref_by_bai_datasize.py "
        "--reference-fai {input.ref_idx} "
        "--target-data-size {params.target_data_size} "
        "--bam-list {input.samples} "
        "> {output} "
        "2> {log} "


rule variant_calling_bcftools_mpileup:
    input:
        alignments=lambda w: set(
            [f for f in read_mapping.get_collect_bams_input(w) if "illumina" in f]
        ),
        ref=config["genome"],
        index=rules.samtools_faidx.output[0],
        all_beds=rules.generate_freebayes_regions.output,
    output:
        pileup=pipe("results/bcftools/mpileup/vcfs/{chrom}/variants.{bp}.vcf"),
    params:
        uncompressed_bcf=True,
        extra=lambda w: f"--region {w.chrom}:{w.bp}",
    log:
        "results/logs/bcftools/mpileup/{chrom}.{bp}.log",
    wrapper:
        "v6.1.0/bio/bcftools/mpileup"


rule variant_calling_bcftools_call:
    input:
        pileup=rules.variant_calling_bcftools_mpileup.output.pileup,
    output:
        calls=temp("results/bcftools/calls/vcfs/{chrom}/variants.{bp}.vcf"),
    params:
        caller="-m",
    log:
        "results/logs/bcftools/call/{chrom}.{bp}.log",
    wrapper:
        "v6.1.0/bio/bcftools/call"


rule concat_bcftools_vcfs:
    input:
        vcfs=get_concat_bcftools_vcfs_input,
    output:
        temp("results/bcftools/calls.unnormalized.vcf"),
        temp("results/bcftools/calls.unnormalized.vcf.idx"),
    params:
        extra="",
    log:
        "results/logs/bcftools/picard-merge/out.log",
    wrapper:
        "v6.1.0/bio/picard/mergevcfs"


rule norm_bcftools_vcf:
    input:
        rules.concat_bcftools_vcfs.output[0],
        ref=config["genome"],
    output:
        temp("results/bcftools/calls.normalized.vcf.gz"),
    log:
        "results/logs/bcftools/norm/out.log",
    params:
        extra="--multiallelics +any",
    wrapper:
        "v6.1.0/bio/bcftools/norm"


rule fill_bcftools_vcf_tags:
    input:
        rules.norm_bcftools_vcf.output[0],
    output:
        "results/bcftools/calls.vcf.gz",
    log:
        "results/logs/bcftools/fill-tags/out.log",
    conda:
        "../envs/bcftools.yaml"
    shell:
        "bcftools +fill-tags "
        "{input} "
        "--output-type z "
        "--write-index "
        "--threads {threads} "
        "--output {output} "
        "-- "
        "--tags AN,AC "
        "2> {log} "


rule variant_calling_freebayes:
    input:
        bams=lambda w: set(
            [f for f in read_mapping.get_collect_bams_input(w) if "illumina" in f]
        ),
        index=lambda w: [
            f"{f}.bai"
            for f in read_mapping.get_collect_bams_input(w)
            if "illumina" in f
        ],
        ref=config["genome"],
        samples=rules.generate_bam_list.output,
        all_beds=rules.generate_freebayes_regions.output,
    output:
        temp("results/freebayes/variants/vcfs/{chrom}/variants.{bp}.vcf"),
    params:
        region=lambda w: f"{w.chrom}:{w.bp}",
    log:
        "results/logs/freebayes/variant_calling_freebayes/{chrom}.{bp}.log",
    conda:
        "../envs/freebayes-env.yaml"
    threads: 1
    shell:
        "freebayes "
        "--fasta-reference {input.ref} "
        "--region {params.region} "
        "--bam-list {input.samples} "
        "> {output} "
        "2> {log}"


rule variant_calling_freebayes_gvcf:
    input:
        bams=lambda w: set(
            [f for f in read_mapping.get_collect_bams_input(w) if "illumina" in f]
        ),
        index=lambda w: [
            f"{f}.bai"
            for f in read_mapping.get_collect_bams_input(w)
            if "illumina" in f
        ],
        ref=config["genome"],
        samples=rules.generate_bam_list.output,
        all_beds=rules.generate_freebayes_regions.output,
    output:
        temp("results/freebayes/variants/vcfs/{chrom}/variants.{bp}.gvcf"),
    params:
        region=lambda w: f"{w.chrom}:{w.bp}",
    log:
        "results/logs/freebayes/variant_calling_freebayes_gvcf/{chrom}.{bp}.log",
    conda:
        "../envs/freebayes-env.yaml"
    threads: 1
    shell:
        "freebayes "
        "--fasta-reference {input.ref} "
        "--region {params.region} "
        "--bam-list {input.samples} "
        "--gvcf "
        "> {output} "
        "2> {log}"


rule concat_vcfs:
    input:
        calls=concat_vcfs_input,
    output:
        vcf=temp("results/freebayes/calls.vcf"),
    log:
        "results/logs/concat_vcfs/log.log",
    conda:
        "../envs/freebayes-env.yaml"
    threads: 4
    shell:
        "bcftools concat {input.calls} 2> {log} | vcfuniq > {output} 2>> {log}"


rule compress_vcf:
    input:
        rules.concat_vcfs.output,
    output:
        "results/freebayes/calls.vcf.gz",
    log:
        "results/logs/compress_vcf/log.log",
    params:
        extra="--write-index",
    wrapper:
        "v5.10.0/bio/bcftools/view"


rule concat_gvcfs:
    input:
        calls=concat_gvcfs_input,
    output:
        vcf=temp("results/freebayes/calls.g.vcf"),
    log:
        "results/logs/concat_vcfs/log.log",
    conda:
        "../envs/freebayes-env.yaml"
    threads: 4
    shell:
        "bcftools concat {input.calls} 2> {log} | vcfuniq > {output} 2>> {log}"


rule compress_gvcf:
    input:
        rules.concat_gvcfs.output,
    output:
        "results/freebayes/calls.g.vcf.gz",
    log:
        "results/logs/compress_vcf/log.log",
    params:
        extra="--write-index",
    wrapper:
        "v5.10.0/bio/bcftools/view"
