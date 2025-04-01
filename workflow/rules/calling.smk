checkpoint samtools_faidx:
    input:
        config["genome"],
    output:
        config["genome"] + ".fai",
    log:
        "results/logs/samtools_faidx/index.log",
    wrapper:
        "v3.7.0/bio/samtools/faidx"


rule generate_freebayes_regions:
    input:
        ref_idx=rules.samtools_faidx.output,
    output:
        directory("results/freebayes/regions")
    log:
        "results/logs/freebayes/generate_freebayes_regions.log"
    params:
        chunks = config["freebayes"]["chunks"],
    localrule: True
    script:
        "../scripts/fasta_generate_regions.py"


rule generate_bam_list:
    input:
        bams=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
    output:
        temp("results/freebayes/bams.txt")
    localrule: True
    run:
        with open(output[0], 'w') as f:
            f.writelines([l + '\n' for l in input.bams])


rule variant_calling_freebayes:
    input:
        bams=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
        index=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f],
        ref=config["genome"],
        samples=rules.generate_bam_list.output,
        all_beds=rules.generate_freebayes_regions.output,
    output:
        temp("results/freebayes/variants/vcfs/{chrom}/variants.{i}.vcf")
    params:
        regions="results/freebayes/regions/genome.{chrom}.region.{i}.bed"
    log:
        "results/logs/freebayes/variant_calling_freebayes/{chrom}.{i}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads: 1
    shell: "freebayes -f {input.ref} -t {params.regions} -L {input.samples} > {output} 2> {log}"


rule variant_calling_freebayes_gvcf:
    input:
        bams=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
        index=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f],
        ref=config["genome"],
        samples=rules.generate_bam_list.output,
        all_beds=rules.generate_freebayes_regions.output,
    output:
        temp("results/freebayes/variants/vcfs/{chrom}/variants.{i}.gvcf")
    params:
        regions="results/freebayes/regions/genome.{chrom}.region.{i}.bed"
    log:
        "results/logs/freebayes/variant_calling_freebayes_gvcf/{chrom}.{i}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads: 1
    shell: "freebayes -f {input.ref} -t {params.regions} -L {input.samples} --gvcf > {output} 2> {log}"



def concat_vcfs_input(w):
    with checkpoints.samtools_faidx.get(**w).output[0].open() as f:
        chroms = [l.split('\t')[0] for l in f.readlines()]
    return expand(
            "results/freebayes/variants/vcfs/{chrom}/variants.{i}.vcf",
            chrom=chroms,
            i=range(1,config["freebayes"]["chunks"]+1)
            )


rule concat_vcfs:
    input:
        calls = concat_vcfs_input,
    output:
        vcf = temp("results/freebayes/calls.vcf"),
    log:
        "results/logs/concat_vcfs/log.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads: 4
    shell:  
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"


rule compress_vcf:
    input:
        rules.concat_vcfs.output,
    output:
        "results/freebayes/calls.vcf.gz"
    log:
        "results/logs/compress_vcf/log.log",
    params:
        extra="--write-index",
    wrapper:
        "v5.10.0/bio/bcftools/view"


def concat_gvcfs_input(w):
    with checkpoints.samtools_faidx.get(**w).output[0].open() as f:
        chroms = [l.split('\t')[0] for l in f.readlines()]
    return expand(
            "results/freebayes/variants/vcfs/{chrom}/variants.{i}.gvcf",
            chrom=chroms,
            i=range(1,config["freebayes"]["chunks"]+1)
            )


rule concat_gvcfs:
    input:
        calls = concat_gvcfs_input,
    output:
        vcf = temp("results/freebayes/calls.g.vcf"),
    log:
        "results/logs/concat_vcfs/log.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads: 4
    shell:
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"


rule compress_gvcf:
    input:
        rules.concat_gvcfs.output,
    output:
        "results/freebayes/calls.g.vcf.gz"
    log:
        "results/logs/compress_vcf/log.log",
    params:
        extra="--write-index",
    wrapper:
        "v5.10.0/bio/bcftools/view"
