import subprocess as sp


checkpoint samtools_faidx:
    input:
        config["genome"],
    output:
        config["genome"] + ".fai",
    log:
        "results/logs/samtools_faidx/index.log",
    wrapper:
        "v3.7.0/bio/samtools/faidx"


#rule freebayes:
#    input:
#        alns=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
#        idxs=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f],
#        ref=config["genome"],
#        ref_idx=rules.samtools_faidx.output,
#    output:
#        vcf = "results/freebayes/calls.vcf",
#    log:
#        "results/logs/freebayes/calls.log",
#    threads: 16
#    resources:
#        mem_mb=1024,
#    wrapper:
#        "v3.7.0/bio/freebayes"

rule GenerateFreebayesRegions:
    input:
        ref_idx=rules.samtools_faidx.output,
        #ref_idx = reference,
        #index = reference + ".fai",
        #bams = expand("resources/alignments/{sample}.bam", sample=samples)
    output:
        directory("results/freebayes/regions")
        #regions = expand(
        #        "results/freebayes/regions/genome.{chrom}.region.{i}.bed",
        #        chrom=chroms,
        #        chrom=chroms,
        #        i=range(1,config["freebayes"]["chunks"]+1)
        #        )
    log:
        "results/logs/freebayes/GenerateFreebayesRegions.log"
    params:
        chunks = config["freebayes"]["chunks"],
    #conda:
    #    "../envs/freebayes-env.yaml"
    script:
        "../scripts/fasta_generate_regions.py"


rule generate_bam_list:
    input:
        bams=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
    output:
        temp("results/freebayes/bams.txt")
    run:
        with open(output[0], 'w') as f:
            f.writelines([l + '\n' for l in input.bams])
        

rule VariantCallingFreebayes:
    input:
        bams=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f]),
        index=lambda w: [f + '.bai' for f in read_mapping.get_collect_bams_input(w) if 'illumina' in f],
        ref=config["genome"],
        samples=rules.generate_bam_list.output,
        all_beds=rules.GenerateFreebayesRegions.output,
        #bams = expand("resources/alignments/{sample}.bam", sample=samples),
        #index = expand("resources/alignments/{sample}.bam.bai", sample=samples),
        #ref = reference,
        #samples = bamlist,
        #regions = "results/freebayes/regions/genome.{chrom}.region.{i}.bed"
    output:
        temp("results/freebayes/variants/vcfs/{chrom}/variants.{i}.vcf")
    params:
        regions="results/freebayes/regions/genome.{chrom}.region.{i}.bed"
    log:
        "results/logs/freebayes/VariantCallingFreebayes/{chrom}.{i}.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads:1
    shell: "freebayes -f {input.ref} -t {params.regions} -L {input.samples} > {output} 2> {log}"


def concat_vcfs_input(w):
    with checkpoints.samtools_faidx.get(**w).output[0].open() as f:
        chroms = [l.split('\t')[0] for l in f.readlines()]
    return expand(
            "results/freebayes/variants/vcfs/{chrom}/variants.{i}.vcf",
            chrom=chroms,
            i=range(1,config["freebayes"]["chunks"]+1)
            )


rule ConcatVCFs:
    input:
        calls = concat_vcfs_input,
    output:
        #"results/freebayes/variants/vcfs/variants.{chrom}.vcf"
        vcf = "results/freebayes/calls.vcf",
    log:
        "results/logs/ConcatVCFs/log.log"
    conda:
        "../envs/freebayes-env.yaml"
    threads: 4
    shell:  
        "bcftools concat {input.calls} | vcfuniq > {output} 2> {log}"
