#     SV calling is done by sample for high-coverage genomes or in small batches for low-coverage genomes
# delly call -g hg19.fa -o s1.bcf -x hg19.excl sample1.bam
rule delly_bcf:
    input:
        ref=config["genome"],
        ref_idx=rules.samtools_faidx.output,
        alns=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if w.sample + '.bam' in f]),
        idx=lambda w: set([f + '.bai' for f in read_mapping.get_collect_bams_input(w) if w.sample + '.bam' in f]),
    output:
        "results/delly/individual_calls/{sample}.bcf",
    params:
        extra="",  # optional parameters for delly (except -g, -x)
    log:
        "results/logs/delly_bcf/{sample}.log",
    threads: 1  # It is best to use as many threads as samples
    wrapper:
        "v5.10.0/bio/delly"



#     Merge SV sites into a unified site list
# delly merge -o sites.bcf s1.bcf s2.bcf ... sN.bcf
rule delly_merge:
    input:
        bcfs=expand(
            "results/delly/individual_calls/{sample}.bcf",
            sample=samples.sample_id.unique()
        ),
    output:
        sites="results/delly/delly_sites/delly_sites.bcf"
    conda: "../envs/delly.yaml"
    shell:
        "delly merge "
        "--outfile {output.sites} "
        "{input.bcfs} "



#     Genotype this merged SV site list across all samples. This can be run in parallel for each sample.
# delly call -g hg19.fa -v sites.bcf -o s1.geno.bcf -x hg19.excl s1.bam
# delly call -g hg19.fa -v sites.bcf -o sN.geno.bcf -x hg19.excl sN.bam
rule delly_genotype:
    input:
        ref=config["genome"],
        sites=rules.delly_merge.output.sites,
        alns=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if w.sample + '.bam' in f]),
    output:
        genotype="results/delly/genotype/{sample}.bcf",
    conda: "../envs/delly.yaml"
    resources:
        mem_mb=9*1024 # See https://github.com/dellytools/delly/issues/172
    shell:
        "delly call "
        "--genome {input.ref} "
        "--vcffile {input.sites} "
        "--outfile {output.genotype} "
        "{input.alns} "


rule delly_bcftools_index:
    input:
        rules.delly_genotype.output.genotype,
    output:
        "results/delly/genotype/{sample}.bcf.csi",
    log:
        "results/logs/delly_bcftools_index/{sample}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v5.10.0/bio/bcftools/index"


#     Merge all genotyped samples to get a single VCF/BCF using bcftools merge
# bcftools merge -m id -O b -o merged.bcf s1.geno.bcf s2.geno.bcf ... sN.geno.bcf
rule delly_bcftools_merge:
    input:
        calls=expand(
            "results/delly/genotype/{sample}.bcf",
            sample=samples.sample_id.unique()
        ),
        idx=expand(
            "results/delly/genotype/{sample}.bcf.csi",
            sample=samples.sample_id.unique()
        ),
    output:
        "results/delly/merged_genotype/merged_genotype.bcf",
    log:
        "results/logs/delly_bcftools_merge/delly_bcftools_merge.log",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v5.10.0/bio/bcftools/merge"


rule delly_copy_number_segmentation:
    input:
        ref=config["genome"],
        ref_idx=rules.samtools_faidx.output,
        alns=lambda w: set([f for f in read_mapping.get_collect_bams_input(w) if w.sample + '.bam' in f]),
        idx=lambda w: set([f + '.bai' for f in read_mapping.get_collect_bams_input(w) if w.sample + '.bam' in f]),
        mappability_map=config['delly']['mappability_map'],
    output:
        covfile="results/delly/read_depth_profile/{sample}.cov.gz",
        vcf="results/delly/read_depth_profile/{sample}.vcf",
    conda: "../envs/delly.yaml"
    log:
        "results/logs/read_depth_profile/{sample}.log",
    shell:
        "delly cnv "
        "--adaptive-windowing "
        "--segmentation "
        "--genome {input.ref} "
        "--mappability {input.mappability_map} "
        "--covfile {output.covfile} "
        "{input.alns} "
        "> {output.vcf} "
        "2> {log} "


rule segmentation_bed:
    input:
        vcf=rules.delly_copy_number_segmentation.output.vcf,
    output:
        bed="results/delly/read_depth_profile/{sample}.bed",
    conda: "../envs/bcftools.yaml"
    log:
        "results/logs/segmentation_bed/{sample}.log",
    shell:
        "bcftools query "
        '--format "%CHROM\t%POS\t%INFO/END\t%ID[\t%RDCN]\n" '
        "{input.vcf} "
        "> {output.bed} "
