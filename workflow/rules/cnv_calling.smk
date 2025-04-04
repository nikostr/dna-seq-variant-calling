#    Call CNVs for each sample and optionally refine breakpoints using delly SV calls
# delly cnv -o c1.bcf -g hg19.fa -m hg19.map -l delly.sv.bcf input.bam
rule delly_cnv:
    input:
        ref=config["genome"],
        sv_calls=rules.delly_genotype.output.genotype,
        alns=lambda w: set(
            [
                f
                for f in read_mapping.get_collect_bams_input(w)
                if f"{w.sample}.bam" in f
            ]
        ),
        mappability_map=rules.convert_mappability_map.output.fa,
    output:
        genotype="results/delly/individual_cnv_calls/{sample}.bcf",
    conda:
        "../envs/delly.yaml"
    log:
        "results/logs/individual_cnv_calls/{sample}.log",
    shell:
        "delly cnv "
        "--outfile {output.genotype} "
        "--genome {input.ref} "
        "--mappability {input.mappability_map} "
        "--svfile {input.sv_calls} "
        "{input.alns} "
        "2> {log} "


#    Merge CNVs into a unified site list
# delly merge -e -p -o sites.bcf -m 1000 -n 100000 c1.bcf c2.bcf ... cN.bcf
rule delly_merge_cnvs:
    input:
        bcfs=expand(
            "results/delly/individual_cnv_calls/{sample}.bcf",
            sample=samples.sample_id.unique(),
        ),
    output:
        sites="results/delly/delly_cnv_sites/delly_sites.bcf",
    conda:
        "../envs/delly.yaml"
    log:
        "results/logs/delly_cnv_sites/delly_cnv_sites.log",
    shell:
        "delly merge "
        "--cnvmode "
        "--pass "
        "--outfile {output.sites} "
        "--minsize 1000 "
        "--maxsize 100000 "
        "{input.bcfs} "
        "2> {log} "


#    Genotype CNVs for each sample
# delly cnv -u -v sites.bcf -g hg19.fa -m hg19.map -o geno1.bcf input.bam
rule delly_genotype_cnvs:
    input:
        ref=config["genome"],
        alns=lambda w: set(
            [
                f
                for f in read_mapping.get_collect_bams_input(w)
                if f"{w.sample}.bam" in f
            ]
        ),
        sites=rules.delly_merge_cnvs.output.sites,
        mappability_map=rules.convert_mappability_map.output.fa,
    output:
        genotype="results/delly/cnv_genotype/{sample}.bcf",
    conda:
        "../envs/delly.yaml"
    log:
        "results/logs/cnv_genotype/{sample}.log",
    shell:
        "delly cnv "
        "--segmentation "
        "--vcffile {input.sites} "
        "--genome {input.ref} "
        "--mappability {input.mappability_map} "
        "--outfile {output.genotype} "
        "{input.alns} "
        "2> {log} "


rule delly_bcftools_cnv_index:
    input:
        rules.delly_genotype_cnvs.output.genotype,
    output:
        "results/delly/cnv_genotype/{sample}.bcf.csi",
    log:
        "results/logs/delly_bcftools_cnv_index/{sample}.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v3.3.3/bio/bcftools/index"


#    Merge genotypes using bcftools
# bcftools merge -m id -O b -o merged.bcf geno1.bcf ... genoN.bcf
rule delly_bcftools_merge_cnv:
    input:
        calls=expand(
            "results/delly/cnv_genotype/{sample}.bcf",
            sample=samples.sample_id.unique(),
        ),
        idx=expand(
            "results/delly/cnv_genotype/{sample}.bcf.csi",
            sample=samples.sample_id.unique(),
        ),
    output:
        "results/delly/merged_cnv_genotype/merged_cnv_genotype.unfiltered.bcf",
    log:
        "results/logs/delly_bcftools_merge_cnv/delly_bcftools_merge_cnv.log",
    params:
        uncompressed_bcf=False,
        extra="",  # optional parameters for bcftools concat (except -o)
    wrapper:
        "v3.3.3/bio/bcftools/merge"


rule delly_bcftools_merged_cnv_index:
    input:
        rules.delly_bcftools_merge_cnv.output,
    output:
        "results/delly/merged_cnv_genotype/merged_cnv_genotype.unfiltered.bcf.csi",
    log:
        "results/logs/delly_bcftools_merged_cnv_index/delly_bcftools_merged_cnv_index.log",
    params:
        extra="",  # optional parameters for bcftools index
    wrapper:
        "v3.3.3/bio/bcftools/index"


#    Filter for germline CNVs
# delly classify -f germline -o filtered.bcf merged.bcf
rule delly_filter_germline_cnvs:
    input:
        bcf=rules.delly_bcftools_merge_cnv.output,
        idx=rules.delly_bcftools_merged_cnv_index.output,
    output:
        filtered="results/delly/merged_cnv_genotype/merged_cnv_genotype.bcf",
    conda:
        "../envs/delly.yaml"
    log:
        "results/logs/delly_filter_germline_cnvs/delly_filter_germline_cnvs.log",
    shell:
        "delly classify "
        "--filter germline "
        "--outfile {output.filtered} "
        "{input.bcf} "
        "2> {log} "
