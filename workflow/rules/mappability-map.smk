rule dicey_chop:
    input:
        ref=config["genome"],
    output:
        r1=temp("results/delly/mappability-map/read1.fq.gz"),
        r2=temp("results/delly/mappability-map/read2.fq.gz"),
    params:
        lambda w, output: f"--fq1 {output.r1.removesuffix('.fq.gz')} --fq2 {output.r2.removesuffix('.fq.gz')}",
    log:
        "results/logs/delly/mappability-map/dicey-chop.log",
    conda:
        "../envs/dicey.yaml"
    shell:
        "dicey chop {params} {input.ref} &> {log}"


rule map_dicey:
    input:
        reads=rules.dicey_chop.output,
        reference=config["genome"],
        idx=multiext(
            config["genome"],
            ".0123",
            ".amb",
            ".ann",
            ".pac",
            ".pos_packed",
            ".suffixarray_uint64",
            ".suffixarray_uint64_L0_PARAMETERS",
            ".suffixarray_uint64_L1_PARAMETERS",
            ".suffixarray_uint64_L2_PARAMETERS",
        ),
    output:
        temp("results/delly/mappability-map/srt.bam"),
    log:
        "results/logs/delly/mappability-map/mapping.log",
    params:
        bwa="bwa-meme",
        extra="-M",
        sort="samtools",  # Can be 'none' or 'samtools or picard'.
        sort_order="coordinate",  # Can be 'coordinate' (default) or 'queryname'.
        sort_extra="",  # Extra args for samtools.
        dedup="mark",  # Can be 'none' (default), 'mark' or 'remove'.
        dedup_extra="-M",  # Extra args for samblaster.
        exceed_thread_limit=True,  # Set threads als for samtools sort / view (total used CPU may exceed threads!)
        embed_ref=True,  # Embed reference when writing cram.
    threads: 16
    wrapper:
        "v3.3.6/bio/bwa-memx/mem"


rule samtools_index_temp:
    input:
        rules.map_dicey.output[0],
    output:
        temp(rules.map_dicey.output[0] + ".bai"),
    log:
        "results/logs/samtools_index/dicey.log",
    threads: 8  # This value - 1 will be sent to -@
    wrapper:
        "v3.3.6/bio/samtools/index"


rule dicey_mappability:
    input:
        bam=rules.map_dicey.output[0],
        bai=rules.samtools_index_temp.output,
    output:
        tmp=temp("results/delly/mappability-map/tmp.fa.gz"),
    conda:
        "../envs/dicey.yaml"
    log:
        "results/logs/delly/mappability-map/dicey-mappability.log",
    shell:
        "dicey mappability2 "
        "--outfile {output} "
        "{input.bam} "
        "&> {log}"


rule convert_mappability_map:
    input:
        tmpmap=rules.dicey_mappability.output.tmp,
    output:
        fa="results/delly/mappability-map/map.fa.gz",
        fai="results/delly/mappability-map/map.fa.gz.fai",
        gzi="results/delly/mappability-map/map.fa.gz.gzi",
    params:
        tmp_in=lambda w, input: os.path.splitext(input.tmpmap)[0],
        tmp_out=lambda w, output: os.path.splitext(output.fa)[0],
    conda:
        "../envs/samtools.yaml"
    log:
        "results/logs/delly/mappability-map/convert-mappability-map.log",
    shell:
        "gunzip {input} "
        "2> {log} "
        "&& mv {params.tmp_in} {params.tmp_out} "
        "2>> {log} "
        "&& bgzip {params.tmp_out} "
        "2>> {log} "
        "&& samtools faidx {output.fa} "
        "2>> {log} "
