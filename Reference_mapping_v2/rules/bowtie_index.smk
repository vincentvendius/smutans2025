#index all the reference genomes
rule index:
    input:
        ref=(
            lambda wildcards: REFTABLE.fa[wildcards.refid]
        ),
    output:
        multiext(config["bowtie2_tp_db"] + "/{refid}",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2")
    params:
        ref_name=config["bowtie2_tp_db"] + "/{refid}"
    threads: 35
    conda:
        "../envs/bam-filter.yaml"
    shell:
        "bowtie2-build  --threads {threads} {input.ref} {params.ref_name}"


