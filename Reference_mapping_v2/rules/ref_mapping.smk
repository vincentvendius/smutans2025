

#samplelist=list(sample_table_read.label_md5["SII"])

def get_fq(wildcards):
    samples=sample_table_read.loc[(wildcards.smp)].label_md5
    if type(samples)==str:
        samplelist=[samples]
    else:   
        samplelist=samples.tolist()

    #print(samplelist)
    pathlist=["/science/willerslev/users-shared/science-snm-willerslev-xvh856/pipelines/aMAW-taxonomy/results/read-noneuk/" + lib + ".read-noneuk.fq.gz" for lib in samplelist]
    #print(pathlist)
    return(pathlist)



rule ref_mapping:
    input:
        fastq=(
            lambda wildcards: SAMPLETABLE.collapsed[wildcards.smp]
        ),
        #fastq=(
        #    lambda wildcards: SAMPLETABLE.fq[wildcards.smp]
        #),
        #fastq=get_fq,
        #print(fastq)
        #fastq=config["rdir"] + "/read-noneuk/{smp}.read-noneuk.fq.gz",
        #fastq= config["strep_fq_dir"] + "NEO137" + config["strep_fq_file"],
        #stats=config["rdir"] + "/stats/{smp}.stats-initial.txt",
        ref = multiext(config["bowtie2_tp_db"] + "/{refid}",".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2"),
    output:
        # filter-bam output 
        bowtie2_bam_dedup_filt=(
            config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.dedup.filtered.bam"
        ),
        bowtie2_bam_dedup_metrics=(
            config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.dedup.metrics"
        ),
        gc=config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.srt.filter.genomecov",
        cov=config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.srt.filter.coverage.txt"

    threads: config["tp_threads"]
    params:
        pdir=config["pdir"] + "/bam",
        wdir=config["wdir"],
        bowtie2_bin=config["bowtie2_bin"],
        samtools_bin=config["samtools_bin"],
        samtools_view_parms=config["samtools_view_parms"],
        samtools_sort_parms=config["samtools_sort_parms"],
        picard_bin=config["picard_bin"],
        picard_jar=config["picard_jar"],
        java_opts=config["picard_java_opts"],
        tp_bowtie2_parms=config["tp_bowtie2_parms"],
        bowtie2_tp_db=config["bowtie2_tp_db"] + "/{refid}",
        label="{smp}",
        bowtie2_bam_tmp=config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.tmp.bam",
        bowtie2_bam_dedup_tmp=config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.dedup.tmp.bam",
        bowtie2_sam_tmp=config["pdir"] + "/bam/{smp}/" + config["prefix"] + "/{refid}.dedup.filtered.sam",
    log:
        config["pdir"] + "/logs/bam/{smp}/{refid}.bam.log",
    benchmark:
        config["pdir"] + "/benchmarks/bam/{smp}/{refid}.bam.bmk"
    conda:
        "../envs/bam-filter.yaml"
    message:
        """--- High-res taxonomic profiling."""
    shell:
        """
        set -x
        cd {params.pdir} || {{ echo "Cannot change dir"; exit 1; }}
        seqkit seq -j {threads} {input.fastq} | \
            {params.bowtie2_bin} -x {params.bowtie2_tp_db} \
            -p {threads} \
            -q -U - \
            {params.tp_bowtie2_parms} \
            | {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} > {params.bowtie2_bam_tmp}
        # Mark Duplicates
        {params.picard_bin} {params.java_opts} \
            -jar {params.picard_jar} MarkDuplicates \
            I={params.bowtie2_bam_tmp} \
            O={params.bowtie2_bam_dedup_tmp} \
            M={output.bowtie2_bam_dedup_metrics} \
            >> {log} 2>&1
        # Sort, filter and index 
        {params.samtools_bin} sort -@ {threads} {params.samtools_sort_parms} {params.bowtie2_bam_dedup_tmp} \
                | {params.samtools_bin} view {params.samtools_view_parms} > {output.bowtie2_bam_dedup_filt}
        {params.samtools_bin} index -@ {threads} {output.bowtie2_bam_dedup_filt}
        # Update timestamp on indexed file from samtools to prevent erroneous reruns in snakemake
        touch -c {output.bowtie2_bam_dedup_filt}.bai
        touch {output.gc} 
        (bedtools genomecov -ibam {output.bowtie2_bam_dedup_filt} > {output.gc}) || true 
        cat {output.gc} | mawk '{{print $1"\\t"$2*$3"\\t"$4}}' | datamash -g1 sum 2 first 3 | mawk '{{print $1,$2/$3}}' > {output.cov}
        rm -rf {params.bowtie2_bam_dedup_tmp} {params.bowtie2_bam_tmp}
        cd {params.wdir} || {{ echo "Cannot change dir"; exit 1; }}
        """
