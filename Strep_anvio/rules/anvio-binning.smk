#profile the bam files
rule sample_profile:
    input:
        cdb=ASSEMBLY_DIR + PREFIX + "/contigs.db",
        bam="data/bam/{sample}.bam",
        profile=expand(ASSEMBLY_DIR + PREFIX + "/{sample}/PROFILE.db",sample=[i.split("\t",3)[3] for i in open(SAMPLE_INFO).read().splitlines()]),
    params:
        profile=ASSEMBLY_DIR + PREFIX + "/{sample}",
        minlen="1000",
        merged=ASSEMBLY_DIR + PREFIX + "/PROFILE_MERGED",
        cluster_algo=BINNING_ALGO,
        collection_name=COLLECTION_NAME,
        out_dir=ASSEMBLY_DIR + PREFIX + "/CLUSTER-SUMMARY",
    threads: 40
    output:
        out_dir=ASSEMBLY_DIR + PREFIX + "/CLUSTER-SUMMARY/bins_summary.txt",
    shell:
        """
        rmdir {params.profile}
        anvi-profile -c {input.cdb} -i {input.bam} --min-contig-length {params.minlen} --num-threads {threads} -o {params.profile} --sample-name {wildcards.sample} --min-coverage-for-variability 5 --cluster-contigs
        echo merge the bam profiles for all samples in the set and do automatic binning
        rmdir {params.merged}
        anvi-merge {input.profile} -o {params.merged} -c {input.cdb}
        anvi-cluster-contigs -p {output.merged} -c {input.cdb} --num-threads {threads} -C {params.collection_name} --driver {params.cluster_algo} --just-do-it
        echo summarize the automatic binning cluster
        rmdir {params.out_dir} 
        anvi-summarize -c {input.cdb} -p {input.merged} -C {params.collection_name} -o {params.out_dir}
        """

