#once a contig database is completed it is time to align the samples against the contigs
#we will use bowtie2, which requires us to index the contigs first
rule index:
    input:
        cdb=ASSEMBLY_DIR + PREFIX + "/contigs.db",
        fa=ASSEMBLY_DIR + PREFIX + "/filter.contigs.fa"
    output:
        multiext("data/ref_index/" + PREFIX,".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2")
    params:
        ref_name="data/ref_index/" + PREFIX,
        query=STREP_FQ_DIR + "{sample}.*.trim.fq.gz",
        ref_name="data/ref_index/" + PREFIX
        profile=ASSEMBLY_DIR + PREFIX + "/{sample}",
        minlen="1000"
    threads:80
    shell:
        """
        bowtie2-build  --threads {threads} {input.fa} {params.ref_name}
        bowtie2 -D 20 -R 3 -N 1 -L 20 -i S,1,0.50 -p {threads} -x {params.ref_name} -U {params.query} | samtools view -F4 -bh > {output.bam}.raw
        samtools sort {output.bam}.raw -o {output.bam} -m 5000000000
        samtools index {output.bam}
        rmdir {params.profile}
        anvi-profile -c {input.cdb} -i {input.bam} --min-contig-length {params.minlen} --num-threads {threads} -o {params.profile} --sample-name {wildcards.sample} --min-coverage-for-variability 5 --cluster-contigs
        """



#Alignment of sample reads against sample set contig database
rule align:
    input:
        ref = multiext("data/ref_index/" + PREFIX,".1.bt2",".2.bt2",".3.bt2",".4.bt2",".rev.1.bt2",".rev.2.bt2")
    params:
    output:
        bam="data/bam/{sample}.bam"
    threads: 40
    shell:

#profile the bam files
rule sample_profile:
    input:
        cdb=ASSEMBLY_DIR + PREFIX + "/contigs.db",
        bam="data/bam/{sample}.bam",
    output:
        profile=ASSEMBLY_DIR + PREFIX + "/{sample}/PROFILE.db"
    params:
    threads: 40
    shell:
        """
        """

