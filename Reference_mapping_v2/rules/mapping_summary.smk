rule mapping_summary:
    input:
        bowtie2_bam_dedup_filt=(
            expand(config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.dedup.filtered.bam",refid=REFID,allow_missing=True)
        ),
        bowtie2_bam_dedup_metrics=(
            expand(config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.dedup.metrics",refid=REFID,allow_missing=True)
        ),
        gc=expand(config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.srt.filter.genomecov",refid=REFID,allow_missing=True),
        cov=expand(config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/{refid}.srt.filter.coverage.txt",refid=REFID,allow_missing=True),
    output: 
        summary_tsv=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.summary.tsv.gz",
        plot=config["pdir"] + "/plots/" + config["prefix"] + "/{smp}.cov.pdf",
        bestmap=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmapping.tsv",
        #pdf=config["pdir"] + "/plots/" + config["prefix"] + "/{smp}.editDist.pdf",
        #editdist_tsv=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.editDist.tsv",
    params:
        Rbin=config["Rbin"],
        editdistdir=config["editdistdir"],
        summarydir=config["summarydir"],
        prefix=config["prefix"],
        libdb=config["libdb"],
        bam_dir=config["pdir"] + "/",
        dist_matrix=config["dist_matrix"],
        label="{smp}",
    conda:
        "../envs/r.yaml"
    threads: 24
    shell:
        """
        {params.Rbin} {params.summarydir} {params.prefix} {params.libdb} {params.label} {threads} {params.bam_dir} {output.plot} {params.dist_matrix} {output.bestmap} {output.summary_tsv}
        """
        #{params.Rbin} {params.editdistdir} {params.prefix} {params.libdb} {output.pdf} {output.editdist_tsv} {params.bam_dir} {params.label}
        


        
