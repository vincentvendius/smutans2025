#use megahit to assembly the raw reads into contigs using de brujin graph architecture
rule generate_contigs:
    params:
        outdir=ASSEMBLY_DIR + PREFIX,
        contiglen="1000",
        minlen="1000"
        samples=(','.join(expand(STREP_FQ_DIR + "Yana_old{sample_lib}.fq.gz",sample_lib=glob_wildcards(os.path.join(STREP_FQ_DIR,"Yana_old{sample_lib}.fq.gz")).sample_lib))),
        tmp=TMPPATH
        contigs=ASSEMBLY_DIR + PREFIX + "/final.contigs.fa",
        contigs_fxd=ASSEMBLY_DIR + PREFIX + "/FXD.contigs.fa"
        contigs_filtered=ASSEMBLY_DIR + PREFIX + "/filter.contigs.fa",
    threads: 80
    output:
        cdb=ASSEMBLY_DIR + PREFIX + "/contigs.db"
    shell:
        """
        rmdir {params.outdir}
        megahit -r {params.samples} -o {params.outdir} --min-contig-len {params.contiglen} -t {threads} --tmp-dir {params.tmp}
        echo Safety change of megahit contig names to be compatible with anvio
        anvi-script-reformat-fasta {params.contigs} --simplify-names -o {params.contigs_fxd} --prefix {PREFIX} 
        echo this rule is to filter out small contigs, but i have choosen a threshold which makes this step redundant at the moment
        anvi-script-reformat-fasta {params.contigs_fxd} --min-len {params.minlen} --simplify-names -o {params.contigs_filtered}
        echo generate a contig database of metadata, which i pad with contig profile information
        anvi-gen-contigs-database -f {params.contigs_filtered} -o {output.cdb} --num-threads {threads} -n {PREFIX}
        anvi-run-hmms -c {output.cdb} --num-threads {threads}
        anvi-run-ncbi-cogs -c {output.cdb} --num-threads {threads}
        """

