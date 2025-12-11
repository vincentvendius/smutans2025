rule genotyping:
    input:
        bestmap=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmapping.tsv",
    output:
        hetsum=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmap.het.tsv",
        genes=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmap.genes.tsv",
        depth=multiext(config["pdir"] + "/stats/" + config["prefix"] + "/mosdepth/{smp}",".per-base.bed.gz",".mosdepth.global.dist.txt",".mosdepth.summary.txt",".per-base.bed.gz.csi"),
        gtf_filter=config["pdir"] + "/stats/" + config["prefix"] + "/{smp}.gtf.filter",
        consensus=config["pdir"] + "/fasta/" + config["prefix"] + "/{smp}.fa",
    params:
        pdir=config["pdir"],
        #refpath=config["ref_dir"] + "*/Streptococcus_*",
        #reffile="*" + config["ref_file"],
        #refpath="/projects/sikora/people/xvh856/projects/mutans_assemblies/bacteria/",
        refpath="/maps/projects/lundbeck/data/dbs/krakenuniq/hum_microbe/library/bacteria/*/*",
        reffile="*.fna",
        bampath=config["pdir"] + "/bam/" + config["prefix"] + "/{smp}/",
        bamfile=".dedup.filtered.bam",
        vcf=config["pdir"] + "/stats/" + config["prefix"] + "/{smp}.vcf.gz",
        mosout=config["pdir"] + "/stats/" + config["prefix"] + "/mosdepth/{smp}",
        depth=config["pdir"] + "/stats/" + config["prefix"] + "/mosdepth/{smp}.per-base.bed.gz",
        gtf=config["pdir"] + "/stats/" + config["prefix"] + "/{smp}.gtf.gz",
        genetab=config["pdir"] + "/tables/" + config["prefix"] + "/{smp}.bestmap.genes.tmp.tsv",
        calls_norm=config["pdir"] + "/stats/" + config["prefix"] + "/{smp}.vcf.norm.gz",
        calls_filt=config["pdir"] + "/stats/" + config["prefix"] + "/{smp}.vcf.filt.gz",
    conda:
        "../envs/genotyping.yaml"
    shell:
        """
        set -x
        rm -rf {params.gtf} {params.mosout} {params.vcf} {params.calls_norm} {params.calls_filt} 
        cd {params.pdir} || {{ echo "Cannot change dir"; exit 1; }}
        #get id for best reference genome for sample
        BEST_REF=$(cat {input.bestmap} | cut -f2)

        #compute variant call for the best mapping reference
        bcftools mpileup -q20 -a FORMAT/AD,FORMAT/DP -f \
        {params.refpath}${{BEST_REF}}{params.reffile} \
        {params.bampath}${{BEST_REF}}{params.bamfile} -Ou \
        | bcftools call -m -A -a FORMAT/GQ -Oz > {params.vcf} \

        #Extract substitutions from the reference genome variant call file, and calculate the heterogenity rate using AWK as (0/0+0/1+1/1)/(0/1)
        bcftools view -i'FORMAT/DP >= 5' -V indels -M2 {params.vcf} \
        | bcftools query -f'[%GT]\n' -e '(REF=\"C\" & ALT=\"T\") | \
        (REF=\"T\" & ALT=\"C\") | (REF=\"G\" & ALT=\"A\") | \
        (REF=\"A\" & ALT=\"G\")' | sort | uniq -c | \
        awk '{{sum+=$1}} FNR==2 {{het=$1/sum}} END{{print het;}}' | paste {input.bestmap} - > {output.hetsum}

        #calculate depth of mapping at different positions
        mosdepth {params.mosout} {params.bampath}${{BEST_REF}}{params.bamfile} 

        #Use wget and entrez tools to retrieve metadata for the reference genomes we are interested in
        wget -r -A gtf.gz -nc -O {params.gtf} \"$(esearch -db \"assembly\" -query ${{BEST_REF}} \
        | esummary | xtract -pattern DocumentSummary -element FtpPath_RefSeq)\"
        cat {params.gtf} | zgrep -P \"\\tRefSeq\\t\" | cut -f1,4,5,6,9 | \
        cut --delimiter=\";\" -f1,4 | sed \'s/ \"/\\t/g\' | sed \'s/\";/\\t/g\' | \
        sed \'s/\"//g\' > {output.gtf_filter}
        rm -rf {params.pdir}/ftp.ncbi.nlm.nih.gov

        #get statistics about coverage of different genes on our sample
        bedtools intersect -a {output.gtf_filter} -b {params.depth} -wa -wb | \
        bedtools groupby -g 6,8,2,3 -c 12 -o mean,min,max | \
        awk \'{{print $0, \"\\t{wildcards.smp}\"}}' > {params.genetab}
        zcat {params.depth} | grep -P -v "\\t0$" | \
        bedtools coverage -a {output.gtf_filter} -b - | cut -f 12 | \
        paste {params.genetab} - > {output.genes} 
        

        #create consensus genomes
        bcftools index {params.vcf}
        bcftools norm -f {params.refpath}${{BEST_REF}}{params.reffile} {params.vcf} -Ob -o {params.calls_norm}
        bcftools filter --IndelGap 5 {params.calls_norm} -Ob -o {params.calls_filt}
        bcftools index {params.calls_filt}
        # apply variants to create consensus sequence
        cat {params.refpath}${{BEST_REF}}{params.reffile} | bcftools consensus {params.calls_filt} > {output.consensus}
        """

