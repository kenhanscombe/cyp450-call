rule plink2_list_duplicates:
    """
    List duplicate variants to exclude.
    """
    input:
        "data/chr{chromosome}.pgen",
        "data/chr{chromosome}.pvar",
        "data/chr{chromosome}.psam"
    output: "data/chr{chromosome}.rmdup.list"
    log: "logs/plink2_list_duplicates_chr{chromosome}.log"
    params: "data/chr{chromosome}"
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    resources: mem_per_cpu='32G', time='0-0:30'
    shell:
        """
        plink2 \
        --pfile {params} \
        --rm-dup 'exclude-mismatch' 'list' \
        --out {params} 2> {log}

        if [ ! -e {params}.rmdup.list ]; then
            touch {params}.rmdup.list
        fi
        """


rule plink2_remove_duplicates:
    """
    Remove duplicate variants to resolve ambiguity in position update.
    """
    input: 
        "data/chr{chromosome}.pgen",
        "data/chr{chromosome}.pvar",
        "data/chr{chromosome}.psam",
        rmdup = "data/chr{chromosome}.rmdup.list"
    output:
        "data/chr{chromosome}_dedup.pgen",
        "data/chr{chromosome}_dedup.pvar",
        "data/chr{chromosome}_dedup.psam"
    log: "logs/plink2_remove_duplicates_chr{chromosome}.log"
    params: 
        pfile = lambda w, input: os.path.splitext(input[0])[0],
        pfile_dedup = lambda w, output: os.path.splitext(output[0])[0]
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    resources: mem_per_cpu='32G', time='0-0:30'
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
        --exclude {input.rmdup} \
        --make-pgen \
        --out {params.pfile_dedup} 2> {log}
        """
        # --export vcf-4.2 vcf-dosage=DS 'bgz' \


rule make_BED:
    """
    Make a UCSC BED format file to lift over.
    """
    input: "data/chr{chromosome}_dedup.pvar",
    output: "data/chr{chromosome}.BED"
    log: "logs/make_BED_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem='64G', time='0-1:30'
    shell:
        """
        awk 'NR>1 {{print "chr" $1, $2, $2 + 1, $3}}' {input} > {output} 2> {log}
        """
        # bcftools query -f '%CHROM %POS %ID\n' {input} > {params.tmp} 2> {log}
        # awk '{{print "chr" $1, $2, $2 + 1, $3}}' {params.tmp} > {output}
        # rm {params.tmp}


rule ucsc_liftover:
    """
    Lift imputed variant positions over from GRCh37 to GRCh38 to match
    exome data.
    """
    input: 
        bed = "data/chr{chromosome}.BED",
        chain = "data/hg19ToHg38.over.chain.gz"
    output:
        bed38 = "data/chr{chromosome}_grch38.BED",
        unlifted = "data/chr{chromosome}_unlifted.BED"
    log: "logs/ucsc_liftover_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem='64G', time='0-1:30'
    shell:
        """
        liftOver {input.bed} {input.chain} {output.bed38} {output.unlifted} &> {log}
        """


rule plink2_update_map:
    """
    Update variant position in pfiles
    """
    input:
        "data/chr{chromosome}_dedup.pgen",
        "data/chr{chromosome}_dedup.pvar",
        "data/chr{chromosome}_dedup.psam",
        bed = "data/chr{chromosome}_grch38.BED"
    output:
        "data/chr{chromosome}_grch38.pgen",
        "data/chr{chromosome}_grch38.pvar",
        "data/chr{chromosome}_grch38.psam"
    log: "logs/plink2_update_map_chr{chromosome}.log"
    params:
        pfile_in = lambda w, input: os.path.splitext(input[0])[0],
        pfile_out = lambda w, output: os.path.splitext(output[0])[0]
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    resources: mem='64G', time='0-1:30'
    shell:
        """
        plink2 \
        --pfile {params.pfile_in} \
        --update-map {input.bed} 2 4 \
        --sort-vars \
        --make-pgen \
        --out {params.pfile_out} 2> {log}
        """
    

rule greedyrelate_samples_to_remove:
    """
    Samples above KING relatedness cut-off to remove.
    """
    input: "data/chr{chromosome}_grch38.psam"
    output: "data/samples_relatedness_cutoff_chr{chromosome}.remove"
    log: "logs/greedyrelate_samples_to_remove_chr{chromosome}.log"
    params: chrom = "{chromosome}"
    conda: "../envs/ukbqc.yaml"
    envmodules: r
    resources: mem='64G', time='0-1:30'
    shell:
        f"""
        Rscript workflow/scripts/king_remove.R {prj_dir} {{input}} {{params.chrom}} &> {{log}}
        """


rule plink2_to_vcf:
    """
    Convert to VCF, removing relateds
    """
    input: 
        pgen = "data/chr{chromosome}_grch38.pgen",
        pvar = "data/chr{chromosome}_grch38.pvar",
        psam = "data/chr{chromosome}_grch38.psam",
        rel = "data/samples_relatedness_cutoff_chr{chromosome}.remove"
    output:
        vcf = "data/chr{chromosome}.vcf.gz"
    params:
        pfile = lambda w, input: os.path.splitext(input.pgen)[0],
        out = lambda w, output: os.path.splitext(os.path.splitext(output.vcf)[0])[0]
    log: "logs/plink2_to_vcf_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    resources: mem='64G', time='0-1'
    shell:
        """
        plink2 \
        --pfile {params.pfile} \
        --remove {input.rel} \
        --export vcf-4.2 vcf-dosage=DS 'bgz' \
        --out {params.out} 2> {log}
        """