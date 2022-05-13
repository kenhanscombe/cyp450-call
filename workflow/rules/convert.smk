
gene_hg19 = config["gene_hg19"]
gene_bgenix_range = {}

def bgenix_range(chr, pos, window=1e6):
    range_list = [str(int(sum(p))) for p in zip(pos, [-window, window])]
    return f"{chr}:" + "-".join(range_list)

# NB. This only works because there is one gene per chr!
for g in gene_hg19:
    chrom = gene_hg19[g][0]
    position = gene_hg19[g][1:]
    gene_bgenix_range[chrom] = bgenix_range(chrom, position)


rule bgenix_get_region:
    """Extract 1MB up- and downstream of pharmacogene"""
    input: imp_dir / "ukb_imp_chr{chromosome}_v3.bgen"
    output: "data/cyp_region_chr{chromosome}.bgen"
    params: incl_range = lambda w: gene_bgenix_range[int(w.chromosome)]
    log: "logs/bgenix_get_region_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    # threads: 4
    # resources: mem='128G', time='0-1', nodes=1, ntasks=4
    resources: mem='64G', time='0-30'
    shell:
        f"""
        {bgenix} -g {{input}} -incl-range {{params.incl_range}} > {{output}} 2> {{log}}
        """


rule sample_exclusions:
    """
    Write gender mismatch, heterozygosity and missingness outliers,
    putative sex chromosome aneuploidy, and excess relatives list for
    removal
    """
    input:
        sqc = imputed / sqc_file,
        fam = genotyped / fam_file,
        vcf = "data/ukb23156_c22_b16_v1.vcf.gz"
    output: 
        excl = "results/samples_fail_ukb_qc.remove",
        keep = "results/samples_with_exome_data.keep"
    params: tmp = "results/temp.tmp"
    conda: "../envs/ukbqc.yaml"
    envmodules: python3, vcftools
    log: "logs/sample_exclusions.log"
    shell:
        """
        vcf-query -l {input.vcf} > {output.keep}
        awk '{{print $1, $1}}' {output.keep} > {params.tmp} && \
        mv {params.tmp} {output.keep}

        python workflow/scripts/sqc_exclusions.py \
        --fam {input.fam} \
        --sqc {input.sqc} \
        --out {output.excl}
        """

# Use: ukb_imp_chr*_v3_MAF0_INFO7.bgen INFO=0.7 (no MAF filter)
# rule variant_exclusions:
#     """
#     Write MAF<0.0001 and INFO<0.5 variant list for exclusion
#     """
#     input:
#         info = imputed / "ukb_mfi_chr{chromosome}.txt"
#     output:
#         excl = "results/variants_fail_maf_info_chr{chromosome}.exclude"
#     conda: "../envs/ukbqc.yaml"
#     log: "logs/variant_excl_chr{chromosome}.log"
#     shell:
#         """
#         awk '{{if($6<0.0001 || $8<0.5) print $2;}}' \
#         {input.info} > {output.excl}
#         """


rule plink2_convert:
    """
    Convert imputed bgen format files to PLINK binary format. Keep
    samples with exome data, remove gender mismatches, heterozygosity
    and missingness outliers, putative sex chromosome aneuploidy, excess
    relatives. Exclude variants with maf<0.01, info<0.7, hwe p <1e-15,
    or missing in greater than 10% of individuals. Keep only bi-allelic
    variants.

    See chrchang comment https://www.biostars.org/p/383652/
    """
    input:
        # gen_file = imputed / "ukb_imp_chr{chromosome}.bgen",
        # gen_file = imp_dir / "ukb_imp_chr{chromosome}_v3_MAF0_INFO7.bgen",
        gen_file = "data/cyp_region_chr{chromosome}.bgen",
        sample_file = imputed / sample_file,
        sample_keep =  "results/samples_with_exome_data.keep",
        sample_excl = "results/samples_fail_ukb_qc.remove"
    output:
        pgen = "data/chr{chromosome}.pgen",
        pvar = "data/chr{chromosome}.pvar",
        psam = "data/chr{chromosome}.psam"
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    params:
        out = lambda w, output: os.path.splitext(os.path.splitext(output[0])[0])[0]
    log: "logs/plink2_convert_chr{chromosome}.log"
    resources: mem='32G', time='0-1'
    # resources: mem='32G', time='0-1', nodes=1, ntasks=1, cpus_per_task=4
    shell:
        """
        plink2 \
        --bgen {input.gen_file} 'ref-first' \
        --sample {input.sample_file} \
        --keep {input.sample_keep} \
        --remove {input.sample_excl} \
        --hwe 1e-25 \
        --make-pgen \
        --out {params.out} 2> {log}
        """
        # --geno 0.1 \
        # --max-alleles 2 \
        # snp_excl = "results/variants_fail_maf_info_chr{chromosome}.exclude"
        # --exclude {input.snp_excl} \
        # --export vcf-4.2 vcf-dosage=DS 'bgz' \


# GATK liftover failed!

# rule create_sequence_dict:
#     # input: "data/GRCh38_full_analysis_set_plus_decoy_hla.fa.gz"
#     # output: "data/GRCh38_full_analysis_set_plus_decoy_hla.dict"
#     input: "data/chr{chromosome}.fa.gz"
#     output: "data/chr{chromosome}.dict"
#     log: "logs/create_sequence_dict_chr{chromosome}.log"
#     resources: mem=12000, time=60
#     conda: "../envs/picard.yaml"
#     shell: f"{gatk} CreateSequenceDictionary -R {{input}} &> {{log}}"


# rule liftover_grch37_to_grch38:
#     # Change imputed data build coordinates to hg38 to match exome data.
#     # Chain: https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz
#     # Reference sequence: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/chromosomes/
#     # wget https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/md5sum.txt
#     # mv md5sum.txt md5sum_chain.txt
#     # md5sum --check md5sum_chain.txt 2>/dev/null | grep -v 'FAILED open or read'
#     input:
#         vcf = "data/chr{chromosome}.vcf.gz",
#         chain = "data/hg19ToHg38.over.chain.gz",
#         ref = "data/chr{chromosome}.fa.gz",
#         lift_dict = "data/chr{chromosome}.dict" 
#     output:
#         vcf = "data/chr{chromosome}_grch38.vcf.gz",
#         reject = "data/chr{chromosome}_grch38_rejected.vcf"
#     log: "logs/liftover_grch19_to_grch38_chr{chromosome}.log"
#     resources: mem=24000, time=60*10, ntasks=10
#     conda: "../envs/picard.yaml"
#     shell:
#         f"""
#         {gatk} LiftoverVcf \
#         -I {{input.vcf}} \
#         -O {{output.vcf}} \
#         -CHAIN {{input.chain}} \
#         -REJECT {{output.reject}} \
#         -R {{input.ref}} &>> {{log}}
#         """
