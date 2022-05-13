gene = {"10": "CYP2C19", "22": "CYP2D6"}
cyp_snps = {"10": ["rs4244285", "rs12248560"],
            "22": ["rs5758550", "rs16947"]}


rule pgxpop_call_imputed:
    """Call PGx alleles in imputed data."""
    input: 
        vcf = "data/chr{chromosome}_phased.vcf.gz",
        tbi = "data/chr{chromosome}_phased.vcf.gz.tbi"
    output: "data/chr{chromosome}.pgxpopout"
    params:
        pgx_gene = lambda wildcards: gene[str(wildcards.chromosome)],
        out_dir = directory("out/chr{chromosome}")
    log: wrk_dir / "logs/pgxpop_call_imputed_chr{chromosome}.log"
    conda: "../envs/pgxpop.yaml"
    # threads: 10
    # resources: mem='164G', time='1-0', ntasks=10
    resources: mem='32G', time='0-3'
    shell:
        f"""
        {pgxpop} \
        --vcf {{input.vcf}} \
        -g {{params.pgx_gene}} \
        --phased \
        --build grch38 \
        -o {{output}} &> {{log}}
        """

rule pgxpop_call_wes:
    """Call PGx alleles in WES data."""
    input:
        vcf = "data/chr{chromosome}_wes_phased.vcf.gz",
        tbi = "data/chr{chromosome}_wes_phased.vcf.gz.tbi"
    output: "data/chr{chromosome}_wes.pgxpopout"
    params:
        pgx_gene = lambda wildcards: gene[str(wildcards.chromosome)]
    log: wrk_dir / "logs/pgxpop_call_wes_chr{chromosome}.log"
    conda: "../envs/pgxpop.yaml"
    resources: mem='164G', time='1-0'
    shell:
        f"""
        {pgxpop} \
        --vcf {{input.vcf}} \
        -g {{params.pgx_gene}} \
        --phased \
        --build grch38 \
        -o {{output}} &> {{log}}
        """


rule pgxpop_call_merged:
    """Call PGx alleles in merged data."""
    input:
        vcf = "data/chr{chromosome}_merged_phased.vcf.gz",
        tbi = "data/chr{chromosome}_merged_phased.vcf.gz.tbi"
    output: "data/chr{chromosome}_merged.pgxpopout"
    params:
        pgx_gene = lambda wildcards: gene[str(wildcards.chromosome)]
    log: wrk_dir / "logs/pgxpop_call_merged_chr{chromosome}.log"
    conda: "../envs/pgxpop.yaml"
    resources: mem='164G', time='1-0'
    shell:
        f"""
        {pgxpop} \
        --vcf {{input.vcf}} \
        -g {{params.pgx_gene}} \
        --phased \
        --build grch38 \
        -o {{output}} &> {{log}}
        """


# Calling from imputed SNPs for comparison
# - Chiara's algorithm

# rule bgenix_get_snps:
#     """Retrieve SNPs used to call CYP2C19/ CYP2D6"""
#     input: imp_dir / "ukb_imp_chr{chromosome}_v3_MAF0_INFO7.bgen",
#     output: 
#         snp_file = "data/cyp_snp_chr{chromsome}.txt",
#         bgen = "data/cyp_snp_chr{chromosome}.bgen"
#     params:
#         snps = lambda wildcards: " ".join(cyp_snps[str(wildcards.chromosome)])
#     log: "logs/bgenix_get_snps_chr{chromosome}.log"
#     shell:
#         f"""
#         printf {{params.snps}} > {{output.snp_file}}
#         {bgenix} -g {{input}} -incl-rsids {{output.snp_file}} > {{output.bgen}} 2> {{log}}
#         """


# rule plink2_cyp_snp_stats:
#     """Calculate SNP statistics for CYP2C19/ CYP2D6 snps"""
#     input: "data/cyp_snp_chr{chromosome}.bgen"
#     output:
#     log: "plink2_cyp_snp_stats_chr{chromsome}.log"
#     conda: "../envs/pgxpop.yaml"
#     envmodules: "apps/plink2/2.0.0a2"
#     shell:
#         """
#         plink2 \
#         --bgen data/cyp_snp_chr22.bgen 'ref-first' \
#         --sample /scratch/datasets/ukbiobank/ukb56514/imputed/wukb56514_imp_chr1_v3_s487296.sample \
#         --freq \
#         --missing variant-only \
#         --hardy 'midp' \
#         --geno-counts
#         """

# Then I coded the phenotypes based on Chiaras cyp2c19 paper using the following code:
 
# awk '{ if ($7 ~/CC/ && $8 ~/AA/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ poor_met/' > cyp2c19_poor_met
# awk '{ if ($7 ~/CC/ && $8 ~/AG/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ intermediate_met/' > cyp2c19_intermediate_met
# awk '{ if ($7 ~/TC/ && $8 ~/AG/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ intermediate_plus_met/' > cyp2c19_intermediate_plus_met
# awk '{ if ($7 ~/CC/ && $8 ~/GG/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ extensive_met/' > cyp2c19_extensive_met
# awk '{ if ($7 ~/TC/ && $8 ~/GG/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ extensive_plus_met/' > cyp2c19_extensive_plus_met
# awk '{ if ($7 ~/TT/ && $8 ~/GG/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ ultrarapid_met/' > cyp2c19_ultrarapid_met
# awk '{ if ($7 ~/00/ || $8 ~/00/) print }' ukb56514_chr10_cyp2c19_snps_genotypes.ped | sed 's/$/ missing_met/' > cyp2c19_missing_met