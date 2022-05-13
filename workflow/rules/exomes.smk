wes_block = {"10": "24", "22": "16"}


rule wget_exome_block_list:
    """
    Get exome data block list.
    """
    output: "data/pvcf_blocks.txt"
    log: "logs/wget_exome_block_list.log"
    conda: "../envs/ukbqc.yaml"
    shell:
        """
        wget  \
        -O data/pvcf_blocks.txt \
        -nd  biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/pvcf_blocks.txt 2> {log}
        """


# rule gfetch_exome_blocks:
#     """
#     Get exome blocks containing PGx gene
#     """
#     input:
#         exome_key = "/scratch/datasets/ukbiobank/ukb56514/raw/k56514x38594.key"
#     output:
#     log: "logs/gfetch_exome_blocks_{locus}.log"
#     params:
        
#         loc = "{locus}"
#     conda: "../envs/ukbqc.yaml"
#     shell:
#         """
#         {gfetch} 23156 -c10 -b24 -a{input.exome_key} 2> {log}
#         """


rule bcftools_exome_qc:
    """
    Sample filter: Use the same sample as used for the imputed data (i.e.,
    remove samples failing UKB SQC, KING cutoff related)
    """
    input:
        vcf = "data/chr{chromosome}.vcf.gz",
        wes_vcf = lambda w: f"data/ukb23156_c{w}_b{wes_block[str(w)]}_v1.vcf.gz",
    output: "data/chr{chromosome}_wes.vcf.gz"
    params: "data/temp_chr{chromosome}.keep"
    log: "logs/bcftools_exome_qc_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem='128G', time='0-5'
    shell:
        """
        bcftools query -l {input.vcf} | sed '{{s/_.*$//}}' > {params}
        bcftools view -Oz -S {params} --force-samples {input.wes_vcf} -o {output} &> {log}

        rm {params}
        """
        # bcftools view -Oz -S {params} --force-samples --threads 10 -m2 -M2 -v snps {input.wes_vcf} -o {output} &> {log}



# PLINK format exome data

rule make_exome_fam:
    """
    Use UKB-supplied bridging file to make 56514 exome fam file
    """
    input:
        fam = "/scratch/datasets/ukbiobank/exomes/ukb23155_c20_b0_v1_s200632.fam",
        bridge = "/scratch/datasets/ukbiobank/ukb56514/raw/bridge_56514_18177.csv"
    output: "data/exome.fam"
    log: "logs/make_exome_fam.log"
    shell:
        """
        python workflow/scripts/make_wes_fam.py -b {input.bridge} -f {input.fam} -o {output} &> {log}
        """


rule write_variants_to_keep:
    """Range to keep for PLINK to VCF conversion"""
    output: "data/variant_plink_wes.keep",
    log: "logs/write_variants_to_keep.log",
    conda: "../envs/ukbqc.yaml"
    resources: mem='3G', time='0-0:10'
    shell:
        """
        echo "10 93795159 95869829 CYP2C19" > {output}
        echo "22 40410785 44186823 CYP2D6" >> {output}
        """

rule plink2_wes_to_vcf:
    """
    Convert PLINK format exome data to VCF

    Writes a temporary sample keep file from the filtered imputed data,
    and a temporary variant range extract file in PLINK 'set range'
    format to include (pVCF) exome block including gene and one block
    either side.
    """
    input:
        imputed_vcf = "data/chr{chromosome}.vcf.gz",
        bed = "/scratch/datasets/ukbiobank/exomes/ukb23155_c{chromosome}_b0_v1.bed",
        bim = "/scratch/datasets/ukbiobank/exomes/UKBexomeOQFE_chr{chromosome}.bim",
        fam = "data/exome.fam",
        extract = "data/variant_plink_wes.keep"
    output:
        vcf = "data/chr{chromosome}_plink_wes.vcf.gz"
    params:
        keep = temp("data/sample_chr{chromosome}.keep"),
        out = lambda w, output: os.path.splitext(os.path.splitext(output.vcf)[0])[0]
    log: "logs/plink2_wes_to_vcf_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    envmodules: plink2
    resources: mem='32G', time='0-1'
    shell:
        """
        bcftools query -l {input.imputed_vcf} | \
        sed '{{s/_.*$//}}' | \
        awk '{{print $1, $1}}' > {params.keep}

        plink2 \
        --bed {input.bed} \
        --bim {input.bim} \
        --fam {input.fam} \
        --keep {params.keep} \
        --extract 'range' {input.extract} \
        --hwe 1e-50 \
        --export vcf-4.2 vcf-dosage=DS 'bgz' \
        --out {params.out} 2> {log}
        """
        # --hwe 1e-25 \


# bcftools_collapse_multiple_alleles
# bcftools norm