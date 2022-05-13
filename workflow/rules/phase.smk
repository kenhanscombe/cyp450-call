

# From PharmGKB (https://www.pharmgkb.org/gene/PA124/overview) 

# CYP2C19
# GRCh37	chr10 : 96522463 - 96612671
# GRCh38	chr10 : 94762681 - 94853205

# CYP2D6
# GRCh37	chr22 : 42522501 - 42526883
# GRCh38	chr22 : : 42125531 - 42130881

# How much of the chromosome needs to be phased?
# Need to phase up and downstream, e.g., rs12248560 used to call CYP2C19
# is at GRCh38 94761900

# Use UKB exome block coordinates:
# chr10, block 24, coordinates 93795159 - 95869829

# chr22, block 15, coordinates 40410785 - 41603020
# chr22, block 16, coordinates 41603021 - 42642420
# chr22, block 17, coordinates 42642421 - 44186823

# For input to SHAPEIT4 --region
# Note. indicates how chromosome is recorded in the file, e.g., 22 or chr22
region = {"10": "10:93795159-95869829", "22": "22:40410785-44186823"}
wes_region = {"10": "chr10:93795159-95869829", "22": "chr22:40410785-44186823"}


rule bcftools_sort:
    """Sort the imputed data"""
    input: "data/chr{chromosome}.vcf.gz"
    output: "data/chr{chromosome}_sorted.vcf.gz"
    log: "logs/bcftools_sort_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    threads: 4
    resources: mem="38G", time="1-0", ntasks=4
    shell:
        """
        bcftools sort {input} --output-type z --output {output} 2> {log}
        """


rule tabix_index_imputed:
    """Index the imputed data VCF files."""
    input: "data/chr{chromosome}_sorted.vcf.gz"
    output: "data/chr{chromosome}_sorted.vcf.gz.tbi"
    log: "logs/tabix_index_imputed_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    shell: "tabix {input} 2> {log}"


rule shapeit_phase_imputed:
    """SHAPEIT4 phase imputed data.
    
    Genetic maps https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
    """
    input:
        vcf = "data/chr{chromosome}_sorted.vcf.gz",
        index = "data/chr{chromosome}_sorted.vcf.gz.tbi",
        maps = "maps/chr{chromosome}.b38.gmap.gz"
    output: "data/chr{chromosome}_phased.vcf.gz"
    params: region = lambda w: region[str(w)]
    log: wrk_dir / "logs/shapeit_phase_imputed_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    envmodules:
        "utilities/use.dev",
        "apps/shapeit4/4.2.0-gcc7.3.0"
    threads: 4
    resources: mem="64G", time="0-10", ntasks=4
    shell:
        f"""
        shapeit4.2 \
        --input {{input.vcf}} \
        --map {{input.maps}} \
        --region {{params.region}} \
        --thread 10 \
        --output {{output}} &> {{log}}
        """


rule tabix_index_phased_imputed:
    """Index the imputed data VCF files."""
    input: "data/chr{chromosome}_phased.vcf.gz"
    output: "data/chr{chromosome}_phased.vcf.gz.tbi"
    log: "logs/tabix_index_phased_imputed_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    shell: "tabix {input} 2> {log}"


# Whole exome sequencing data

# rule bcftools_sort_wes:
#     """Sort the WES data"""
#     input: "data/chr{chromosome}_wes.vcf.gz"
#     output: "data/chr{chromosome}_wes_sorted.vcf.gz"
#     log: "logs/bcftools_sort_chr{chromosome}_wes.log"
#     conda: "../envs/ukbqc.yaml"
#     resources: mem="164G", time="1-0"
#     shell:
#         """
#         bcftools sort {input} -Oz -o {output} 2> {log}
#         """

# ALT: PLINK format
rule bcftools_sort_plink_wes:
    """Sort the PLINK format WES data"""
    input: "data/chr{chromosome}_plink_wes.vcf.gz"
    output: "data/chr{chromosome}_wes_sorted.vcf.gz"
    log: "logs/bcftools_sort_plink_wes_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem="164G", time="1-0"
    shell:
        """
        bcftools sort {input} --output-type z --output {output} 2> {log}
        """


rule tabix_index_wes:
    """Index the WES data VCF files."""
    input: "data/chr{chromosome}_wes_sorted.vcf.gz"
    output: "data/chr{chromosome}_wes_sorted.vcf.gz.tbi"
    log: "logs/tabix_index_wes_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    shell: "tabix {input} 2> {log}"


rule shapeit_phase_wes:
    """SHAPEIT4 phase WES data.
    
    Genetic maps https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
    """
    input:
        vcf = "data/chr{chromosome}_wes_sorted.vcf.gz",
        tbi = "data/chr{chromosome}_wes_sorted.vcf.gz.tbi",
        maps = "maps/chr{chromosome}.b38.gmap.gz"
    output: "data/chr{chromosome}_wes_phased.vcf.gz"
    params: region = lambda w: region[str(w)]
    log: wrk_dir / "logs/shapeit_phase_wes_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    envmodules:
        "utilities/use.dev",
        "apps/shapeit4/4.2.0-gcc7.3.0"
    threads: 10
    resources: mem='164G', time='1-0', ntasks=10
    shell:
        f"""
        shapeit4.2 \
        --input {{input.vcf}} \
        --map {{input.maps}} \
        --region {{params.region}} \
        --thread 10 \
        --output {{output}} &> {{log}}
        """


rule tabix_index_phased_wes:
    """Index the WES data VCF files."""
    input: "data/chr{chromosome}_wes_phased.vcf.gz"
    output: "data/chr{chromosome}_wes_phased.vcf.gz.tbi"
    log: "logs/tabix_index_phased_wes_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    shell: "tabix {input} 2> {log}"



