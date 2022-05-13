# rule get_wes_regions:
#     """Get hg38 exome regions from PGxPOP"""
#     output:
#         wes_regions = "data/exome_target_regions.bed",
#         chr10 = "data/exome_chr10_regions.bed",
#         chr22 = "data/exome_chr22_regions.bed"
#     log: "logs/get_wes_regions.log"
#     conda: "../envs/ukbqc.yaml"
#     shell:
#         """
#         wget -nc -O {output.wes_regions} \
#         https://raw.githubusercontent.com/PharmGKB/PGxPOP/master/definition/exome_target_regions.bed \
#         2> {log}
        
#         awk '($4 == "CYP2C19")' {output.wes_regions} > {output.chr10} 2>> {log}
#         awk '($4 == "CYP2D6")' {output.wes_regions} > {output.chr22} 2>> {log}
#         """


# rule bcftools_trim_wes_regions:
#     """Remove exome target regions from imputed data"""
#     input:
#         imp = "data/chr{chromosome}_sorted.vcf.gz",
#         regions = "data/exome_chr{chromosome}_regions.bed"
#     output:
#         trim = "data/chr{chromosome}_trim.vcf.gz"
#     log: "logs/bcftools_trim_wes_regions_chr{chromosome}.log"
#     conda: "../envs/ukbqc.yaml"
#     envmodules:
#         "utilities/use.dev",
#         "apps/bedtools2/2.29.0"
#     resources: mem='64G', time='0-1'
#     shell:
#         """
#         bedtools intersect -a {input.imp} -b {input.regions} -lof > {output.trim} 2> {log}
#         """


rule bcftools_merge_imputed_wes:
    """Merge imputed and WES data"""
    input:
        imp = "data/chr{chromosome}_sorted.vcf.gz",
        wes = "data/chr{chromosome}_wes_sorted.vcf.gz"
    output:
        merged = "data/chr{chromosome}_merged.vcf.gz"
    params:
        temp = "data/temp_chr{chromosome}"
    log: "logs/bcftools_merge_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    threads: 4
    resources: mem='164G', time='0-2'
    shell:
        """
        mkdir -p {params.temp}
        
        bcftools isec -p {params.temp} {input.imp} {input.wes} --threads 4 -Oz 2> {log}
        bcftools query -l {params.temp}/0000.vcf.gz | awk 'gsub ("_.*", "")' > {params.temp}/0000.ids
        bcftools reheader {params.temp}/0000.vcf.gz -s {params.temp}/0000.ids --threads 4 -o {params.temp}/0000_id_update.vcf.gz
        bcftools concat {params.temp}/0000_id_update.vcf.gz {input.wes} --threads 4 -Oz -o {output.merged} 2>> {log}

        rm -rf {params.temp}
        """


rule bcftools_sort_merged:
    """Sort the merged data"""
    input: "data/chr{chromosome}_merged.vcf.gz"
    output: "data/chr{chromosome}_merged_sorted.vcf.gz"
    log: "logs/bcftools_sort_chr{chromosome}_merged.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem='164G', time='1-0'
    shell:
        """
        bcftools sort {input} --output-type z --output {output} 2> {log}
        """
        # --max-mem 1000 


rule tabix_index_merged:
    """Index the merged data VCF files."""
    input: "data/chr{chromosome}_merged_sorted.vcf.gz"
    output: "data/chr{chromosome}_merged_sorted.vcf.gz.tbi"
    log: "logs/tabix_index_merged_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    resources: mem='256G', time='1-0'
    shell: "tabix {input} 2> {log}"


rule shapeit_phase_merged:
    """SHAPEIT4 phase merged data.
    
    Genetic maps https://github.com/odelaneau/shapeit4/raw/master/maps/genetic_maps.b38.tar.gz
    """
    input:
        vcf = "data/chr{chromosome}_merged_sorted.vcf.gz",
        tbi = "data/chr{chromosome}_merged_sorted.vcf.gz.tbi",
        maps = "maps/chr{chromosome}.b38.gmap.gz"
    output: "data/chr{chromosome}_merged_phased.vcf.gz"
    params: region = lambda w: region[str(w)]
    log: wrk_dir / "logs/shapeit_phase_merged_chr{chromosome}.log"
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


rule tabix_index_phased_merged:
    """Index the merged data VCF files."""
    input: "data/chr{chromosome}_merged_phased.vcf.gz"
    output: "data/chr{chromosome}_merged_phased.vcf.gz.tbi"
    log: "logs/tabix_index_phased_merged_chr{chromosome}.log"
    conda: "../envs/ukbqc.yaml"
    shell: "tabix {input} 2> {log}"