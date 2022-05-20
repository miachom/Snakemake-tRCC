import os

configfile: "config/samples.yaml"
configfile: "config/config.yaml" 

rule all:
    input:
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"])

rule mutect2:
    input:
        tumor_filepath = config["samples"]
        
    output:
        vcf = temp("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz"),
        tbi = temp("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi"),
        tar = temp("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = temp("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats")
    params:
        # Edited these to match my config.yaml file
        reference_genome = config["reference_genome"],
        germline_resource = config["germline_resource"],
        gatk = config["gatk_path"],
        panel_of_normals = config["panel_of_normals"],
        normals = config["normals"]
        
    log:
        "logs/mutect2/{base_file_name}_{chromosomes}_mutect2.txt"
    shell:
        "({params.gatk} Mutect2 \
        -reference {params.reference_genome} \
        -input {input.tumor_filepath} \
        -normal {params.normals} \
        -intervals {wildcards.chromosomes} \
        --germline-resource {params.germline_resource} \
        --f1r2-tar-gz {output.tar} \
        --panel-of-normals {params.panel_of_normals} \
        -output {output.vcf}) 2> {log}"


rule gather_mutect_calls:
    input:
        chr1_calls = "results/{base_file_name}/unfiltered_chr1.vcf.gz",
        chr2_calls = "results/{base_file_name}/unfiltered_chr2.vcf.gz",
        chr3_calls = "results/{base_file_name}/unfiltered_chr3.vcf.gz",
        chr4_calls = "results/{base_file_name}/unfiltered_chr4.vcf.gz",
        chr5_calls = "results/{base_file_name}/unfiltered_chr5.vcf.gz",
        chr6_calls = "results/{base_file_name}/unfiltered_chr6.vcf.gz",
        chr7_calls = "results/{base_file_name}/unfiltered_chr7.vcf.gz",
        chr8_calls = "results/{base_file_name}/unfiltered_chr8.vcf.gz",
        chr9_calls = "results/{base_file_name}/unfiltered_chr9.vcf.gz",
        chr10_calls = "results/{base_file_name}/unfiltered_chr10.vcf.gz",
        chr11_calls = "results/{base_file_name}/unfiltered_chr11.vcf.gz",
        chr12_calls = "results/{base_file_name}/unfiltered_chr12.vcf.gz",
        chr13_calls = "results/{base_file_name}/unfiltered_chr13.vcf.gz",
        chr14_calls = "results/{base_file_name}/unfiltered_chr14.vcf.gz",
        chr15_calls = "results/{base_file_name}/unfiltered_chr15.vcf.gz",
        chr16_calls = "results/{base_file_name}/unfiltered_chr16.vcf.gz",
        chr17_calls = "results/{base_file_name}/unfiltered_chr17.vcf.gz",
        chr18_calls = "results/{base_file_name}/unfiltered_chr18.vcf.gz",
        chr19_calls = "results/{base_file_name}/unfiltered_chr19.vcf.gz",
        chr20_calls = "results/{base_file_name}/unfiltered_chr20.vcf.gz",
        chr21_calls = "results/{base_file_name}/unfiltered_chr21.vcf.gz",
        chr22_calls = "results/{base_file_name}/unfiltered_chr22.vcf.gz",
        chrX_calls = "results/{base_file_name}/unfiltered_chrX.vcf.gz",
        chrY_calls = "results/{base_file_name}/unfiltered_chrY.vcf.gz"
    output:
        temp("results/{base_file_name}/gathered_unfiltered.vcf.gz")
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/gather_mutect_calls/{base_file_name}_gather_mutect_calls.txt"
    shell:
        "({params.java} -jar {params.picard_jar} GatherVcfs \
        I={input.chr1_calls} \
        I={input.chr2_calls} \
        I={input.chr3_calls} \
        I={input.chr4_calls} \
        I={input.chr5_calls} \
        I={input.chr6_calls} \
        I={input.chr7_calls} \
        I={input.chr8_calls} \
        I={input.chr9_calls} \
        I={input.chr10_calls} \
        I={input.chr11_calls} \
        I={input.chr12_calls} \
        I={input.chr13_calls} \
        I={input.chr14_calls} \
        I={input.chr15_calls} \
        I={input.chr16_calls} \
        I={input.chr17_calls} \
        I={input.chr18_calls} \
        I={input.chr19_calls} \
        I={input.chr20_calls} \
        I={input.chr21_calls} \
        I={input.chr22_calls} \
        I={input.chrX_calls} \
        I={input.chrY_calls} \
        O={output}) 2> {log}"

