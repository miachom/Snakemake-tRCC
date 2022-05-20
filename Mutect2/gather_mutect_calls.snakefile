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

rule learn_read_orientation_model:
    input:
        f1r2 = config["f1r2"]
    output:
        protected("results/{tumors}/read_orientation_model.tar.gz")
    params:
        gatk = config["gatk_path"]
    log:
        "logs/learn_read_orientation_model/{tumors}_learn_read_orientation_model.txt"
    shell:
        "({params.gatk} LearnReadOrientationModel \
        -I {input.f1r2} \
        -O {output}) 2> {log}"

rule gather_mutect_calls:
    input:
        vcf_files = config["vcf_files"]
    output:
        temp("results/{base_file_name}/gathered_unfiltered.vcf.gz")
    params:
        java = config["java"],
        picard_jar = config["picard_jar"]
    log:
        "logs/gather_mutect_calls/{base_file_name}_gather_mutect_calls.txt"
    shell:
        "({params.java} -jar {params.picard_jar} GatherVcfs \
        I={input.vcf_files} \
        O={output}) 2> {log}"
 

