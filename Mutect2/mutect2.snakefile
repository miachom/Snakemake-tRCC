configfile: "config/samples.yaml"
configfile: "config/config.yaml" 

rule all:
    input:
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),

rule mutect2:
    input:
        tumor_filepath = list(config["samples"].values())
        
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
