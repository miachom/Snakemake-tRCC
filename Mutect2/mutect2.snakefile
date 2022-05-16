configfile: "config/samples.yaml"
configfile: "config/config.yaml" 
 
rule all:
    input:
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz",tumors=config["samples"],chromosomes=config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.tbi",tumors=config["samples"],chromosomes=config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}_f1r2.tar.gz",tumors=config["samples"],chromosomes=config["chromosomes"]),
        expand("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.stats",tumors=config["samples"],chromosomes=config["chromosomes"]),

rule mutect2:
    input:
        tumor_filepath = lambda wildcards: config["samples"][wildcards.tumors]
    output:
        vcf = temp("results/{tumors}/unfiltered_{chromosomes}.vcf.gz"),
        tbi = temp("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.tbi"),
        tar = temp("results/{tumors}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = temp("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.stats")
    params:
        # Edited these to match my config.yaml file
        reference_genome = config["reference_genome"],
        germline_resource = config["germline_resource"],
        gatk = config["gatk_path"],
        panel_of_normals = config["panel_of_normals"],
        normals = lambda wildcards: config["normals"][wildcards.tumors]
    log:
        "logs/mutect2/{tumors}_{chromosomes}_mutect2.txt"
    shell:
        "({params.gatk} Mutect2 \
        -reference {params.reference_genome} \
        -input {input.tumor_filepath} \
        -normal {params.normals} \
        -intervals {wildcards.chromosomes} \
        --germline-resource {params.mutect2_germline_resource} \
        --f1r2-tar-gz {output.tar} \
        --panel-of-normals {params.panel_of_normals} \
        -output {output.vcf}) 2> {log}"
