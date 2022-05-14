configfile: "config/samples.yaml"
configfile: "config/config.yaml"

rule all:
    input: 
        expand(),
        
rule mutect2:
