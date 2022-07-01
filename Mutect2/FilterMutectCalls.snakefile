import os

configfile: "config/config.yaml"
configfile: "config/samples.yaml"

rule all:
    input:
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.tbi",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}_f1r2.tar.gz",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/unfiltered_{chromosomes}.vcf.gz.stats",base_file_name=config["base_file_name"],chromosomes=config["chromosomes"]),
        expand("results/{base_file_name}/gathered_unfiltered.vcf.gz",base_file_name=config["base_file_name"]),
        expand("results/{base_file_name}/mutect_merged.stats", base_file_name = config["base_file_name"]),
        expand("results/{base_file_name}/read_orientation_model.tar.gz", base_file_name = config["base_file_name"]),
        expand("results/{base_file_name}/tumor_{base_file_name}_pileup.tab",base_file_name=config["base_file_name"]),
        expand("results/{base_file_name}/normal_{base_file_name}_pileup.tab", base_file_name = config["base_file_name"]),
        expand("results/{base_file_name}/{base_file_name}_tum_segments.tab", base_file_name = config["base_file_name"]),
        expand("results/{base_file_name}/{base_file_name}_matched_contamination.tab", base_file_name = config["base_file_name"]),
        expand("results/{base_file_name}/{base_file_name}_f1r2_filtered_somatic_vcf.gz", base_file_name = config["base_file_name"])

rule Mutect2:
     input:
	tumor_filepath = lambda wildcards: config["samples"][wildcards.tumor],
	normal_filepath = lambda wildcards: config["samples"][config["pairings"][wildcards.tumor]]
     output:
        vcf = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz"),
        tbi = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.tbi"),
        tar = protected("results/{tumors}/unfiltered_{chromosomes}_f1r2.tar.gz"),
        stats = protected("results/{tumors}/unfiltered_{chromosomes}.vcf.gz.stats")
     params:
        reference_genome = config["reference_genome"],
        mutect2_germline_resource = config["mutect2_germline_resource"],
        gatk = config["gatk"],
        panel_of_normals = config["panel_of_normals"],
        normals = lambda wildcards: config["samples"][wildcards.tumors][1]
     log:
        "logs/mutect2/{tumors}_{chromosomes}_mutect2.txt"
     shell:
        "({params.gatk} Mutect2 \
        -reference {params.reference_genome} \
        -input {input.tumor_filepath} \
        -input {input.normal_filepath} \
        -normal {params.normals} \
        -intervals {wildcards.chromosomes} \
        --germline-resource {params.mutect2_germline_resource} \
        --f1r2-tar-gz {output.tar} \
        --panel-of-normals {params.panel_of_normals} \
        -output {output.vcf}) 2> {log}"
     

rule MergeMutectStats:
     output:
        protected("results/{tumors}/mutect_merged.stats")
     params:
        gatk = config["gatk_path"]
     logs:
        "logs/MergeMutectStats/{tumors}_merge_mutect_stats.txt"
     shell:
        "
	all_stat_inputs=`for chromosome in {chromosomes}; do
        printf -- "-stats results/{tumors}/unfiltered_${chromosome}.vcf.gz.stats "; done`

	({params.gatk} MergeMutectStats \
        $all_stat_inputs \
        -O {output}) 2> {log}"

rule LearnReadOrientationModel:
      output:
        protected("results/{tumors}/read_orientation_model.tar.gz")
      params:
        gatk = config["gatk_path"]
      log:
        "logs/LearnReadOrientationModel/{tumors}_learn_read_orientation_model.txt"
      shell:
        "
	all_f1r2_inputs=`for chromosome in {chromosomes}; do
        printf -- "-stats results/{tumors}/unfiltered_${chromosome}_f1r2.tar.gz "; done`
	
	({params.gatk} LearnReadOrientationModel \
	$all_f1r2_inputs \
	-O {output}) 2> {log}"

rule GetPileupSummariesTumor:
     input:
        tumor_filepath = config["samples"]
     output:
        table = protected("results/{tumors}/tumor_{tumors}_pileup.tab")
     params:
        gatk = config["gatk_path"],
        small_exac_common = config["small_exac_common"],
        reference_genome = config["reference_genome"]
     log:
        "logs/GetPileupSummaries/{tumors}_GetPileupSummaries_tumor.txt"
     shell:
         "({params.gatk} GetPileupSummaries \
         -reference {params.reference_genome} \
         -I /mnt/scratch/DTRCC10/DTRCC10_Rkidney.bam  \
         -V {params.small_exac_common} \
         -L {params.small_exac_common} \
         -O {output.table}) 2> {log}"

rule GetPileupSummariesNormal:
     input:
        tumor_filepath = config["samples"]
     output:
        table = protected("results/{tumors}/normal_{tumors}_pileup.tab")
     params:
        gatk = config["gatk_path"],
        small_exac_common = config["small_exac_common"],
        reference_genome = config["reference_genome"]
     log:
        "logs/GetPileupSummaries/{tumors}_GetPileupSummaries_tumor.txt"
     shell:
         "({params.gatk} GetPileupSummaries \
         -reference {params.reference_genome} \
         -I /mnt/scratch/DTRCC10/DTRCC10_Blood.bam  \
         -V {params.small_exac_common} \
         -L {params.small_exac_common} \
         -O {output.table}) 2> {log}"


rule CalculateContamination:
      input:
         tumor = expand("results/{base_file_name}/tumor_{base_file_name}_pileup.tab",base_file_name=config["base_file_name"]),
         matched = expand("results/{base_file_name}/normal_{base_file_name}_pileup.tab", base_file_name = config["base_file_name"])
      output:
         table = protected("results/{tumors}/{tumors}_matched_contamination.tab"),
         tum_seg = protected("results/{tumors}/{tumors}_tum_segments.tab")
      params:
        gatk = config["gatk_path"]
      log:
          "logs/CalculateContamination/{tumors}.log"
      shell:
          "({params.gatk} CalculateContamination \
          -I {input.tumor} \
          -matched {input.matched} \
          -tumor-segmentation {output.tum_seg} \
          -O {output.table}) 2> {log}"

rule FilterMutectCalls:
      input:
           tum_seg = expand("results/{base_file_name}/{base_file_name}_tum_segments.tab",base_file_name=config["base_file_name"]),
           matched_contamination =  expand("results/{base_file_name}/{base_file_name}_matched_contamination.tab",base_file_name=config["base_file_name"])
      output:
           filtered_f1r2 = protected("results/{tumors}/{tumors}_f1r2_filtered_somatic_vcf.gz")
      log:
           "logs/FilterMutectCalls/{tumors}.log"
      params:
           gatk = config["gatk_path"],
           reference_genome = config["reference_genome"],
           tum_seg = "results/{tumors}/{tumors}_tum_segments.tab"
      shell:
           "({params.gatk} FilterMutectCalls \
           -R {params.reference_genome} \
           -V /home/mi724/Tools/Snakemake/Mutect2/results/DTRCC_10/gathered_unfiltered.vcf.gz \
           -tumor-segmentation {params.tum_seg} \
           --contamination-table {input.matched_contamination} \
           --ob-priors /home/mi724/Tools/Snakemake/Mutect2/results/DTRCC_10/read_orientation_model.tar.gz \
           -O {output.filtered_f1r2}) 2> {log}"
