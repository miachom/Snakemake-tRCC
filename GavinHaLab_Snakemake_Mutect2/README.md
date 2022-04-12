# How to run the Snakemake workflow for [Mutect2 from Gavin Ha Lab on github](https://github.com/GavinHaLab/mutect2_snakemake/tree/788cb735f5198195e21b392f11034b4506fa39f7)

Clone the repo to your local directory
```
git clone https://github.com/GavinHaLab/mutect2_snakemake.git
```

Samples to be analyzed should be in the config/samples.yaml file. Need to edit this file to include the samples you are working on. Need to be in the same format as the samples.yaml that comes with the repo. Samples need to be in .bam format. The indexed bam files need to have .bam.bai extension for the program to work and need to be in the same directory as your .bam files. 

If you have cram files, use samtools to make them into bam files:
```
samtools view -T <fasta_path> -b -o <bam_file_path> <cram_file_path>
```

What to update in the config/samples.yaml
```
samples: #format is tumor_name: [tumor_filepath, matched_normal_name, matched_normal_filepath]
  tumor_name1:  [/tumor1/filepath.bam, normal_name1, /normal1/filepath.bam]
  tumor_name2:  [/tumor2/filepath.bam, normal_name2, /normal2/filepath.bam]
```

Parameters need to be adjusted in the config/config.yaml file. This includes path to java, gatk, tabix, picard.

What to update in the config/config.yaml
```
reference_genome:
    /path/to/genome.fa
    
mutect2_germline_resource:
    /path/to/somatic-hg38_af-only-gnomad.hg38.vcf.gz
    
known_polymorphic_sites:
    /path/to/somatic-hg38_small_exac_common_3.hg38.vcf.gz

gatk: /path/to/gatk
java: /path/to/java
tabix: /path/to/tabix
picard_jar: /path/to/picard.jar
```

To find Java path use:
```
which java
```
