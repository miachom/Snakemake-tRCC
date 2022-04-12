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

