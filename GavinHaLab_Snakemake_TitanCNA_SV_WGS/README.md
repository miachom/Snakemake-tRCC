# How to run the Snakemake workflow for [TitanCNA_SV_WGS from GavinHaLab on git](https://github.com/GavinHaLab/TitanCNA_SV_WGS)

Clone the repo to your local directory
```
git clone https://github.com/GavinHaLab/TitanCNA_SV_WGS.git
```

Samples to be analyzed should be in the config/samples.yaml file. Need to edit this file to include the samples you are working on. Need to be in the same format as the samples.yaml that comes with the repo.

samples.yaml
```
samples: #format is sample_name: /path/to/sample.bam for both tumors and normals
  tumor_sample_name1:  /path/to/tumor/sample1.bam
  tumor_sample_name2:  /path/to/tumor/sample2.bam
  normal_sample_name1:  /path/to/normal/sample1.bam
  normal_sample_name2:  /path/to/normal/sample2.bam

pairings: #format is tumor: matched_normal
  tumor_sample_name1: normal_sample_name1
  tumor_sample_name2: normal_sample_name2
```

Parameters need to be adjusted in the config/config.yaml file. These include the paths to the TitanCNA and ichorCNA repositories, paths to your NGS tools (samtools, bcftools, svaba) and readCounterScript.

Need to load these modules in your terminal before running
```
ml snakemake/5.19.2-foss-2019b-Python-3.7.4
ml R/3.6.2-foss-2019b-fh1
ml Python/3.7.4-foss-2019b-fh1
ml BCFtools/1.9-GCC-8.3.0
ml Pysam/0.15.4-GCC-8.3.0-Python-3.7.4
ml PyYAML/5.1.2-GCCcore-8.3.0-Python-3.7.4
```
