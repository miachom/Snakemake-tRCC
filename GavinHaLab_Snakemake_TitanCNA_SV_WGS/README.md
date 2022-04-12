# How to run the Snakemake workflow for [TitanCNA_SV_WGS from GavinHaLab on git](https://github.com/GavinHaLab/TitanCNA_SV_WGS)

Clone the repo to your local directory
```
git clone https://github.com/GavinHaLab/TitanCNA_SV_WGS.git
```

Need also to clone the following repos:
+ TitanCNA
```
git clone https://github.com/gavinha/TitanCNA.git
```
+ ichorCNA
```
git clone https://github.com/GavinHaLab/ichorCNA.git
```
+ svaba
```
git clone --recursive https://github.com/walaj/svaba
cd svaba
./configure
make
make install
```
+ HMMcopy
```
git clone https://github.com/shahcompbio/hmmcopy_utils.git
cd hmmcopy_utils
cmake .
make
```
This compiles the C++ script hmmcopy_utils/src/bin/readCounter.cpp
The readCounterScript comes from hmmcopy_utils/bin/readCounter (Above script needs to be compiled using cmake . and make)

Samples to be analyzed should be in the config/samples.yaml file. Need to edit this file to include the samples you are working on. Need to be in the same format as the samples.yaml that comes with the repo. Samples need to be in .bam format.
If you have cram files, use samtools to make them into bam files:
```
samtools view -T <fasta_path> -b -o <bam_file_path> <cram_file_path>
```

What to update in the config/samples.yaml
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

What to update in the config/config.yaml
```
## path to tools
samTools:  /path/to/samtools
bcfTools:  /path/to/bcftools
svaba_exe:  /path/to/svaba #/home/ma1111/tools/svaba

readCounterScript:  /path/to/readCounter
ichorCNA_rscript:  /path/to/ichorCNA/scripts/runIchorCNA.R #/home/ma1111/tools/ichorCNA/scripts/runIchorCNA.R

TitanCNA_rscript: /path/to/TitanCNA/scripts/R_scripts/titanCNA.R

TitanCNA_selectSolutionRscript: /path/to/TitanCNA/scripts/R_scripts/selectSolution.R

# include if ichorCNA/TitanCNA R source if files changed but package has not been installed within R
TitanCNA_libdir:  /path/to/TitanCNA
ichorCNA_libdir:  /path/to/ichorCNA 


refFasta:  /path/to/fasta/genome.fa
snpVCF:  /path/to/hapmap_3.3.hg38.vcf.gz
cytobandFile:  /path/to/ichorCNA/inst/extdata/cytoBand_hg38.txt # only need if hg38
centromere:  /path/to/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
sex:  male   # use None if both females and males are in sample set

# must use gc wig file corresponding to same binSize (required)
ichorCNA_gcWig: /path/to/ichorCNA/inst/extdata/gc_hg38_10kb.wig
# must use map wig file corresponding to same binSize (required)
ichorCNA_mapWig:  /path/to/ichorCNA/inst/extdata/map_hg38_10kb.wig
# use bed file if sample has targeted regions, eg. exome data (optional)
ichorCNA_repTimeWig: /path/to/ichorCNA/inst/extdata/RepTiming_hg38_10kb.wig
ichorCNA_centromere:  /path/to/ichorCNA/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt

svaba_dbSNPindelVCF:  /path/to/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.vcf.gz
```

