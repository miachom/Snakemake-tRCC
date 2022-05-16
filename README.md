# Snakemake
Snakemake is a workflow management system tool. It is python based.
This repository has syntax and explanations of Snakemake file code snippets and structure. It also contains many of the best practice applications required for production pipelines. 

Snakemake workflows are executed in 3 phases:
1. In the initialization phase, the files defining the workflow are parsed and all rules are instantiated. 
2. In the DAG phase, the directed acyclic dependency graph of all jobs is built by filling wildcards and matching input files to output files.
3. In the scheduling phase, the DAG of jobs is executed, with jobs started according to available resources.

## Breakdown of the README.md file
1. [Base structure of the Snakemake file](#Basics-of-Snakemake-Files)
2. [Customizing You Snakemake File](#Customize-Snakemake-File)
    + [expand function](#The-expand-function)
3. [Input Functions](#Input-Functions)
4. [Lambda Expressions and wildcards](#Lambda-Expressions-and-wildcards)
5. [Errors and How-to Solve Them](#Errors-and-How-to-Solve-them)


## Basics of Snakemake Files
Snakemake workflows are defined as rules. Rules decompose the workflow into smaller steps (in WDL these are the tasks).

### rule
The rule contains all the information necessary to "do something" in the script. It contains *input* and *output* file information. It can also run as a *shell* command, using a pythong *script*, etc.

Example:
```
rule rule_name:
  input:
    "path/to/input/file1",
    "path/to/input/file2"
  output:
    "path/to/output/file"
  shell:
    "terminal/command/here"
```

An example of a rule that is using bwa_mem:
```
rule bwa_mem:
  input:
    "path/to/fastq/file1.fq",
    "path/to/reference/genome.fa"
  output:
    "path/to/mapped/reads.bam"
  shell:
    "bwa mem {input} | samtools view -Sb - > {output}"
```

When a workflow is executed, Snakemake tries to generate given target files. Target files can be specified via the command line by executing:
```
snakemake -np path/to/mapped/reads.bam
```
-n is for --dry_run, will show the execution plan instead of performing the steps.
-p flag instructs to print results to ther shell command.

To generate the target files, Snakemake applies the rules given in the Snakemake file in a top-down way.

Apart from filenames, Snakemake also accepts rule names as targets if the requested rule does not have wildcards. Hence, it is possible to write target rules collecting particular subsets of the desired results or all results. Moreover, if no target is given at the command line, Snakemake will define the first rule of the Snakefile as the target. Hence, it is best practice to have a rule **all** at the top of the workflow which has all typically desired target files as input files.    


Can execute the code using:
```
snakemake --cores 1 path/to/mapped/reads.bam
```
Need to specify the number of cores used when executing a workflow. Can allow also for parallelization. 
Snakemake will not try to create path/to/mapped/reads.bam again, because it is already present in the file system. 
Snakemake only re-runs jobs if one of the input files is newer than one of the output files or one of the input files will be updated by another job.

Gerneralize the Snakemake workflow by using named wildcards.

Long Snakemake File
```
SAMPLES = ["A", "B"]


rule all:
    input:
        "plots/quals.svg"


rule bwa_map:
    input:
        "data/genome.fa",
        "data/samples/{sample}.fastq"
    output:
        "mapped_reads/{sample}.bam"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"


rule samtools_sort:
    input:
        "mapped_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam"
    shell:
        "samtools sort -T sorted_reads/{wildcards.sample} "
        "-O bam {input} > {output}"


rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    shell:
        "samtools index {input}"


rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"


rule plot_quals:
    input:
        "calls/all.vcf"
    output:
        "plots/quals.svg"
    script:
        "scripts/plot-quals.py"
```

## Customize Snakemake File

### The expand function

The expand() function allows us to resolve and combine different variables.
```
bam=expand("sorted_reads/{sample}.bam", sample=config["samples"])
```
where
```
samples:
  A: data/samples/A.fastq
  B: data/samples/B.fastq
```
Would provide us with the output:[ A.bam, B.bam]

We can also define the final target files in the *run all* workflow.
```
run all:
    input:
        expand("sorted_reads/{sample}.bam", sample=config["samples"])
```

### Config files
Previously saw adding a SAMPLE list into the above Snakemake file for passing sample information for the workflow. However, need the workflow to be customizable, easily adapted to new data. This is done with the config file mechanism which can be in .json or .yaml and passed with the *configfile* directive. This will be at the top of the Snakemake file and will look as seen below.

Example in workflow:
```
configfile: "config.yaml"
```

When passing the *configfile* directive,use the expand() function and can update the rules in the above sample to:
```
rule bcftools_call:
    input:
        fa="data/genome.fa",
        bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
        bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
    output:
        "calls/all.vcf"
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | "
        "bcftools call -mv - > {output}"
```
This is how to have multiple samples run from the configuration file.

### YAML File Format

Example of adding multiple samples that works for the above example:
```
samples:
  A: data/samples/A.fastq
  B: data/samples/B.fastq
```

## Input Functions

Since we stored the path to our FASTQ files in the cofig file, we can generalize other rules to use these paths.
The expand() function in the list of input files of the rule bcftools_call are executed during the initialization phase. In this phase, we don't know about jobs, wildcard values and rule dependencies. Need to push it to the DAG phase.
We are doing this because we want the paths to the files not just the [A.bam and B.bam file names](#The-expand-function).
Add this function to get all the paths and be able to run BWA on all the files in the configfile

```
def get_bwa_map_input_fastqs(wildcards):
    return config["samples"][wildcards.sample]
    
rule bwa_map:
    input:
        "data/genome.fa",
        get_bwa_map_input_fastqs
    output:
        "mapped_reads/{sample}.bam"
    threads: 8
    shell:
        "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

```

## Lambda Expressions and wildcards
Small functions can be written by using the **lambda** keyword. 
General syntax: lambda *arguments* : *expression*
wildcards can be used in the input, output and params but not in the shell commands inside the rule.

Example:
```
lambda wildcards: config["samples"][wildcards.sample]
```
The input function takes the argument a *wildcards* object, this allows access to the wildcards values. Can return a string or a list of strings that are interpreted as paths to input files.

config file for the example:
```
samples:
    A: data/samples/A.fastq
    B: data/samples/B.fastq
```

## Errors and How-to Solve Them



