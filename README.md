# Snakemake
Snakemake is a workflow management system tool. It is python based.
This repository has syntax and explanations of Snakemake file code snippets and structure. It also contains many of the best practice applications required for production pipelines. 

## Breakdown of the README.md file
1. [Base structure of the Snakemake file](#Basics)


## Basics of Snakemake files
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

When a workflow is executed, Snakemake trie to generate given target files. Target files can be specified via the command line by executing:
```
snakemake -np path/to/mapped/reads.bam
```
-n is for --dry_run, will show the execution plan instead of performing the steps.
-p flag instructs to print results to ther shell command.

To generate the target files, Snakemake applies the rules given in the Snakemake file in a top-down way.

Can execute the code using:
```
snakemake --cores 1 path/to/mapped/reads.bam
```
Need to specify the number of cores used when executing a workflow. Can allow also for parallelization. 
Snakemake will not try to create path/to/mapped/reads.bam again, because it is already present in the file system. 
Snakemake only re-runs jobs if one of the input files is newer than one of the output files or one of the input files will be updated by another job.

Gerneralize the Snakemake workflow by using named wildcards.
