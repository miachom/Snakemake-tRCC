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
