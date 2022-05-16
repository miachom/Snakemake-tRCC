# Mutect2

## Run
```
/mnt/storage/apps/anaconda3/bin/snakemake -s /home/mi724/Tools/snakefile --cluster-config /home/mi724/Tools//config/cluster_qsub.yaml --cluster "qsub -l h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -pe {cluster.pe} -binding {cluster.binding}" --jobs 30 --rerun-incomplete
```
