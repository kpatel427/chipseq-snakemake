# chipseq-snakemake

### To create an enviornment
`conda create -y -n snakemake-env python=3.6 && source activate snakemake-env && conda install -y -c bioconda -c conda-forge snakemake`

### How to run:
Activate enviornment
`source activate snakemake-env`
To run it locally
`snakemake -s Snakefile -j4 all`
To run it on cluster
`snakemake -s Snakefile --latency-wait 20 --cluster "qsub -cwd -V -l h_vmem=20G -l m_mem_free=16G" -j4 all`
