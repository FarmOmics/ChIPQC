#!/bin/bash -l

module load snakemake/5.32.2
source activate snakemake-5.32.2

mkdir -p Logs
snakemake -j 1000 --cluster-config /group/zhougrp3/zhoulab3/liqi/fish_testD/ChIPQC/cluster.yaml \
	-s Snakefile \
	--configfile config.yaml \
	--cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} -c {cluster.cpus} --mem={cluster.mem} -J {cluster.name} -o {cluster.output} -e {cluster.output}" \
	--nolock --latency-wait 480 -p -k --use-conda $@

