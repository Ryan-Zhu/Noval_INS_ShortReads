#!/usr/bin/env bash
#SBATCH --export ALL
#SBATCH --partition hotel
#SBATCH --account htl189
#SBATCH --qos hotel
#SBATCH --job-name gencDNA
#SBATCH --nodes 1
#SBATCH --ntasks 5
#SBATCH --cpus-per-task 4
#SBATCH --mem-per-cpu 4G   
#SBATCH --time 24:00:00
#SBATCH --output /tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/logs_0907/slurm-%j.out-%N-2100


module load cpu/0.17.3

launch_script=/tscc/nfs/home/yzhu/projects/gencDNA/scripts/Find_insertions-4-25.py

samples=/tscc/nfs/home/yzhu/projects/gencDNA/2024june_TSCC_bwaPaths_V2.txt
target_genes=/tscc/nfs/home/yzhu/projects/gencDNA/analysis_0907/target_genes_0824.tsv

outdir=/tscc/lustre/ddn/scratch/yzhu/gencDNA_0924/analysis_1111/2100

mkdir -p $outdir

head -n 2100 $samples  | tail -n 700 | parallel python3 $launch_script {} -o $outdir -f $target_genes -m 20 -c 30
