#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=22
#SBATCH --mem=100G 
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.vcf_n90.%J.out"
#SBATCH --job-name="jointGeno_n90"
#SBATCH --array=1-20

module load r

K_values=${SLURM_ARRAY_TASK_ID}
echo "Running LEA with K=$K_values"
mkdir -p "snmf_results/snmf_K${K_values}"

cd snmf_results/snmf_K${K_values}

cp /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/02_AD1FL_n166/02_FLn196_LEA/LEA/FLn196.geno ./
Rscript /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/02_AD1FL_n166/02_FLn196_LEA/LEA/step2_run_snmfparallele.R "$K_values"

for (K in 2:20) { 
project_path <- file.path(paste0("snmf_K", K), "FLn196.snmfProject") 
combinedProject <- combine.snmfProject("snmf_K1/FLn196.snmfProject", project_path)}

#cut -d ' ' -f 2 *.ped > samplename.txt
