## Genomic Diversity Comparison between Wild Cottons using Pixy and Vcftools

### Pixy
#### Step1 filter fixed heterrozygous site by population in each group
```
#SBATCH --job-name="mergevcf"
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

module purge

ml vcftools bcftools

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/04_Pixy_n380/00_vcf_nofixed

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_FL_n166/AD1_FL_n166.Ah_$seq.combined.vcf.gz -Oz -o AD1_FL_n166.Ah_$seq.combined.nofixhe.vcf.gz
bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_FL_n166/AD1_FL_n166.Dh_$seq.combined.vcf.gz -Oz -o AD1_FL_n166.Dh_$seq.combined.nofixhe.vcf.gz

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_GD_n21/AD1_GD_n21.Ah_$seq.combined.vcf.gz -Oz -o AD1_GD_n21.Ah_$seq.combined.nofixhe.vcf.gz
bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_GD_n21/AD1_GD_n21.Dh_$seq.combined.vcf.gz -Oz -o AD1_GD_n21.Dh_$seq.combined.nofixhe.vcf.gz

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_PR_n5/AD1_PR_n5.Ah_$seq.combined.vcf.gz -Oz -o AD1_PR_n5.Ah_$seq.combined.nofixhe.vcf.gz
bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_PR_n5/AD1_PR_n5.Dh_$seq.combined.vcf.gz -Oz -o AD1_PR_n5.Dh_$seq.combined.nofixhe.vcf.gz

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_Yuan_n30/AD1_Yuan_n30.Ah_$seq.combined.vcf.gz -Oz -o AD1_Yuan_n30.Ah_$seq.combined.nofixhe.vcf.gz
bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_Yuan_n30/AD1_Yuan_n30.Dh_$seq.combined.vcf.gz -Oz -o AD1_Yuan_n30.Dh_$seq.combined.nofixhe.vcf.gz

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_YUC_n158/AD1_YUC_n158.Ah_$seq.combined.vcf.gz -Oz -o AD1_YUC_n158.Ah_$seq.combined.nofixhe.vcf.gz
bcftools view --exclude "F_PASS(GT='het')=1" --threads 5 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_YUC_n158/AD1_YUC_n158.Dh_$seq.combined.vcf.gz -Oz -o AD1_YUC_n158.Dh_$seq.combined.nofixhe.vcf.gz


module load parallel/20220522-sxcww47
parallel tabix -f {} ::: *.*h_$seq.combined.nofixhe.vcf.gz

output=AD1_n380
mkdir -p outputn380_combined

bcftools merge --threads 5 \
AD1_FL_n166.Ah_$seq.combined.nofixhe.vcf.gz \
AD1_GD_n21.Ah_$seq.combined.nofixhe.vcf.gz \
AD1_PR_n5.Ah_$seq.combined.nofixhe.vcf.gz \
AD1_Yuan_n30.Ah_$seq.combined.nofixhe.vcf.gz \
AD1_YUC_n158.Ah_$seq.combined.nofixhe.vcf.gz \
-Oz -o outputn380_combined/$output.Ah_$seq.combined.vcf.gz

bcftools merge --threads 5 \
AD1_FL_n166.Dh_$seq.combined.nofixhe.vcf.gz \
AD1_GD_n21.Dh_$seq.combined.nofixhe.vcf.gz \
AD1_PR_n5.Dh_$seq.combined.nofixhe.vcf.gz \
AD1_Yuan_n30.Dh_$seq.combined.nofixhe.vcf.gz \
AD1_YUC_n158.Dh_$seq.combined.nofixhe.vcf.gz \
-Oz -o outputn380_combined/$output.Dh_$seq.combined.vcf.gz

cd outputn380_combined 
parallel tabix -f {} ::: $output.*h_$seq.combined.vcf.gz
```

#### Step2 prepare input population txt files for each population
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz

ml bcftools
bcftools query -l /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz \
| awk -F'_' '{grp=$2; if($0~/^AD1_YUC_RiCa/||$0~/^AD1_YUC_RiCh/) grp="YUCE"; else if($0~/^AD1_YUC/) grp="YUCW"; print $0 "\t" grp}' > pixy_populationlist_8group.txt

cp pixy_populationlist_8group.txt pixy_populationlist_domwild.txt
awk 'BEGIN{OFS="\t"} 
$2=="FL" || $2=="GD" || $2=="PR" || $2=="YUCE" || $2=="YUCW" {$2="wild"} 
{print}' pixy_populationlist_domwild.txt > tmp && mv tmp pixy_populationlist_domwild.txt
```

#### Step4 running pixy Global, i.e., Cultivar, Landrace1, Landrace2 and Wild (all wild cotton populations as one) 
```
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})

seqA="Ah_"$seq
seqD="Dh_"$seq

thr=30 #NUMBER_THREADS
mkdir -p 02_output_n380_domwild

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz 
outputfolder=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/02_output_n380_domwild
output1=AD1380_domwild

module purge

module load py-numpy/1.26.3-py310-gntgk2n
module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/pixy_env

echo $seqA
echo $seqD

pixy --stats pi fst dxy --bypass_invariant_check --chromosomes $seqA --vcf $vcf --populations pixy_populationlist_domwild.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqA
pixy --stats pi fst dxy --bypass_invariant_check --chromosomes $seqD --vcf $vcf --populations pixy_populationlist_domwild.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqD

micromamba deactivate
```

#### Step5 running pixy Regionally, seperate each wild cotton population  
```
#SBATCH --array=1-13

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})

seqA="Ah_"$seq
seqD="Dh_"$seq

thr=30 #NUMBER_THREADS
mkdir -p 01_output_n380_8groups

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz 
outputfolder=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/08_Pixy_n380/01_output_n380_8groups
output1=AD1380_8groups

module purge

module load py-numpy/1.26.3-py310-gntgk2n
module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/pixy_env

echo $seqA
echo $seqD

pixy --stats pi fst dxy watterson_theta tajima_d --bypass_invariant_check --chromosomes $seqA --vcf $vcf --populations pixy_populationlist_8group.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqA
pixy --stats pi fst dxy watterson_theta tajima_d --bypass_invariant_check --chromosomes $seqD --vcf $vcf --populations pixy_populationlist_8group.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqD

micromamba deactivate
```
#### Step6 tabluating output of pixy use [caculated.sh](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/03_GenomicDiversityComparisonbetweenWildCottons/caculated.sh) 
