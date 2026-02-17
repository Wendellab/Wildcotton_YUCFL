## Genetic structure within YUC and FL

### YUC 
#### subset samples for n158 (Yucatan only) and n188 (Yucatan + 30 domesticated cotton) 
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/09_HeInbreedRoh/AD1_n380.AhDh.combined.rehead.recode.vcf

ml vcftools bcftools

bcftools query -l $vcf | grep -E 'YUC|Cultivar|LR' > YUCn188.txt
bcftools query -l $vcf | grep -E 'YUC' > YUCn158.txt

vcftools \
  --vcf $vcf \
  --keep YUCn188.txt \
  --maf 0.01 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out YUCn188.maf001

vcftools \
  --vcf $vcf \
  --keep YUCn158.txt \
  --maf 0.01 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out YUCn158.maf001

bcftools annotate -Oz \
  -o YUCn188.maf001.id.vcf.gz \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  YUCn188.maf001.recode.vcf
  
bcftools annotate -Oz \
  -o YUCn158.maf001.id.vcf.gz \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  YUCn158.maf001.recode.vcf
```
#### PCA and PI_HAT relatedness 
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/01_AD1YUC_n158/YUCn158.maf001.id.vcf.gz 
plink=/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/plink
output=YUCn158


#Filtering LD by Plink 
$plink --threads 5 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

$plink --threads 5 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

$plink --threads 5 --vcf $vcf \
--allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --genome --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```
#### NJ tree and LEA 
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/01_AD1YUC_n158/YUCn188.maf001.id.vcf.gz
plink=/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/plink
output=YUCn188


#Filtering LD by Plink 
$plink --threads 5 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

$plink --threads 5 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

$plink --threads 5 --vcf $vcf \
--allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --genome --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```
#
### FL
#### subset samples for N166 (Florida only) and N196 (Florida + 30 domesticated cotton) 
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/09_HeInbreedRoh/AD1_n380.AhDh.combined.rehead.recode.vcf
ml vcftools bcftools

bcftools query -l $vcf | grep -E 'FL|Cultivar|LR' > FLn196.txt
bcftools query -l $vcf | grep -E 'FL' > FLn166.txt

vcftools \
  --vcf $vcf \
  --keep FLn196.txt \
  --maf 0.01 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out FLn196.maf001

vcftools \
  --vcf $vcf \
  --keep FLn166.txt \
  --maf 0.01 \
  --min-alleles 2 \
  --max-alleles 2 \
  --recode \
  --recode-INFO-all \
  --out FLn166.maf001

ml vcftools bcftools

bcftools annotate -Oz \
  -o FLn196.maf001.id.vcf.gz \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  FLn196.maf001.recode.vcf
 
bcftools annotate -Oz \
  -o FLn166.maf001.id.vcf.gz \
  --set-id '%CHROM:%POS:%REF:%ALT' \
  FLn166.maf001.recode.vcf 
```
#### PCA and PI_HAT relatedness 
```
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/02_AD1FL_n166/FLn166.maf001.id.vcf.gz 
plink=/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/plink
output=FLn166


#Filtering LD by Plink 
$plink --threads 5 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

$plink --threads 5 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

$plink --threads 5 --vcf $vcf \
--allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --genome --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```
#### NJ tree and LEA 
```

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/02_AD1FL_n166/FLn196.maf001.id.vcf.gz
plink=/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/plink
output=FLn196


#Filtering LD by Plink 
$plink --threads 5 --vcf $vcf \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

$plink --threads 5 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf

$plink --threads 5 --vcf $vcf \
--allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --genome --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```
#
### Pixy (dxy tabulate by site)
#### tabluating dxy by site (total 25 sites)
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=100G 
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="joboutput/job.pixy_YUC158.%J.out"
#SBATCH --job-name="pixy_YUC158"
#SBATCH --array=1-13

mkdir -p 

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})

seqA="Ah_"$seq
seqD="Dh_"$seq

thr=30 #NUMBER_THREADS
vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/04_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz 
outputfolder=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/04_Pixy_n380/01_output_n350_25Pop
output1=AD1350_25Pop

ml bcftools
bcftools query -l /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/04_Pixy_n380/00_vcf_nofixed/outputn380_combined/AD1_n380.AhDh.combined.rehead.vcf.gz  | awk -F'_' '{print $0 "\t" $2}' > pixy_populationlist_7group.txt
awk -F'_' '{print $0 "\t" $2"_"$3}' pixy_populationlist_7group.txt |  cut -f1,3 |  sed 's/\tGD_[^\t]*/\tGD/' |  grep -v -E 'LR1|LR2|Cultivar' > pixy_populationlist_25Poplist.txt

module purge

module load py-numpy/1.26.3-py310-gntgk2n
module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/pixy_env

echo $seqA
echo $seqD

pixy --stats pi dxy --bypass_invariant_check --chromosomes $seqA --vcf $vcf --populations pixy_populationlist_25Poplist.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqA
pixy --stats pi dxy --bypass_invariant_check --chromosomes $seqD --vcf $vcf --populations pixy_populationlist_25Poplist.txt --window_size 10000 --n_cores $thr --output_folder $outputfolder --output_prefix $output1.$seqD

micromamba deactivate
```

