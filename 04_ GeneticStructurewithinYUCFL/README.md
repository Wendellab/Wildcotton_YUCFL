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

### FL
#### 
```

```
