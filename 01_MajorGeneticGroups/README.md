## Gh(AD1) Gb(AD2) Gm(AD4) Joint Genotyping 

### 1. Getting four lists of Gvcfs
1. AD1 wild: 166 Florida (FL) (141 newly collected + 25 MK);
2. AD1 wild: 158 Yucatan (YUC);
3. AD1 wild: 5 Puerto Rico (PR);
4. AD1 wild: 21 Guadeloupe (GD);
5. AD1 germplasm collections: 10 Cultivar, 10 Landrace1 (LR1), 10 Landrace2 (LR2);
6. AD2 wild: 9 wild from [this manuscript](https://www.biorxiv.org/content/10.1101/2025.06.02.657498v1.abstract)
7. AD4 outgroup: 3 AD4.

### 2.Reads trimming, mapping to reference genome TEX2094, and calling individual gVCFs. 
```
DIR=/work/LAS/jfw-lab/weixuan/03_Ggvcf/AD1_Yucatan/01_readsrename
tDir=/work/LAS/jfw-lab/weixuan/03_Ggvcf/AD1_Yucatan/01_readsrename/trimmedReads
ref=/work/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=20

file1=$(ls -1 $DIR/*_1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(ls -1 $DIR/*_1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p| sed 's/_1[.]fq/_2\.fq/')

name=$(basename $file1 _1.fq.gz)

echo "the first file is " $file1
echo "the second file is " $file2

#trim reads with Trimmomatic:
module load  trimmomatic/0.39-zwxnnrx
cd $Dir
trimmomatic PE -threads $thr $file1 $file2 $tDir/$name.R1.fq.gz $tDir/$name.U1.fq.gz $tDir/$name.R2.fq.gz $tDir/$name.U2.fq.gz ILLUMINACLIP:Adapters.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:75
mv $name.R[12].fq.gz $tDir
module purge

#map reads with bwa:
module load bwa
cd $tDir
bwa mem -M -R "@RG\tID:$name \tSM:$name \tPL:ILLUMINA" -t $thr $ref $tDir/$name.R1.fq.gz $tDir/$name.R2.fq.gz > $name.sam
module purge 

#scoring with sentieon-genomics:
module load sentieon-genomics
sentieon util sort -o $name.sort.bam -t $thr --sam2bam -i $name.sam

#extimate stats with sentieon-genomics:
sentieon driver -t $thr -r $ref -i $name.sort.bam --algo GCBias --summary $name.GC.summary $name.GC.metric --algo MeanQualityByCycle $name.MQ.metric --algo QualDistribution $name.QD.metric --algo InsertSizeMetricAlgo $name.IS.metric --algo AlignmentStat $name.ALN.metric 
sentieon plot metrics -o $name.metric.pdf gc=$name.GC.metric mq=$name.MQ.metric qd=$name.QD.metric isize=$name.IS.metric 

sentieon driver -t $thr -i $name.sort.bam --algo LocusCollector --fun score_info $name.score 
sentieon driver -t $thr -i $name.sort.bam --algo Dedup --rmdup --score_info $name.score --metrics $name.dedup.metric $name.dedup.bam 

#realign with sentieon-genomics:
sentieon driver -t $thr -r $ref -i $name.dedup.bam --algo Realigner $name.realign.bam

#create gvcf files with sentieon-genomics:
#sentieon driver -t $thr -r $ref -i $name.realign.bam  --algo Haplotyper $name.gVCF --emit_mode gvcf 

#coverage caculation:
sentieon driver -t $thr -r $ref -i $name.realign.bam --algo CoverageMetrics coverageoutput/$name

#Base quality score recalibration:
sentieon driver -t $thr -r $ref -i $name.realign.bam --algo QualCal $name.RECAL_DATA.TABLE

#create gvcf files with sentieon-genomics:
sentieon driver -t $thr -r $ref -i $name.realign.bam -q $name.RECAL_DATA.TABLE --algo Haplotyper $name.gVCF --emit_mode gvcf 

#cleanups:
rm $name.sam 
mv $name*.dedup.bam* dedupBam 
mv $name*.metric* metrics
mv $name*.summary* metrics 
mv $name*.realign.bam* realignBam 
mv $name*.sort.bam* sortBam 
mv $name*.score* score
mv $name*.gVCF* gvcf
```

### 3. Using individual gVCFs we joint-called the variant and invariant site from each population/group that we defined in step1, and filtered seven VCFs.

```
seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

Dir=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/
ref=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=30 #NUMBER_THREADS

output1=AD1_YUC_n158

cd $Dir
#mkdir $output1

module load sentieon-genomics/202308.02-e2gz6fb
export SENTIEON_LICENSE=reimu.las.iastate.edu:8990

#joint SNP calling:

cat list_gvcf_YUC_n158.txt | sentieon driver --interval Ah_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Ah_$seq.vcf -
mv $TMPDIR/$output1.Ah_$seq.vcf* $Dir/$output1/

cat list_gvcf_YUC_n158.txt | sentieon driver --interval Dh_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode all $TMPDIR/$output1.Dh_$seq.vcf -
mv $TMPDIR/$output1.Dh_$seq.vcf* $Dir/$output1/

#echo "Filtering VCF file using vcftools and rename samples using bcftools"
cd $Dir/$output1/

ml vcftools bcftools

vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Ah_$seq.variant
vcftools --vcf $output1.Ah_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Ah_$seq.invariant

vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-missing-count 0 --max-alleles 2 --min-meanDP 10 --max-meanDP 100 --mac 2 --recode --recode-INFO-all --out  $output1.Dh_$seq.variant
vcftools --vcf $output1.Dh_$seq.vcf --remove-indels --max-maf 0 --min-meanDP 10 --max-meanDP 100 --recode --out  $output1.Dh_$seq.invariant

module load parallel/20220522-sxcww47

parallel bgzip {} ::: $output1.*h_$seq.*variant.recode.vcf
parallel tabix {} ::: $output1.*h_$seq.*variant.recode.vcf.gz

bcftools concat --allow-overlaps --threads $thr $output1.Ah_$seq.variant.recode.vcf.gz $output1.Ah_$seq.invariant.recode.vcf.gz -Oz -o $output1.Ah_$seq.combined.vcf.gz
bcftools concat --allow-overlaps --threads $thr $output1.Dh_$seq.variant.recode.vcf.gz $output1.Dh_$seq.invariant.recode.vcf.gz -Oz -o $output1.Dh_$seq.combined.vcf.gz

parallel tabix {} ::: $output1.*h_$seq.*combined.vcf.gz
```

### 4. Using individual gVCFs we joint-called the variant and invariant site from each population/group that we defined in step1, and filtered jointed VCFs by chromosome.
```
seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

module purge

ml vcftools bcftools

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_FL_n166/AD1_FL_n166.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_Yuan_n30/AD1_Yuan_n30.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_YUC_n158/AD1_YUC_n158.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_GD_n21/AD1_GD_n21.Ah_$seq.combined.vcf.gz \
 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_PR_n5/AD1_PR_n5.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD2_n9/AD2_n9.Ah_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD4_n3/AD4_n3.Ah_$seq.combined.vcf.gz \
-Oz -o YUCFLAD2AD4_n392.Ah_$seq.combined.vcf.gz

bcftools merge --threads 10 \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_FL_n166/AD1_FL_n166.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_Yuan_n30/AD1_Yuan_n30.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_YUC_n158/AD1_YUC_n158.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_GD_n21/AD1_GD_n21.Dh_$seq.combined.vcf.gz \
 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD1_PR_n5/AD1_PR_n5.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD2_n9/AD2_n9.Dh_$seq.combined.vcf.gz \
/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/AD4_n3/AD4_n3.Dh_$seq.combined.vcf.gz \
-Oz -o YUCFLAD2AD4_n392.Dh_$seq.combined.vcf.gz

bcftools view YUCFLAD2AD4_n392.Ah_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING==0' -q 0.001:minor -Oz -o YUCFLAD2AD4_n392.Ah_$seq.combined.bi.vcf.gz

bcftools view YUCFLAD2AD4_n392.Dh_$seq.combined.vcf.gz --threads 10 \
-m2 -M2 -i 'F_MISSING==0' -q 0.001:minor -Oz -o YUCFLAD2AD4_n392.Dh_$seq.combined.bi.vcf.gz

module load parallel/20220522-sxcww47
parallel tabix {} -f ::: YUCFLAD2AD4_n392.*h_$seq.combined.bi.vcf.gz
```

### 5. Put all chromosomes together back to one single VCF for all samples.
```
module load picard/2.27.4

picard GatherVcfs \
$(for vcf in *.combined.bi.vcf.gz; do echo -I "$vcf"; done) \
-O YUCFLAD2AD4_n392.AhDh.combined.bi.vcf.gz 

ml vcftools bcftools

bcftools reheader -s rename_YUCFLAD2AD4_n392.txt YUCFLAD2AD4_n392.AhDh.combined.bi.vcf.gz -o YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.vcf
bcftools annotate --set-id +"%CHROM:%POS:%REF:%ALT" YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.vcf >  YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
```

### 6. Major genetic group by PCA and NJ tree
```
module load plink/1.90b6.21

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
output=YUCFLAD2AD4_n392

bcftools view $vcf \
-m2 -M2 -i 'F_MISSING=0' -q 0.01:minor \
-o $output.AhDh.combined.bi.rehead.id.maf001.vcf

vcf2=$output.AhDh.combined.bi.rehead.id.maf001.vcf

grep -v "#" $vcf | wc -l 
grep -v "#" $vcf2 | wc -l 

#Filtering LD by Plink 
plink --threads 30 --vcf $vcf2 \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 30 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf2

plink --threads 30 --vcf $vcf2 --allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R
```

### 7. Major genetic group by LEA
To run LEA in parallel we divided LEA into three steps: 

[Step1 bash](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/01_MajorGeneticGroups/step1PCA_lea.sh)  
[Step1 Rscript file convert
](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/01_MajorGeneticGroups/step1PCA_LEA3.R)
```
module load r
Rscript step1PCA_LEA3.R

cut -d ' ' -f 2 *.ped > samplename.txt
```

[Step2 LEA bash in parallel
](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/01_MajorGeneticGroups/step2_lea_parallele.sh)
[Step2 R script 
](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/01_MajorGeneticGroups/step2_run_snmfparallele.R)
```
#SBATCH --output="job.vcf_n90.%J.out"
#SBATCH --job-name="jointGeno_n90"
#SBATCH --array=1-40

module load r

K_values=${SLURM_ARRAY_TASK_ID}
echo "Running LEA with K=$K_values"
mkdir -p "snmf_results/snmf_K${K_values}"

cd snmf_results/snmf_K${K_values}

cp /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/LEA_n392/YUCFLAD2AD4_n392.geno ./
Rscript /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/LEA_n392/step2_run_snmfparallele.R "$K_values"
```

Step 3 in R merge all LEA Ks into one
```
module load r
library(LEA)
for (K in 2:40) { 
project_path <- file.path(paste0("snmf_K", K), "YUCFLAD2AD4_n392.snmfProject") 
combinedProject <- combine.snmfProject("snmf_K1/YUCFLAD2AD4_n392.snmfProject", project_path)}
```

Step4 make the plots with [LEA_plot.R](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/01_MajorGeneticGroups/LEA_plot.R) 

### 8. Major genetic group without outgroup
```

ml bcftools

module load plink/1.90b6.21

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
output=AD1380

bcftools query -l $vcf | grep -vE "AD2|AD4" > subset_$output.txt


bcftools view -S subset_$output.txt $vcf \
-m2 -M2 -i 'F_MISSING=0' -q 0.01:minor \
-o $output.AhDh.combined.maf001.vcf

vcf2=$output.AhDh.combined.maf001.vcf

#Filtering LD by Plink 
plink --threads 30 --vcf $vcf2 \
--indep-pairwise 50 10 0.1 --allow-extra-chr --const-fid \
--out $output

plink --threads 30 --extract $output.prune.in  \
--make-bed --allow-extra-chr --const-fid --out $output \
--recode vcf-iid --vcf $vcf2

plink --threads 30 --vcf $vcf2 \
--allow-extra-chr --const-fid \
--extract $output.prune.in --recode \
--make-bed --pca 20 var-wts --genome --distance square 1-ibs --out $output

paste -d '\t'  $output.mdist.id $output.mdist >  $output.distmatrix

module load r

Rscript PCA_plot.R

```
