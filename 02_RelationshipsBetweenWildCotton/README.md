## Genetic relationship inference between wild cotton populations

### 1. Estimating Gene Flow among Groups using Treemix and Dsuite
```
ml vcftools bcftools
module load plink/1.90b6.21

output=YUCFLAD2AD4_n392; vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/01_Plink_n392/YUCFLAD2AD4_n392.vcf

vcftools --vcf $vcf --plink-tped --out $output

awk '{split($1,a,"_"); $1=a[1]"_"a[2]; if($2~/^AD1_YUC_RiCa/||$2~/^AD1_YUC_RiCh/){$1="AD1_YUCE"} else if($2~/^AD1_YUC/){$1="AD1_YUCW"} print $0}' $output.tfam > tmp && mv tmp $output.tfam

awk '{print $1"\t"$2"\t"$1}' $output.tfam > sample.pop.cov

/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/plink --threads 10 --tfile $output --freq --allow-extra-chr --within sample.pop.cov
gzip plink.frq.strat
python /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/plink2treemix.py plink.frq.strat.gz sample.treemix.in.gz

mkdir -p Treemixoutput
module purge
module load gsl/2.7.1-uuykddp python/3.10.10-zwlkg4l boost/1.81.0-zwxu2hi py-numpy/1.26.3-py310-gntgk2n

for m in {0..8}; do for i in {1..5}; do /lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/src/treemix -se -bootstrap -i sample.treemix.in.gz -root AD4_mus -o Treemixoutput/TreeMix.${m}.${i} -m ${m} -k 1000; done; done

mkdir -p Dsuite

/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/treemix-1.13/src/treemix -i sample.treemix.in.gz -o sample.ML.tree -root AD4_mus -k 1000
zcat sample.ML.tree.treeout.gz | sed 's/AD4_mus/Outgroup/g' > Dsuite/sp.tre
cut -f 2,3 sample.pop.cov | awk '{gsub("AD4_mus","Outgroup",$2); print}' > Dsuite/sets.txt
sed -i 's/ /\t/g' Dsuite/sets.txt

cd Dsuite
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/Build/Dsuite Dtrios -c -n $output -t sp.tre $vcf sets.txt
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/Build/Dsuite Fbranch sp.tre sets_${output}_tree.txt > species_sets_${output}_Fbranch.txt
/lustre/hdd/LAS/jfw-lab/weixuan/07_PRGD_popgene/03_AD1AD2AD4/Dsuite/utils/dtools.py species_sets_${output}_Fbranch.txt sp.tre --dpi 300 --tree-label-size 18 --ladderize

sed -i "s/AD1_/Gh_/g;s/AD2_Wild/Gb/g;s/YUCW/YUC-W/g;s/YUCE/YUC-E/g" sp.tre
sed -i "s/AD1_/Gh_/g;s/AD2_Wild/Gb/g;s/YUCW/YUC-W/g;s/YUCE/YUC-E/g" species_sets_${output}_Fbranch.txt

#R part for TreeMix model selection

library(OptM)
linear <- optM("./Treemixoutput")
plot_optM(linear, plot = FALSE, pdf = "Treemodelselection_file.pdf")
```

### 2. Introgression Estimation using [Flare](https://github.com/browning-lab/flare)
#### step1 subset VCF to exclude AD4, remove fixed heterozygous sites, and phase VCF via Beagle
```
# subset and filter
module load vcftools bcftools
chr=$(printf %02d ${SLURM_ARRAY_TASK_ID})
output=YUCFLAD2AD4_n389.bi

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
bcftools query -l $vcf | grep -v "AD4_" > keep_samples.txt
vcftools --gzvcf $vcf --keep keep_samples.txt --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $output
bcftools view --include "F_PASS(GT='het')=1" $output.recode.vcf -o $output.nofixed.vcf

# phase by chromosome
#SBATCH --array=1-13

newvcf=$output.nofixed.vcf
module load openjdk/21.0.3_9-vngib7s

java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$newvcf out=$output.Ah_$chr.phased chrom=Ah_$chr
java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$newvcf out=$output.Dh_$chr.phased chrom=Dh_$chr

# joint everything together
ml vcftools bcftools
for vcf in *phased*.gz; do tabix "$vcf"; done
bcftools concat -Oz -o YUCFLAD2AD4_n389.biphased.vcf.gz $(for vcf in *phased*.gz; do echo "$vcf"; done)
bcftools index YUCFLAD2AD4_n389.biphased.vcf.gz
```

#### step2 prepare FLARE input
```
#wget https://faculty.washington.edu/browning/flare.jar

module load openjdk/21.0.3_9-vngib7s

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/05_Beagle_n389_nooutgroup/YUCFLAD2AD4_n389.biphased.vcf.gz

ml bcftools

## we first extract AD2 and YUC as reference, to find out the relationships between other wild cotton populations vs these two.

bcftools query -l $vcf | grep -E 'AD2|YUC' > set1_refgrouplist_AD2YUC.txt
bcftools view -S set1_refgrouplist_AD2YUC.txt $vcf -Oz -o set1_refgrouplist_AD2YUC.vcf.gz
awk -F'_' '{print $0 "\t" $1 "_"$2 }' set1_refgrouplist_AD2YUC.txt > set1_refpanel_AD2YUC.txt

bcftools query -l $vcf | grep -E 'FL|GD|Ph' > set1_testgrouplist_FLGDPh.txt
bcftools view -S set1_testgrouplist_FLGDPh.txt $vcf -Oz -o set1_testgrouplist_FLGDPh.vcf.gz

java -Xmx100g -jar flare.jar ref=set1_refgrouplist_AD2YUC.vcf.gz \
	ref-panel=set1_refpanel_AD2YUC.txt \
	gt=set1_testgrouplist_FLGDPh.vcf.gz \
	map=Ghirsutum_map_sort_clean.txt \
	out=set1_refAD2YUC_testFLGDPh


## we second test all domesticated cottons vs all AD1 wild cotton populations and AD2

bcftools query -l $vcf | grep -E 'AD2|FL|YUC|GD|Ph' > set2_refgrouplist_AD2YUCFLPRGDPh.txt
bcftools view -S set2_refgrouplist_AD2YUCFLPRGDPh.txt $vcf -Oz -o set2_refgrouplist_AD2YUCFLPRGDPh.vcf.gz
awk -F'_' '{print $0 "\t" $1 "_"$2 }' set2_refgrouplist_AD2YUCFLPRGDPh.txt > set2_refpanel_AD2YUCFLPRGDPh.txt

bcftools query -l $vcf | grep -E 'Cultivar|LR1|LR2' > set2_testgrouplist_CulLR1LR2.txt
bcftools view -S set2_testgrouplist_CulLR1LR2.txt $vcf -Oz -o set2_testgrouplist_CulLR1LR2.vcf.gz


java -Xmx100g -jar flare.jar ref=set2_refgrouplist_AD2YUCFLPRGDPh.vcf.gz \
	ref-panel=set2_refpanel_AD2YUCFLPRGDPh.txt \
	gt=set2_testgrouplist_CulLR1LR2.vcf.gz \
	map=Ghirsutum_map_sort_clean.txt \
	out=set2_refAD2YUCFLPRGDPh_testCulLR1LR2
```

### 2. Kmer-based Genetic Relationship Inference between Major Genetic Groups
```
```

### 3. Plastome Variation 
```
```
