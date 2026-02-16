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
```
```

### 2. Kmer-based Genetic Relationship Inference between Major Genetic Groups
```
```

### 3. Plastome Variation 
```
```
