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
# step1 subset and filter
module load vcftools bcftools
chr=$(printf %02d ${SLURM_ARRAY_TASK_ID})
output=YUCFLAD2AD4_n389.bi

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
bcftools query -l $vcf | grep -v "AD4_" > keep_samples.txt
vcftools --gzvcf $vcf --keep keep_samples.txt --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out $output
bcftools view --include "F_PASS(GT='het')=1" $output.recode.vcf -o $output.nofixed.vcf

# step2 phase by chromosome
#SBATCH --array=1-13

newvcf=$output.nofixed.vcf
module load openjdk/21.0.3_9-vngib7s

java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$newvcf out=$output.Ah_$chr.phased chrom=Ah_$chr
java -Xmx300g -jar beagle.27Feb25.75f.jar gt=$newvcf out=$output.Dh_$chr.phased chrom=Dh_$chr

# step3 joint everything together
ml vcftools bcftools
for vcf in *phased*.gz; do tabix "$vcf"; done
bcftools concat -Oz -o YUCFLAD2AD4_n389.biphased.vcf.gz $(for vcf in *phased*.gz; do echo "$vcf"; done)
bcftools index YUCFLAD2AD4_n389.biphased.vcf.gz
```

#### step2 prepare FLARE input. Genetic map from [Zhang et al 2019](https://link.springer.com/article/10.1186/s12864-019-6214-z) 
```
## FLARE assumes that genetic distance increases monotonically along the chromosome. Let's first fix the genetic map file, by removing markers where genetic distance decreases compared to the previous marker in the file. 
awk '{chr=$1; dist=$3;
if(chr!=last_chr){last_dist=-1}
if(dist>=last_dist){print; last_dist=dist} last_chr=chr}' \
Ghirsutum_map_sort.txt > Ghirsutum_map_sort_clean.txt


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
#### step1 array clean the reads data using [kraken2](https://benlangmead.github.io/aws-indexes/k2) and [mirabait](https://sourceforge.net/projects/mira-assembler/files/MIRA/development/)
```
#SBATCH --array=253-325

DIR=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/00_readsfilter
thr=20

file1=$(sed -n ${SLURM_ARRAY_TASK_ID}p list_AD1AD2AD4AD5n223.txt)
file2=$(sed -n ${SLURM_ARRAY_TASK_ID}p list_AD1AD2AD4AD5n223.txt | sed 's/R1[.]fq/R2\.fq/')
name=$(basename $file1 .R1.fq.gz)

echo "the first file is " $file1
echo "the second file is " $file2
echo "sample name is" $name

#cat /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/00_TEX2094/ncbiAD1mtDNA.fasta /lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle/00_TEX2094/plastome_output/TX2094-USDA/embplant_pt.K105.complete.graph1.1.path_sequence.fasta > chloro_mito.fasta
#wget "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz"
#mkdir -p kraken2db
#tar -xzf k2_standard_08_GB_20251015.tar.gz -C kraken2db

module load kraken2/2.1.2-2kmlv3r

kraken2 $file1 $file2 \
    --paired \
	--db kraken2db/ \
    --classified-out ${name}_contaminated#.fastq \
    --unclassified-out ${name}_clean#.fastq \
    --use-names \
    --gzip-compressed \
    --threads $thr \
    --output ${name}_output.txt \
    --report ${name}_report.txt

mkdir -p 00_cleanread
mkdir -p 01_cleanreadreport
mv ${name}_clean_*.fastq 00_cleanread/
mv ${name}_report.txt 01_cleanreadreport
rm ${name}_output.txt
rm ${name}_contaminated_*.fastq

mkdir -p 02_chloro_mito

#wget "https://sourceforge.net/projects/mira-assembler/files/MIRA/development/mira_4.9.6_linux-gnu_x86_64_static.tar.bz2/download"
#tar -xvjf download

chloro_mito=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/00_readsfilter/chloro_mito.fasta

/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/00_readsfilter/mira_V5rc1_linux-gnu_x86_64_static/bin/mirabait \
-I \
-o 02_chloro_mito \
-N $name \
-b $chloro_mito \
-p 00_cleanread/${name}_clean_1.fastq 00_cleanread/${name}_clean_2.fastq \
-t $thr

mkdir -p 03_nuclear
module load pigz/2.8-vt5ps6w

pigz -p $thr -c "${name}_miss_${name}_clean_1.fastq" > "03_nuclear/${name}_clean_R1.fq.gz"
pigz -p $thr -c "${name}_miss_${name}_clean_2.fastq" > "03_nuclear/${name}_clean_R2.fq.gz"

rm ${name}_miss_${name}_clean_*.fastq
```
#### step2 set up KmerCity configuration file 
```
Data:
    FastQ_table: "/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/fastq_table_n158.txt"

Outputs:
    Output_directory: "/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n158_results"
    Output_prefix: "n158"

Parameters:
    kat: "--mer_len 50 -H 1000000000"
    cut_count_threshold: "30"
    kmer_db_chunk_size: "1000000"

Envmodules:
    seqtk: "seqtk/1.3-r106"
    kat: "kat/2.4.2-conda"
    jellyfish: "jellyfish/2.3.0"

Conda: "workflow/envs/KmerCity_env.yaml"

Container:
    seqtk: "https://depot.galaxyproject.org/singularity/fusioncatcher-seqtk%3A1.2--h5bf99c6_2"
    kat: "https://depot.galaxyproject.org/singularity/kat%3A2.4.2--py36h873903e_2"
    jellyfish: "https://depot.galaxyproject.org/singularity/kmer-jellyfish%3A2.3.0--h7d875b9_2"
```

#### step3 also `/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/kmercity/slurm/config.yaml` needs to update to system
```
#Rename runtime → walltime in default-resources and set-resources.
#Update your cluster template to use {resources.walltime}
#Leave SLURM --time=72:00:00 as-is.
#Valid partition (nova)  --partition=nova

cluster:
  mkdir -p logs/{rule} &&
  sbatch
    --partition={resources.partition}
    --mem={resources.mem_mb}
    --time={resources.walltime}
    --job-name={rule}
    --output=logs/{rule}/{rule}_%j.log
    --error=logs/{rule}/{rule}_%j.log

default-resources:
  - mem_mb=4096
  - walltime="72:00:00"
  - partition=nova

set-resources:
  - subset_fastq:mem_mb=15000
  - count_kmers:mem_mb=40000
  - build_kmer_db:mem_mb=20000

set-threads:
  - count_kmers=6

restart-times: 3
max-jobs-per-second: 1
max-status-checks-per-second: 1
local-cores: 1
latency-wait: 60
jobs: 100
keep-going: True
rerun-incomplete: True
```
#### step3 prepare input `fastq_table_n158.txt` table use reads information in `/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/00_readsfilter/03_nuclear/*.gz`
```
#SBATCH --output="job.KmerCity.%J.out"
#SBATCH --job-name="KmerCity"

(
echo -e "Sample\tFastQ\tNbReads\tSeed"
ls /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/00_readsfilter/03_nuclear/*.gz | sort | while read f; do
    sample=$(basename "$f" | sed 's/_clean.*//')
    # If the sample starts with "YUC", keep only the part before the first underscore
    if [[ $sample == YUC* ]]; then
        sample=$(echo "$sample" | sed 's/_.*//')
    fi
    seed=$(od -vAn -N2 -tu2 < /dev/urandom | tr -d ' ')
    echo -e "${sample}\t${f}\t6000000\t${seed}"
done
) | head -n317 > fastq_table_n158.txt


module load snakemake/7.22.0-py310-y4hvhdl
module load singularity/1.1.9-py310-wsbt4ge

cd /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/kmercity

srun --time=72:00:00 \
     --partition=nova \
     --job-name="KmerCity_controller" \
     --output="job.KmerCity_controller.%J.out" \
	 --cpus-per-task=10 \
     --mem=150G \
     --cpu-bind=none \
     snakemake --configfile config/AD1AD2config.yaml \
	 --use-singularity \
	 --profile slurm \
	 --jobs 40 \
	 --rerun-incomplete


module load perl/5.40.0-l2sxfqz
module load zlib/1.2.13-wcef3x6
module load perl-list-moreutils/0.430-6om3ljo
module load perl-parallel-forkmanager/2.02-vh7wtb2
module load perl-gd/2.53-u7fgagv


perl /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/kmercity/workflow/scripts/KCounts_2_intersections.pl \
    -database /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/KMER_database/n137_50mer_cut30x_counts.tab.gz \
    -outdir /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/ \
    -prefix n158 
```

#### step4 clean Kmercity output
```
################################################################################
################################################################################

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/FastQ
ls -1 -d * > ../INTERSEC-n158/samplelist_n158.txt

cp /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/rename_YUCFLAD2AD4_n392.txt /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/INTERSEC-n158
sed 's/\\ //g' rename_YUCFLAD2AD4_n392.txt > samplerename.txt

awk 'NR==FNR{keep[$1];next} $1 in keep' samplelist_n158.txt samplerename.txt > matched_n158.txt

cut -d" " -f2 matched_n158.txt | awk -F"_" '{group=$1"_"$2; if(!(group in c)){colors[0]="red";colors[1]="blue";colors[2]="green";colors[3]="orange";colors[4]="purple";colors[5]="cyan";colors[6]="magenta";colors[7]="gray";colors[8]="brown";colors[9]="pink"; c[group]=colors[++n%10]} print $0 "\t" c[group]}' | sort > graph_accession.tab

awk 'NR==FNR{map[$1]=$2;next}{for(i=1;i<=NF;i++){if($i in map)$i=map[$i]}}1' samplerename.txt n158.intersections > n158.intersections.rename.txt

################################################################################
################################################################################

perl /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/kmercity/workflow/scripts/Draw_KGraph.pl \
-in /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/INTERSEC-n158/n158.intersections.rename.txt \
-list /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/INTERSEC-n158/graph_accession.tab \
-outprefix /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/01_KmerCity/n137_results/INTERSEC-n158/n158_graph \
-classes 10,20,30
```

#### step 5 count kmer results in R and group by population using [sortingKmer.R](https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/02_RelationshipsBetweenWildCotton/sortingKmer.R)

### 3. Plastome Variation 
#### step1 assemble plastomes
```
#SBATCH --output=joblog/job.makePR.%A_%a.out 
#SBATCH --job-name="getorganelle"
#SBATCH --array=1-158

module load getorganelle/1.7.7.0-py310-u45ybv3

DIR=/lustre/hdd/LAS/jfw-lab/weixuan/03_Ggvcf/AD1_Yucatan/01_readsrename/trimmedReads
outputDIR=/lustre/hdd/LAS/jfw-lab/weixuan/09_getorgenelle
thr=10

file=$(ls -1 $DIR/*.R1.fq.gz | sed -n ${SLURM_ARRAY_TASK_ID}p)
name=$(basename $file .R1.fq.gz)

get_organelle_from_reads.py -t $thr -1 $DIR/"${name}".R1.fq.gz -2 $DIR/"${name}".R2.fq.gz -o $outputDIR/plastome_output/"${name}" -R 15 -k 21,45,65,85,105 -F embplant_pt
```

#### step2 align sequences and build a tree
```
module load mafft
module load iqtree2
file=plastome_subset_n392_sort_A1outgroup_noAD567outgroup.fasta

# MAFFT realign
realign="${file%.fasta}_realign.fasta"
mafft --thread 80 "$file" > "$realign"

# trimAl
trimmed="${realign%.fasta}_trim.fasta"
/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/trimal/source/trimal -in "$realign" -out "$trimmed" -nogaps

# IQ-TREE2
iqtree2 -s "$trimmed" --prefix "$trimmed" -T AUTO -m MFP -B 1000

micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/newick_env
nw_ed plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.treefile 'i & b<=80' o > plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.bs80.treefile

grep -E '>A1' plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta | sed 's/[>]//g' > A1outgrouplist.txt
module load py-ete3
python reroot_trees.py plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.bs80.treefile A1outgrouplist.txt > plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.bs80.reroot.treefile
```
Whole plastome sequences: https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/02_RelationshipsBetweenWildCotton/plastome_subset_n392_sort_A1outgroup_noAD567outgroup.zip
Aligned whole plastome sequences: https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/02_RelationshipsBetweenWildCotton/plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign.zip
Raw plastome tree: https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/02_RelationshipsBetweenWildCotton/plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.treefile
Rerooted and node collapsed plastome tree: https://github.com/Wendellab/Wildcotton_YUCFL/blob/main/02_RelationshipsBetweenWildCotton/plastome_subset_n392_sort_A1outgroup_noAD567outgroup_realign_trim.fasta.bs80.reroot.treefile
