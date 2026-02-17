## SIFT4G

### Build database

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G 
#SBATCH --time=3-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.sift4g.%J.out"
#SBATCH --job-name="sift4g"

module load singularity/1.3.6-py311-nvfjdsj
module load perl-bioperl/1.7.6-niqnczv

singularity pull docker://juliahoglund/sift4g
wget "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/05_SIFT4G/SIFT4G_Create_Genomic_DB

singularity exec \
    -B $PWD \
    sift4g.sif \
    perl make-SIFT-db-all.pl \
    -config test_files/gossypium_hirsutum.txt
```
### Configuration file 
```
GENETIC_CODE_TABLE=1
GENETIC_CODE_TABLENAME=Standard

PARENT_DIR=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/05_SIFT4G/SIFT4G_Create_Genomic_DB/test_files/gossypium_hirsutum
ORG=gossypium_hirsutum
ORG_VERSION=TEX2094v2
DBSNP_VCF_FILE=

#Running SIFT4G, this path works for the Dockerfile
SIFT4G_PATH=/sift4g/bin/sift4g

#PROTEIN_DB needs to be uncompressed
PROTEIN_DB=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/05_SIFT4G/SIFT4G_Create_Genomic_DB/test_files/gossypium_hirsutum/uniprot_sprot.fasta

# Sub-directories, don't need to change
GENE_DOWNLOAD_DEST=gene-annotation-src
CHR_DOWNLOAD_DEST=chr-src
LOGFILE=Log.txt
ZLOGFILE=Log2.txt
FASTA_DIR=fasta
SUBST_DIR=subst
ALIGN_DIR=SIFT_alignments
SIFT_SCORE_DIR=SIFT_predictions
SINGLE_REC_BY_CHR_DIR=singleRecords
SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores
DBSNP_DIR=dbSNP

# Doesn't need to change
FASTA_LOG=fasta.log
INVALID_LOG=invalid.log
PEPTIDE_LOG=peptide.log
ENS_PATTERN=ENS
SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord
```

### Subset n381 from `YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf` and run SIFT4G
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=300G 
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.sift4g.%J.out"
#SBATCH --job-name="sift4g"

module load vcftools bcftools
VCF=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/00_AD1AD2AD4_n392/YUCFLAD2AD4_n392.AhDh.combined.bi.rehead.id.vcf
bcftools query -l $VCF | grep -v -E "AD2|AD4_mus_B40AD416|AD4_mus_P401AD4" > listn381.txt

vcftools \
  --vcf $VCF \
  --keep listn381.txt \
  --min-alleles 2 \
  --max-alleles 2 \
  --maf 0.001 \
  --recode \
  --recode-INFO-all \
  --out AD1AD4_n381_bi

bcftools view --exclude "F_PASS(GT='het')=1" --threads 5  AD1AD4_n381_bi.recode.vcf  -o  AD1AD4_n381_bi.recode.nofixed.vcf

module purge
module load openjdk/21.0.3_9-vngib7s

unset DISPLAY
java -Djava.awt.headless=true -jar SIFT4G_Annotator.jar \
    -c -i AD1AD4_n381_bi.recode.nofixed.vcf \
    -d /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/05_SIFT4G/SIFT4G_Create_Genomic_DB/test_files/gossypium_hirsutum/TEX2094v2 \
    -r outputsift \
    -t
```

### Filter the high confident deleterious mutation sites from `xls` file, and extract genotype only from n381 VCF, then combine SIFT4G score with VCF. Finally keep only the genotype of outgroup AD4 is `0/0`.
```
ml vcftools bcftools
awk 'NR>1 && /DELETERIOUS/ && !/Low confidence/ {print $1"\t"$2-1"\t"$2"\t"$13}' AD1AD4_n381_bi.recode.nofixed_SIFTannotations.xls > highconf_deleterious.bed

vcftools --vcf AD1AD4_n381_bi.recode.nofixed_SIFTpredictions.vcf \
  --bed highconf_deleterious.bed \
  --extract-FORMAT-info GT \
  --out highconf_gt
  
awk 'NR>1 && /DELETERIOUS/ && !/Low confidence/ { print $1":"$2"\t"$13}' AD1AD4_n381_bi.recode.nofixed_SIFTannotations.xls > sift_score_map.txt

awk 'BEGIN{OFS="\t"} FNR==NR{sift[$1]=$2; next} /^CHROM/{print $0,"SIFT_SCORE"; next}
{pos=$1":"$2; print $0,(pos in sift ? sift[pos] : "NA")}' sift_score_map.txt highconf_gt.GT.FORMAT > highconf_gt_with_score.tsv

rm sift_score_map.txt highconf_gt.GT.FORMAT

awk 'NR==1 || $383 == "0/0"' highconf_gt_with_score.tsv > highconf_gt_with_score_outgroup0.tsv
cut -f 1,2 highconf_gt_with_score_outgroup0.tsv > sift4g_highconf_gt_with_score_outgroup0_positionsONLY.txt
```
