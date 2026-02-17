## SnpEff

### build database for TX20094, annotate SNP, filter high_impact SNPs, and position for outgroup to be `0/0`
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=300G 
#SBATCH --time=7-02:30:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.snpEff.%J.out"
#SBATCH --job-name="snpEff"

Dir=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/06_snpeff
cd $Dir

module load snpeff/2017-11-24-zbh42lj

#cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/06_snpeff/data/TEX2094
#cp /lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.gtf genes.gtf
#cp /lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa sequences.fa
#cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/06_snpeff
#snpEff build -gtf22 -v TEX2094

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/04_vcfgerpcount/AD1AD4_n381.vcf.gz

snpEff TEX2094 $vcf > AD1_n381.AhDh.combined.bi.snpeff.vcf

module load openjdk/21.0.3_9-vngib7s
java -jar snpEff/SnpSift.jar filter "(ANN[*].IMPACT = 'HIGH')" AD1_n381.AhDh.combined.bi.snpeff.vcf > high_impact_n381.vcf

echo -e "CHROM\tPOS\t$(bcftools query -l high_impact_n381.vcf | tr '\n' '\t')" >  high_impact_n381.vcf.txt
bcftools query -f '%CHROM\t%POS\t[%GT\t]\n' high_impact_n381.vcf >> high_impact_n381.vcf.txt

awk 'NR>1 && $NF == "0/0"' high_impact_n381.vcf.txt | cut -f1,2 > high_impact_n381_outgroup0_positiononly.txt
```
