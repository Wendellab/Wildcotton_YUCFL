## Selection signal 

### Prepare input 
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=150G
#SBATCH --time=1-00:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.YUC_C.%J.out"
#SBATCH --job-name="YUC_C"
#SBATCH --array=1-13

#cat /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/list_AD1_Cultivar_n109.txt \
#/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/00_interval/list_gvcf_YUC_n158.txt |\
#grep -v -E "RiCh|RiCa" > list_gvcf_YUCCultivar_n227.txt

seq=$(printf %02d ${SLURM_ARRAY_TASK_ID})
echo "$seq"

Dir=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/07_selectivesweep
ref=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
thr=30 #NUMBER_THREADS

output1=YUCCultivar_n227

cd $Dir
#mkdir $output1

module load sentieon-genomics/202308.02-e2gz6fb
export SENTIEON_LICENSE=reimu.las.iastate.edu:8990

#joint SNP calling:

cat list_gvcf_YUCCultivar_n227.txt | sentieon driver --interval Ah_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode variant $TMPDIR/$output1.Ah_$seq.vcf -
mv $TMPDIR/$output1.Ah_$seq.vcf* $Dir

cat list_gvcf_YUCCultivar_n227.txt | sentieon driver --interval Dh_$seq -t $thr -r $ref --algo GVCFtyper --emit_mode variant $TMPDIR/$output1.Dh_$seq.vcf -
mv $TMPDIR/$output1.Dh_$seq.vcf* $Dir

#echo "Filtering VCF file using vcftools and rename samples using bcftools"
cd $Dir

ml vcftools bcftools

for g in Ah Dh; do
  vcftools --vcf $output1.${g}_$seq.vcf \
    --remove-indels \
    --min-meanDP 10 --max-meanDP 100 \
    --mac 2 \
    --min-alleles 2 --max-alleles 2 \
    --max-missing 1.0 \
    --recode --recode-INFO-all \
    --out $output1.${g}_$seq.bisnps
	

bcftools view --exclude "F_PASS(GT='het')=1" \
 $output1.${g}_$seq.bisnps.recode.vcf \
 -o  $output1.${g}_$seq.bisnps.nofixed.vcf

done



for chr_prefix in Ah Dh; do
  for i in $(seq -w 01 13); do
    input="YUCCultivar_n227.${chr_prefix}_${i}.bisnps.nofixed.vcf"
    output="YUCCultivar_n227.${chr_prefix}_${i}.bisnps.nofixed.rename.vcf"
    
    bcftools reheader -s rename_YUCFLAD2AD4_n491.txt $input -o $output
  done
done

```
#
### XPCLR
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --time=2-20:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.xpclr.%J.out"
#SBATCH --job-name="xpclr"
#SBATCH --array=1-26


ml bcftools
bcftools query -l YUCCultivar_n227.Ah_01.bisnps.nofixed.rename.vcf | grep "Cultivar"  > list_Cultivar.txt
bcftools query -l YUCCultivar_n227.Ah_01.bisnps.nofixed.rename.vcf | grep -v "Cultivar"  > list_Wild.txt

# Load modules / activate environment
module load py-numpy/1.24.4-py38-6fb52jw
module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/xpclr_env

Dir=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/07_selectivesweep
mkdir -p 01_xpclr_results

# Map SLURM array to chromosome
if [ ${SLURM_ARRAY_TASK_ID} -le 13 ]; then
    chr="Ah_$(printf "%02d" ${SLURM_ARRAY_TASK_ID})"
else
    dh_id=$((SLURM_ARRAY_TASK_ID - 13))
    chr="Dh_$(printf "%02d" $dh_id)"
fi

echo "Running chromosome $chr"

xpclr --format vcf --input YUCCultivar_n227.${chr}.bisnps.nofixed.rename.vcf \
--samplesA $Dir/list_Cultivar.txt --samplesB $Dir/list_Wild.txt \
--ld 0.7 \
--minsnps 10 \
--size 100000 \
--step 20000 \
--chr $chr \
--verbose 10 \
--out ./01_xpclr_results/Cultivar_Wild_${chr}
```
#
### VCFtools
```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100G
#SBATCH --time=2-20:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.vcftools.%J.out"
#SBATCH --job-name="vcftools"
#SBATCH --array=1-26

module load vcftools

mkdir -p 02_Fst_Pi_results

# Map SLURM array to chromosome
if [ ${SLURM_ARRAY_TASK_ID} -le 13 ]; then
    chr="Ah_$(printf "%02d" ${SLURM_ARRAY_TASK_ID})"
else
    dh_id=$((SLURM_ARRAY_TASK_ID - 13))
    chr="Dh_$(printf "%02d" $dh_id)"
fi

vcftools --vcf YUCCultivar_n227.${chr}.bisnps.nofixed.rename.vcf \
  --chr $chr \
  --weir-fst-pop list_Cultivar.txt \
  --weir-fst-pop list_Wild.txt \
  --fst-window-size 100000 \
  --fst-window-step 20000 \
  --out 02_Fst_Pi_results/fst_Cultivar_Wild_$chr

# Pi per population
vcftools --vcf YUCCultivar_n227.${chr}.bisnps.nofixed.rename.vcf \
  --chr $chr \
  --window-pi 100000 \
  --window-pi-step 20000 \
  --keep list_Cultivar.txt \
  --out 02_Fst_Pi_results/pi_Cultivar_$chr

vcftools --vcf YUCCultivar_n227.${chr}.bisnps.nofixed.rename.vcf \
  --chr $chr \
  --window-pi 100000 \
  --window-pi-step 20000 \
  --keep list_Wild.txt \
  --out 02_Fst_Pi_results/pi_Wild_$chr
```
