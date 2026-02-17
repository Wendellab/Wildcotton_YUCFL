## SWIFTE

### 1. Subset reads data and run SWIFTE 
```
#!/bin/bash
#SBATCH --nodes=1 
#SBATCH --cpus-per-task=20 
#SBATCH --mem=300G 
#SBATCH --time=2-20:30:00
#SBATCH --open-mode=append
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output=joblog/job.sam.%A_%a.out 
#SBATCH --job-name="sam"
#SBATCH --array=1-116%50

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/02_swift

ref=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
refTE=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa.mod.EDTA.TElib.fa
Dir=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/02_swift/01_subsetreads

thr=20
file1=$(cat samplelist_subset_n116.txt | sed -n ${SLURM_ARRAY_TASK_ID}p)
file2=$(echo $file1 | sed 's/R1[.]fq/R2\.fq/')
name=$(basename $file1 _clean_R1.fq.gz)

echo "the file1 is " $file1
echo "the file2 is " $file2
echo "the sample is " $name

seqtk=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/02_swift/seqtk/seqtk

$seqtk sample -s 100  $file1 131000000 | gzip > $TMPDIR/$name.R1.fq.gz
$seqtk sample -s 100  $file2 131000000 | gzip > $TMPDIR/$name.R2.fq.gz

mv $TMPDIR/$name.R*.fq.gz $Dir

module load bwa
bwa mem -M -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" -t $thr $ref $Dir/$name.R1.fq.gz $Dir/$name.R2.fq.gz > $TMPDIR/$name.ref.sam
bwa mem -M -R "@RG\tID:$name\tSM:$name\tPL:ILLUMINA" -t $thr $refTE $Dir/$name.R1.fq.gz $Dir/$name.R2.fq.gz > $TMPDIR/$name.te.sam

cd $TMPDIR
/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/SWIFTE/SWIFTEv35 $TMPDIR/$name.ref.sam $TMPDIR/$name.te.sam
mv /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/SWIFTE/$name*.txt /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/04_TE/02_swift/01_subsetreads/02_SWIFTEoutput/
```

### 1. Subset reads data and run SWIFTE 
```
#!/bin/bash

cd 02_SWIFTEoutput
for i in *.txt; do  echo $i; bash /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/SWIFTE/clean_SWIF-TE.sh -i $i -c 3; done
# we subset the reads here, so lets not use default -5 to filter the reads count for TE.

mkdir -p temp
for file in *.bed; do
    base="${file%%.*}"; sep="_"; [[ "$base" == *-* ]] && sep="-"
    prefix="${base%%$sep*}"; suffix="${base#*$sep}"; [[ "$prefix" == "$suffix" ]] && suffix=""
    case "$prefix" in
        CeCo) mapped="YUC_CeCo" ;; CeDo) mapped="YUC_CeDo" ;; CeDr) mapped="YUC_CeDr" ;; CePr) mapped="YUC_CePr" ;;
        CPH) mapped="FL_CPH" ;; CPT) mapped="FL_CPT" ;; FMY) mapped="FL_FMY" ;; NP) mapped="FL_NP" ;;
        OPB2) mapped="FL_OPB2" ;; OPB4) mapped="FL_OPB4" ;; OPB5) mapped="FL_OPB5" ;;
        Ph) mapped="PR_Ph" ;; PK) mapped="FL_PK" ;; Pop1) mapped="FL_Pop1" ;; Pop2) mapped="FL_Pop2" ;; Pop3) mapped="FL_Pop3" ;;
        RBD) mapped="FL_RBD" ;; RBT) mapped="FL_RBT" ;; RiCa) mapped="YUC_RiCa" ;; RiCh) mapped="YUC_RiCh" ;; RNRB) mapped="FL_RNRB" ;;
        SiPa) mapped="YUC_SiPa" ;; SiPr) mapped="YUC_SiPr" ;; SR) mapped="FL_SR" ;; TC) mapped="FL_TC" ;; VKPA) mapped="FL_VKPA" ;; YUC) mapped="GD_YUC" ;; *) mapped="$prefix" ;;
    esac
    label="$mapped"; [[ -n "$suffix" && "$suffix" != "$prefix" ]] && label="${mapped}_$suffix"
    awk -v name="$label" '{print $0, name}' "$file" > "temp/$base.bed"
done

cd temp/
cat * > final_n116.bed
```
