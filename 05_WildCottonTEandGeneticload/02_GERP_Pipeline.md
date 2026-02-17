## Detelerious mutation scanning through genome using [gerp++](https://github.com/tvkent/GERPplusplus) 
##### reference to: https://github.com/HuffordLab/NAM-genomes/tree/master/gerp 

### 1. Prepare the input alignment file for `grepcol`

#### 1.1 select the reference genome using site PubPlant (https://www.plabipd.de/pubplant_main.html), based on phylogeny (https://www.plabipd.de/pubplant_cladogram1.html) and published date.
 
| No. | Family | Species | Accession | Link |
|:------- |:------- |:------- | :------- | :------- |
| 1 | Malvaceae | Gossypium hisutum | TEX2094 v.2 ISU  | ubpublished  |
| 2 | Fabaceae | Medicago truncatula | MtrunA17r5.0-ANR  | [NCBI GCF_003473485.1 ](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003473485.1/)  |
| 3 | Brassicaceae | Arabidopsis thaliana | TAIR10.1  | [NCBI GCF_000001735.4](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001735.4/)  |
| 4 | Asteraceae | Helianthus annuus | HanXRQr2.0-SUNRISE  | [NCBI GCF_002127325.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002127325.2/)   |
| 5 | Poaceae | Oryza sativa Japonica Group | AGIS1.0  | [NCBI GCA_034140825.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_034140825.1/)  |
| 6 | Bromeliaceae | Ananas comosus | ASM154086v1  | [NCBI GCF_001540865.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001540865.1/)   |
| 7 | Apiaceae | Daucus carota | DH1 v3.0  | [NCBI GCF_001625215.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001625215.2/)  |
| 8 | Poaceae | Zea mays | Zm-B73-REFERENCE-NAM-5.0  | [NCBI GCA_902167145.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_902167145.1/) |
| 9 | Solanaceae | Solanum lycopersicum | SLM_r2.1 | [NCBI GCF_036512215.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_036512215.1/)  |
| 10 | Theaceae | Camellia pitardii | ASM5162307v1  | [NCBI GCA_051623075.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_051623075.1/)  |
| 11 | Vitaceae | Vitis davidii | V112.hap2_v1.0 | [NCBI GCA_044588485.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_044588485.1/)  |
| 12 | Salicaceae | Salix dunnii | FNU-M-1-Hap-a | [NCBI GCA_040801865.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_040801865.1/)  |
| 13 | Rutaceae | Citrus sinensis | DVS_A1.0  | [NCBI GCA_022201045.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_022201045.2/)  |
| 14 | Malvaceae | Hibiscus yunnanensis | ASM4900470v1  | [NCBI GCA_049004705.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_049004705.1/) |
| 15 | Malvaceae | Gossypium raimondii | ASM2569854v1  | [NCBI GCA_025698545.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_025698545.1/) |

Download the genoems from NCBI use `datasets`
```
https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/

curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat

[weixuan@nova23-16 ncbidatasets]$ nl genome_accessions.txt
     1  GCF_003473485.1
     2  GCF_000001735.4
     3  GCF_002127325.2
     4  GCA_034140825.1
     5  GCF_001540865.1
     6  GCF_001625215.2
     7  GCF_902167145.1
     8  GCF_000226075.1
     9  GCA_051623075.1
    10  GCA_044588485.1
    11  GCA_040801865.1
    12  GCF_022201045.2
    13  GCA_049004705.1
    14  GCF_025698545.1

mkdir -p genomes
while read acc; do
    echo "Downloading $acc ..."
    ./datasets download genome accession $acc --filename genomes/${acc}.zip --include genome
    if [ $? -ne 0 ]; then
        echo "ERROR: $acc failed to download!"
    else
        echo "$acc downloaded successfully."
        unzip -o genomes/${acc}.zip -d genomes/${acc}
    fi
done < genome_accessions.txt
```

Remove the plastid DNA and mtDNA from the reference genomes, and make sure all their reference genomes are TE masked, and remove all the unplaced scaffolds.
```
seqkit mutate -w0 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_034140825.1/ncbi_dataset/data/GCA_034140825.1/GCA_034140825.1_AGIS1.0_genomic.fna > Oryza_sative.fasta

awk '/^>/{flag=0} 
     /^>CP161459.1 Salix dunnii isolate FNU-M-1 mitochondrion, complete genome/{flag=1; next}
     /^>CP161458.1 Salix dunnii isolate FNU-M-1 chloroplast, complete genome/{flag=1; next}
     flag==0' /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_040801865.1/ncbi_dataset/data/GCA_040801865.1/GCA_040801865.1_FNU-M-1-Hap-a_genomic.fna | \
seqkit mutate -w0 > Salix_dunnii.fasta

seqkit mutate -w0 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_049004705.1/ncbi_dataset/data/GCA_049004705.1/GCA_049004705.1_ASM4900470v1_genomic.fna > Hibiscus_yunnanensis.fasta

seqkit mutate -w0 /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_049004705.1/ncbi_dataset/data/GCA_049004705.1/GCA_049004705.1_ASM4900470v1_genomic.fna > Hibiscus_yunnanensis.fasta

awk '/^>/{flag=0}
     /^>JBPSJW/{flag=1; next}
     flag==0' /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_051623075.1/ncbi_dataset/data/GCA_051623075.1/GCA_051623075.1_ASM5162307v1_genomic.fna | \
	 seqkit mutate -w0 >  Camellia_pitardii.fasta
	 
awk '/^>/{flag=0} 
     /^>NC_037304.1 Arabidopsis thaliana ecotype Col-0 mitochondrion, complete genome/{flag=1; next}
     /^>NC_000932.1 Arabidopsis thaliana chloroplast, complete genome/{flag=1; next}
     flag==0' /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_000001735.4/ncbi_dataset/data/GCF_000001735.4/GCF_000001735.4_TAIR10.1_genomic.fna | \
seqkit mutate -w0 > Arabidopsis_thaliana.fasta

#grep '>' /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_000226075.1/ncbi_dataset/data/GCF_000226075.1/GCF_000226075.1_SolTub_3.0_genomic.fna  | grep 'complete'

awk '/^>/{flag=0}
     /unplaced/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_001540865.1/ncbi_dataset/data/GCF_001540865.1/GCF_001540865.1_ASM154086v1_genomic.fna | \
	 seqkit mutate -w0 > Ananas_comosus.fasta

awk '/^>/{flag=0}
     /mitochondrion/{flag=1; next}
	 /plastid/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_001625215.2/ncbi_dataset/data/GCF_001625215.2/GCF_001625215.2_DH1_v3.0_genomic.fna | \
	 seqkit mutate -w0 > Daucus_sativus.fasta
	 
awk '/^>/{flag=0}
     /unplaced/{flag=1; next}
	 /mitochondrion/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_003473485.1/ncbi_dataset/data/GCF_003473485.1/GCF_003473485.1_MtrunA17r5.0-ANR_genomic.fna | \
	 seqkit mutate -w0 > Medicago_truncatula.fasta

awk '/^>/{flag=0}
	 /mitochondrion/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_022201045.2/ncbi_dataset/data/GCF_022201045.2/GCF_022201045.2_DVS_A1.0_genomic.fna | \
	 seqkit mutate -w0 > Citrus_sinensis.fasta
	 
awk '/^>/{flag=0}
     /unplaced/{flag=1; next}
	 /mitochondrion/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_025698545.1/ncbi_dataset/data/GCF_025698545.1/GCF_025698545.1_ASM2569854v1_genomic.fna | \
	 seqkit mutate -w0 > Gossypium_raimondii.fasta

awk '/^>/{flag=0}
	 /mitochondrion/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/data/GCF_036512215.1/GCF_036512215.1_SLM_r2.1_genomic.fna | \
	 seqkit mutate -w0 > Solanum_lycopersicum.fasta

awk '/^>/{flag=0}
     /unplaced/{flag=1; next}
	 /chloroplast/{flag=1; next}
	 /mitochondrion/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_002127325.2/data/GCF_002127325.2/GCF_002127325.2_HanXRQr2.0-SUNRISE_genomic.fna | \
	 seqkit mutate -w0 > Helianthus_annuus.fasta

awk '/^>/{flag=0}
     /unplaced/{flag=1; next}
	 /mitochondrion/{flag=1; next}
	 /chloroplast/{flag=1; next}
     flag==0'  /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCF_902167145.1/data/GCF_902167145.1/GCF_902167145.1_Zm-B73-REFERENCE-NAM-5.0_genomic.fna | \
	 seqkit mutate -w0 > Zea_mays.fasta	 

awk '/^>/{flag=0}
     /JBIHUZ/{flag=1; next}
     flag==0'   /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/ncbidatasets/genomes/GCA_044588485.1/data/GCA_044588485.1/GCA_044588485.1_V112.hap2_v1.0_genomic.fna | \
	 seqkit mutate -w0 > Vitis_davidii.fasta
```
#
#### 1.2 Send the selected taxa into TimeTree (https://timetree.org/) to build a newick phylogeny with manual adding Malvaceae interspecific relationships
Cactus requires the input phylogeny to be strict binary: A binary tree is a tree where every internal node has exactly two children.
A rooted tree can be problematic, because this creates a multifurcation at the root, which sonLib.nxnewick cannot parse (it expects strictly binary or fully bracketed Newick).
I think the easiest way to fix this is Chatgpt, which gives you a tree to start off.

![alt text](https://github.com/WeixuanPlant/JingPangenome/blob/main/Geneticload/01_GenomeSelection/trailtree2.nwk.jpg)
```
(((Solanum_lycopersicum,(Daucus_carota,Helianthus_annuus)),Camellia_pitardii),((Vitis_davidii,((Arabidopsis_thaliana,((Gossypium_hirsutum,Gossypium_raimondii),Hibiscus_yunnanensis)),(Citrus_sinensis,(Salix_dunnii,Medicago_truncatula)))),(Ananas_comosus,(Oryza_sativa,Zea_mays))));
```
# 
#### 1.3 Mask the repeat regions in TEX2094 reference genome
```
# Install repeatmasker via micromamba

cd /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs

module load micromamba/1.4.2-7jjmfkf
eval "$(micromamba shell hook --shell=bash)"
micromamba create \
  --root-prefix /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/ \
  -n repeatmasker_env \
  -c conda-forge -c bioconda \
  python=3.10 repeatmasker=4.2.0 perl rmblast=2.14 trf=4.09 wget

micromamba activate /lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/envs/repeatmasker_env

# Soft-mask (Not to NNNNN but to lowercase) of repeat regions of TEX2094 reference genome

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/TEX2094rm_out
ref=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa
refTE=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/AD1.TX2094.v2.fa.mod.EDTA.TElib.fa
output=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes
thr=8

RepeatMasker -lib $refTE -pa $thr -xsmall -gff -dir $output/rm_out $ref
```
### 2. Prepare the input alignment file for `grepcol`

#### 2.1 build the multispecies alignment using `cactus`

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=30
#SBATCH --mem=1000G
#SBATCH --time=5-01:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.cactus.%J.out"
#SBATCH --job-name="cactus"

# Set directories
CACTUS_HOME=/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/cactus
WORKDIR=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign
JOBSTORE=$WORKDIR/jb
SEQFILE=$WORKDIR/cactus.txt
HAL_OUT=$WORKDIR/n15ref.hal
LOGFILE=$WORKDIR/n15ref.log

# Make sure jobStore exists
# Run Cactus inside Singularity

singularity exec --home $CACTUS_HOME \
$CACTUS_HOME/cactus_v2.9.9.sif cactus \
$JOBSTORE \
$SEQFILE \
$HAL_OUT \
--logFile $LOGFILE \
--maxCores 30
```

#### 2.2 Converting the aligned `HAL` format into `MAF` format

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=1-20:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --output="job.cactus.%J.out"
#SBATCH --job-name="cactus"

module load apptainer

apptainer exec docker://quay.io/comparative-genomics-toolkit/cactus:v3.0.0 cactus-hal2maf \
./js n15ref.hal cactus_maf/n15ref_target.maf.gz \
--refGenome Gossypium_hirsutum \
--noAncestors \
--coverage \
--onlyOrthologs \  #  --onlyOrthologs       Run hal2maf with --onlyOrthologs option which attempts to keep only duplications that are also separate in ancestor (default: False)
--outType single \ 
--maxCores 50 \
--targetGenomes Solanum_lycopersicum,Daucus_carota,Helianthus_annuus,Camellia_pitardii,Vitis_davidii,Arabidopsis_thaliana,Gossypium_raimondii,Hibiscus_yunnanensis,Citrus_sinensis,Salix_dunnii,Medicago_truncatula,Ananas_comosus,Oryza_sativa,Zea_mays

# --outType {raw,norm,single,consensus} [{raw,norm,single,consensus} ...]
   #     Select which kind of postprocessing to apply to the hal2maf output. Multiple selections allowed. raw: return hal2maf output as-is; norm: run taffy normalization to merge adjacent blocks where possible;
   #     single: heuristically choose single, most similar homolog for each species for each (normalized) block using mafDuplicateFilter; consensus: squish all duplicate rows for each (normalized) block into a
   #     single conensus row using maf_stream. [default=norm] (default: ['norm'])
				   
```

#### 2.3 Spliting the MAF into individual chromosomes (via `mafSplit`), then use `msa_view` convert MAF to final `FASTA` format 

##### 2.3.1 install mafSplit and msa_view
```
wget https://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/mafSplit ./
chmod +x ./mafSplit

wget http://www.netlib.org/clapack/clapack.tgz
tar -xvzf clapack.tgz
cd CLAPACK-3.2.1/
cp make.inc.example make.inc && make f2clib && make blaslib && make lib

wget http://compgen.cshl.edu/phast/downloads/phast.v1_5.tgz
tar -xvzf phast.v1_5.tgz
cd phast/src/
make CLAPACKPATH=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/CLAPACK-3.2.1
```

#### 2.3.2 Split the MAF output into individual chromosomes AND split the reference genome of Gossypium hirsutum TEX2094 into individual chromosomes for msa_view
```
cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf

#seperate by chromosome from maf file 
/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/usuc/mafSplit -byTarget dummy.bed ./ ../n15ref_target.maf.gz -useFullSequenceName

MSA_VIEW=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/phast/bin/msa_view
REF=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/TEX2094rm_out/AD1.TX2094.v2.mask.fa

# Split the reference FASTA into one file per sequence (chromosome)
seqkit split -s 1 -O 00_ref_chroms/ $REF

chrom_names=(Ah_{01..13} Dh_{01..13})
cd 00_ref_chroms
i=0
for f in AD1.TX2094.v2.mask.part_*.fa; do
    mv "$f" "${chrom_names[i]}.fa"
    i=$((i+1))
done
```

#### 2.3.4 msa_view first convert MAF into FASAT using reference genomes `--refseq`; then we cleaned the fasta format by (1) removing the space in fasta header; (2) change all * into -. 
```
cd ../
mkdir -p 01_msa_dir
REF_DIR=00_ref_chroms  # per-chromosome reference folder
for chr in Ah_{01..13} Dh_{01..13}; do
    echo "Processing $chr ..."

    # Convert MAF to FASTA
    $MSA_VIEW ./${chr}.maf \
        --in-format MAF \
        --out-format FASTA \
        -G 1 -f \
        --refseq ${REF_DIR}/${chr}.fa \
        > ${chr}.fa

    # Clean FASTA: replace * with - and trim spaces in headers
	sed 's/^> */>/' ${chr}.fa | sed 's/\*/-/g' | awk '/^>/ {gsub(/ +$/, "", $0)} {print}' > 01_msa_dir/${chr}.clean.fa
    echo "$chr done."
done
```

### 3. Prepare the input alignment file for `grepcol`

#### 3.1 Build a netural tree use 4d sites
```
module load py-biopython/1.81-py310-wt5scj7
db4_bed=/lustre/hdd/LAS/jfw-lab/weixuan/00_RefTX0294/degenotate-out/degeneracy-all-sites-sort.bed

# Make sure your BED coordinates are 0-based for bedtools
python msa_from_bed_multi.py --bed $db4_bed --msa_dir 01_msa_dir/ --out_dir 02_sub_alignment/

mkdir -p iqtree

seqkit concat -o concatenated.fasta ../*.fasta

module load iqtree2
iqtree2 -s concatenated.fasta -m GTR -nt AUTO
```

#### 3.2 Caculate RS (Rejected subsitution score) for each site of reference genome TEX2094 
```
cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/01_msa_dir

neturaltree=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/02_sub_alignment/iqtree/concatenated.fasta.treefile

mkdir -p 03_gerpscore


for chr in Ah_{01..13} Dh_{01..13}; do
    echo "Processing $chr ..."
	
	sed 's/^> */>/' ${chr}.fa | sed 's/\*/-/g' | awk '/^>/ {gsub(/ +$/, "", $0)} {print}' > 01_msa_dir/${chr}.clean.fa
	
	/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/gerp++KRT/gerpcol -t $neturaltree -f 01_msa_dir/${chr}.clean.fa -a  -e Gossypium_hirsutum -j -v
	
	mv 01_msa_dir/${chr}.clean.fa.rates 03_gerpscore/
	
	/lustre/hdd/LAS/jfw-lab/weixuan/00_BioinformaticTools/gerp++KRT/gerpelem -f 03_gerpscore/${chr}.clean.fa.rates
	
    # Make per-position BED from rates
    awk -v chr="$chr" 'BEGIN {OFS="\t"} {print chr, NR-1, NR, $1, $2}' 03_gerpscore/${chr}.clean.fa.rates > 03_gerpscore/${chr}.bed

    # BED of positive scores
    awk '$5 > 0' 03_gerpscore/${chr}.bed > 03_gerpscore/${chr}.pos.bed

    # BED of elements
    awk -v chr="$chr" 'BEGIN {OFS="\t"} {print chr, $2-1, $3, $4, $5, $6, $7, $8}' \
        03_gerpscore/${chr}.clean.fa.rates.elems > 03_gerpscore/${chr}.elems.bed
    sort -k1,1 -k2,2n 03_gerpscore/${chr}.elems.bed -o 03_gerpscore/${chr}.elems.bed

	echo "$chr done."
done
```

#### 3.3 Adding the RS scores as an annotation to the candidate VCF
```
mkdir -p 04_vcfgerpcount

cat 03_gerpscore/*.pos.bed > 04_vcfgerpcount/all_GERP_pos.bed

cd /lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/03_gerpRefGenomes/01_cactusalign/cactus_maf/04_vcfgerpcount

module load bedtools2/2.31.1-py311-6kemgt3

vcf=/lustre/hdd/LAS/jfw-lab/weixuan/08_YUCFL_popgene/01_mergedVCFs/04_AD1AD2AD4_n380/YUCFLAD2AD4_n380.AhDh.combined.bi.rehead.id.vcf
gerpbed=all_GERP_pos.bed

bedtools intersect -a $vcf -b $gerpbed -header > YUCFLAD2AD4_n380.varintgerp.vcf
```
