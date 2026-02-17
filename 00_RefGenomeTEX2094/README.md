Need Josh and Tony's input here for the reference genome assmbely and annotation

##
### BUSCO score
```
# BUSCO version is: 5.8.2 
# The lineage dataset is: eudicotyledons_odb12 (Creation date: 2025-04-11, number of genomes: 76, number of BUSCOs: 2805)
# Summarized benchmarking in BUSCO notation for file /local/scratch/maa146/5557867/AD1.TX2094.v2.fa
# BUSCO was run in mode: euk_genome_min
# Gene predictor used: miniprot

	***** Results: *****

	C:99.3%[S:2.1%,D:97.1%],F:0.4%,M:0.4%,n:2805,E:0.7%	   
	2784	Complete BUSCOs (C)	(of which 20 contain internal stop codons)		   
	59	Complete and single-copy BUSCOs (S)	   
	2725	Complete and duplicated BUSCOs (D)	   
	11	Fragmented BUSCOs (F)			   
	10	Missing BUSCOs (M)			   
	2805	Total BUSCO groups searched		   

Assembly Statistics:
	26	Number of scaffolds
	48	Number of contigs
	2289074707	Total length
	0.000%	Percent gaps
	107 MB	Scaffold N50
```

## Count chromosome length, TE content and gene content from index, gff and gtf files.
### chromosome length
```
[weixuan@nova-login-1 00_RefTX0294]$ awk '
$1 ~ /^Ah_/ {at += $2}
$1 ~ /^Dh_/ {dt += $2}
END {
    print "At genome size:", at
    print "Dt genome size:", dt
    print "Total genome size:", at + dt
}
' AD1.TX2094.v2.fa.fai
At genome size: 1435390769
Dt genome size: 853683938
Total genome size: 2289074707
```

### TE
```
module load bedtools2/2.31.1-py311-6kemgt3

[weixuan@nova-login-1 00_RefTX0294]$ grep -v '^#' AD1.TX2094.v2.fa.mod.EDTA.TEanno.gff3 \
| awk '$1~/^Ah_/{print $1"\t"$4-1"\t"$5}' \
| sort -k1,1 -k2,2n \
| bedtools merge \
| awk '{sum+=$3-$2} END{printf "At TE: %d bp, %.2f%%\n", sum, sum/1435390769*100}'
At TE: 1101407703 bp, 76.73%

[weixuan@nova-login-1 00_RefTX0294]$ grep -v '^#' AD1.TX2094.v2.fa.mod.EDTA.TEanno.gff3 \
| awk '$1~/^Dh_/{print $1"\t"$4-1"\t"$5}' \
| sort -k1,1 -k2,2n \
| bedtools merge \
| awk '{sum+=$3-$2} END{printf "Dt TE: %d bp, %.2f%%\n", sum, sum/853683938*100}'
Dt TE: 521810461 bp, 61.12%
```

### Gene
```
[weixuan@nova-login-1 00_RefTX0294]$ awk '
$3=="gene" {
    if($1~/^Ah_/) at++
    else if($1~/^Dh_/) dt++
}
END {
    print "At genes:", at
    print "Dt genes:", dt
}
' AD1.TX2094.v2.gtf
At genes: 41720
Dt genes: 43651


[weixuan@nova-login-1 00_RefTX0294]$ grep -v '^#' AD1.TX2094.v2.gtf \
| awk '$3=="gene" && $1~/^Ah_/{print $1"\t"$4-1"\t"$5}' \
| sort -k1,1 -k2,2n \
| bedtools merge \
| awk '{sum+=$3-$2} END{printf "At genes: %d bp, %.2f%% of genome\n", sum, sum/1435390769*100}'
At genes: 89194706 bp, 6.21% of genome
[weixuan@nova-login-1 00_RefTX0294]$ grep -v '^#' AD1.TX2094.v2.gtf \
| awk '$3=="gene" && $1~/^Dh_/{print $1"\t"$4-1"\t"$5}' \
| sort -k1,1 -k2,2n \
| bedtools merge \
| awk '{sum+=$3-$2} END{printf "Dt genes: %d bp, %.2f%% of genome\n", sum, sum/853683938*100}'
Dt genes: 92910453 bp, 10.88% of genome
```
