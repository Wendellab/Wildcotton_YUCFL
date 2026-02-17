Need Josh and Tony's input here for the reference genome assmbely and annotation


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
