#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G 
#SBATCH --time=10:00:00
#SBATCH --mail-user=weixuan@iastate.edu
#SBATCH --mail-type=ALL
#SBATCH --open-mode=append
#SBATCH --output="job.pixy_n90.%J.out"
#SBATCH --job-name="pixy_caculate"

mkdir -p output_temp
awk 'FNR>1 || NR==1' AD1380_domwild.*h*_dxy.txt > output_temp/AD1380_domwild.dxy.txt
awk 'FNR>1 || NR==1' AD1380_domwild.*h*_pi.txt > output_temp/AD1380_domwild.pi.txt
awk 'FNR>1 || NR==1' AD1380_domwild.*h*_fst.txt > output_temp/AD1380_domwild.fst.txt

cd output_temp

cut -f 1,2,5 AD1380_domwild.pi.txt > AD1380_domwild.pi2.txt
cut -f 1,2,3,6 AD1380_domwild.dxy.txt > AD1380_domwild.dxy2.txt
cut -f 1,2,3,6 AD1380_domwild.fst.txt > AD1380_domwild.fst2.txt

awk 'FNR > 1 {
    if($3 != "NA") {                 # Skip NA values
        group = $1"\t"$2             # Create a key using columns 1 and 2
        count[group]++               # Increment count
        sum[group] += $3             # Sum of column 3
        sumsq[group] += $3*$3        # Sum of squares
    }
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]
        variance = (sumsq[group] / count[group]) - (mean * mean)
        sd = (variance > 0) ? sqrt(variance) : 0
        print group, mean, sd
    }
}' AD1380_domwild.pi2.txt | sort > AD1380_domwild.pi3.txt



awk 'FNR > 1 {
    if($4 != "NA") {                 # Skip NA values in column 4
        group = $1"\t"$2"\t"$3       # Group by first 3 columns
        count[group]++
        sum[group] += $4
        sumsq[group] += $4*$4
    }
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]
        variance = (sumsq[group] / count[group]) - (mean * mean)
        sd = (variance > 0) ? sqrt(variance) : 0
        print group, mean, sd
    }
}' AD1380_domwild.dxy2.txt | sort > AD1380_domwild.dxy3.txt


awk 'FNR > 1 {
    if($4 != "NA") {                 # Skip NA values in column 4
        group = $1"\t"$2"\t"$3       # Group by first 3 columns
        count[group]++
        sum[group] += $4
        sumsq[group] += $4*$4
    }
}
END {
    for (group in sum) {
        mean = sum[group] / count[group]
        variance = (sumsq[group] / count[group]) - (mean * mean)
        sd = (variance > 0) ? sqrt(variance) : 0
        print group, mean, sd
    }
}' AD1380_domwild.fst2.txt | sort > AD1380_domwild.fst3.txt


awk '{key = $1 OFS $2; sum[key] += $4; count[key]++} END {for (k in sum) print k, sum[k] / count[k]}' AD1380_domwild.fst3.txt | sort > AD1380_domwild.fst4.txt
awk '{key = $1 OFS $2; sum[key] += $4; count[key]++} END {for (k in sum) print k, sum[k] / count[k]}' AD1380_domwild.dxy3.txt | sort > AD1380_domwild.dxy4.txt
awk '{key = $1 ; sum[key] += $3; count[key]++} END {for (k in sum) print k, sum[k] / count[k]}' AD1380_domwild.pi3.txt | sort > AD1380_domwild.pi4.txt
