# MSA BED Extractor

This README explains how to use the `msa_from_bed_multi.py` script to extract regions from multiple sequence alignments (MSAs) based on a BED file.

---

## Overview

`msa_from_bed_multi.py` is a Python script that:

- Extracts specific columns from chromosome-wise MSA files using a BED file.
- Handles gaps correctly in the MSA.
- Automatically detects FASTA files (`.fa` or `.fasta`) for each chromosome.
- Outputs sub-alignments per chromosome in FASTA format.

This is useful when working with multi-chromosome alignments and you want to focus on specific genomic regions.

---

## Requirements

- Python 3
- Biopython library

Install Biopython if not already installed:

```bash
pip install biopython
```

---

## Input Files

1. **BED File**

The BED file specifies the regions to extract. It should follow the standard BED format (tab-delimited, 0-based, end-exclusive):

```
chrom  start  end
Ah_01  0      5
Ah_01  10     12
Dh_01  0      3
Dh_02  4      8
```

- `chrom` must match the prefix of your FASTA filename.
- `start` is 0-based (first base is 0).
- `end` is exclusive (the base at this position is not included).

2. **MSA Directory**

- Directory containing per-chromosome MSA files in FASTA format (`.fa` or `.fasta`).
- Filenames should start with the chromosome name, e.g., `Ah_01.clean.fa` or `Dh_02.fasta`. The script ignores the extension and matches based on the prefix.

---

## Output

- Directory specified by `--out_dir` will contain one sub-alignment per chromosome.
- Files are named `<chrom>_sub_alignment.fasta`.
- Each sub-alignment contains all sequences from the original MSA but only columns corresponding to the BED-specified regions.

---

## Usage

```bash
python msa_from_bed_multi.py --bed <path_to_bed> --msa_dir <msa_directory> --out_dir <output_directory>
```

### Example

```bash
python msa_from_bed_multi.py --bed sites.bed --msa_dir 01_msa_dir/ --out_dir sub_alignment/
```

- `sites.bed` → BED file with regions for extraction.
- `01_msa_dir/` → directory containing chromosome FASTA MSAs.
- `sub_alignment/` → directory to save extracted sub-alignments.

---

## Notes

- Only `.fa` or `.fasta` files are considered. All other files are ignored.
- The script handles gaps in the MSA automatically.
- If a chromosome listed in the BED file does not have a corresponding FASTA file, the script will print a warning and skip that chromosome.
- The output sub-alignments maintain the same sequence order as the original MSAs.


---

## Troubleshooting

- **MSA file not found**: Ensure that the BED chromosome names match the prefix of the FASTA filenames.
- **Incorrect columns extracted**: Make sure the BED file uses 0-based, end-exclusive coordinates.
- **Biopython ImportError**: Install Biopython using `pip install biopython`.

---

## Author

Script by [Your Name].