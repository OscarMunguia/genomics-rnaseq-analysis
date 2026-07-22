# Downloading the RNA-seq dataset

This project uses supplementary FPKM expression data from:

> [DOI: 10.1016/j.canlet.2026.218605](https://doi.org/10.1016/j.canlet.2026.218605)

## Expected file

Place the Excel matrix in the project root as:

```text
rnaseq-01.xlsx
```

The workbook should contain FPKM columns for:

- `Control-1 (FPKM)` … `Control-3 (FPKM)`
- `APOEhi-1 (FPKM)` … `APOEhi-3 (FPKM)`
- `APOElo-1 (FPKM)` … `APOElo-3 (FPKM)`
- Gene annotation columns such as `GeneSymbol` and `Exonic.gene.sizes`

## How to obtain the data

1. Open the article page via the DOI link above.
2. Download the supplementary materials that include the THP-1 macrophage RNA-seq FPKM table.
3. Rename/copy the table to `rnaseq-01.xlsx` in this repository root if needed.
4. Confirm the file opens and contains the nine FPKM sample columns listed above.

If the repository already includes `rnaseq-01.xlsx`, you can skip the download and proceed to the analysis.

## Usage notes

- This dataset is used here for **educational and exploratory** reanalysis.
- Please cite the original publication and follow the journal / publisher terms for reuse.
- For publication-grade differential expression, prefer raw counts with DESeq2 or edgeR rather than FPKM + t-tests.
