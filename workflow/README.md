# Snakemake workflow — bulk RNA-seq (APOEhi vs Control)

Reproducible companion to [`Analysis.ipynb`](../Analysis.ipynb) and [`scripts/run_de_analysis.py`](../scripts/run_de_analysis.py).

## Pipeline

```text
rnaseq-01.xlsx → de_analysis → figures/ (volcano, PCA, heatmap, CSV)
```

| Rule | Outputs |
|------|---------|
| `de_analysis` | All files under `../figures/` |

## Quick start

From the project root, create the venv and install deps (includes Snakemake):

```bash
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python scripts/verify_data.py
```

Then run the workflow:

```bash
cd workflow
..\.venv\Scripts\python.exe -m snakemake -n
..\.venv\Scripts\python.exe -m snakemake -j 1 --cores 1
```

Dry run (`-n`) lists rules without executing. On Windows, use the venv Python path as shown above (no `--use-conda` required).

## Config

Edit `config/config.yaml` to change paths or document thresholds (`de.log2fc_threshold`, `de.fdr_threshold`, etc.). The script currently uses the same defaults; wire additional CLI flags in `run_de_analysis.py` if you extend the config.

## Verify data only

```bash
python scripts/verify_data.py
```
