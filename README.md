# Prioritizing Causal Variants at an Estimated BMD GWAS Locus with AlphaGenome

This repository provides a reproducible pipeline for using [AlphaGenome](https://www.nature.com/articles/s41586-025-10014-0) (Google DeepMind) to prioritize candidate causal variants at the **FHL3** locus on chromosome 1p34.3, a genome-wide significant locus for **estimated bone mineral density (eBMD)** identified by [Morris et al. (2018)](https://doi.org/10.1038/s41588-018-0302-x).

Starting from GWAS summary statistics and a lead SNP, the pipeline identifies LD proxies, lifts coordinates to hg38, runs AlphaGenome variant effect predictions across all GTEx tissues and output modalities, and produces ranked variant prioritization tables suitable for downstream functional follow-up.

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 01 | `01_prepare_snps_for_alphagenome.py` | Extract genome-wide significant SNPs in a +/-200 kb window around lead SNP **rs28428561**, identify independent signals by LD clumping (r^2 >= 0.10), and collect LD proxies (r^2 >= 0.80) using 1000 Genomes Phase 3 EUR samples. All coordinates are GRCh37/hg19. |
| 01B | `01B_liftover_hg19_hg38.py` | Lift over proxy variant coordinates from hg19 to GRCh38 using `pyliftover`. Handles reverse-complement alleles on strand-flipped liftover hits. |
| 02 | `02_run_alphagenome_alltiss_allmodal.py` | Run AlphaGenome `predict_variant()` for every proxy variant across **all GTEx tissues** and **all available output modalities** (RNA-seq, chromatin accessibility, histone marks, TF binding, contact maps, splicing). Produces a long-format parquet and a collapsed per-variant/tissue/modality summary CSV. Includes checkpointing for resume support. |
| 02B | `02B_recollapse_from_long.py` | Utility to re-collapse the long parquet output into the summary CSV without re-calling the API. |
| 03 | `03_score_variant.py` | Run AlphaGenome's built-in `RECOMMENDED_VARIANT_SCORERS` via `score_variant()` to produce calibrated per-gene and per-track effect scores comparable across modalities. |
| 04 | `04_rank_prioritization.py` | Aggregate variant scorer rankings into overall SNP prioritization: per-scorer rankings, cross-scorer consensus ranking, and tissue-specific rankings. |
| 05 | `05_collapse_eqtl_sqtl.py` | Collapse scores to per-SNP-per-gene tables for eQTL-like (RNA-seq `GeneMask` scorers) and sQTL-like (splicing scorer/output) signals. |

## Required Support Files

Several external data files are needed to run the pipeline. They are **not included** in this repository due to size and licensing. Below is how to obtain each one.

### 1. Morris et al. (2018) eBMD GWAS Summary Statistics

The pipeline expects the UK Biobank British eBMD summary statistics file:

```
Biobank2-British-Bmd-As-C-Gwas-SumStats.txt
```

**How to obtain:**

Download from the [GEFOS Consortium](http://www.gefos.org/?q=content/data-release-2018) or the [GWAS Catalog](https://www.ebi.ac.uk/gwas/studies/GCST006979). The file is distributed as part of the `Morrisetal2018.NatGen.SumStats.tar` archive.

> Morris, J.A. *et al.* An atlas of genetic influences on osteoporosis in humans and mice. *Nature Genetics* **51**, 258-266 (2019). https://doi.org/10.1038/s41588-018-0302-x

After downloading, update the `SUMSTATS` path in `01_prepare_snps_for_alphagenome.py` to point to the extracted file.

### 2. hg19-to-hg38 LiftOver Chain File

The liftover step requires the UCSC chain file:

```
hg19ToHg38.over.chain.gz
```

**How to obtain:**

```bash
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

Place it in a local directory and update the `LIFTOVER_DIR` / `CHAIN_19_TO_38` path in `01B_liftover_hg19_hg38.py`.

### 3. 1000 Genomes Phase 3 Panel and VCF Data

Step 01 automatically downloads the sample panel file:

```
integrated_call_samples_v3.20130502.ALL.panel
```

from the [1000 Genomes FTP](https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/). The per-chromosome VCF files are streamed remotely via `bcftools` from the [UCSC 1000 Genomes mirror](https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/), so no manual download is needed. You must have `bcftools` installed and on your `PATH`.

### 4. AlphaGenome API Key

Steps 02 and 03 call the AlphaGenome API and require an API key set as an environment variable:

```bash
export ALPHA_GENOME_API_KEY="your-key-here"
```

Request access through the [AlphaGenome page](https://deepmind.google/technologies/alphagenome/).

## Installation

### Prerequisites

- Python 3.11+
- [bcftools](http://www.htslib.org/) (for remote VCF slicing in step 01)

### Python Dependencies

```bash
pip install pandas numpy cyvcf2 pyliftover alphagenome pyarrow
```

| Package | Used In |
|---------|---------|
| `pandas`, `numpy` | All steps |
| `cyvcf2` | 01 (reading 1000G VCF) |
| `pyliftover` | 01B (hg19 -> hg38 coordinate liftover) |
| `alphagenome` | 02, 03 (AlphaGenome API client and variant scorers) |
| `pyarrow` | 02, 02B (parquet I/O) |

## Usage

Run the scripts in order. Each script reads the output of the previous step from the current working directory.

```bash
# 1. Extract GWAS variants and LD proxies (hg19)
python 01_prepare_snps_for_alphagenome.py

# 2. Liftover to GRCh38
python 01B_liftover_hg19_hg38.py

# 3. Run AlphaGenome predictions (all tissues x all modalities)
python 02_run_alphagenome_alltiss_allmodal.py

# 3b. (Optional) Re-collapse long output without re-calling API
python 02B_recollapse_from_long.py

# 4. Score variants with AlphaGenome recommended scorers
python 03_score_variant.py

# 5. Rank and prioritize SNPs
python 04_rank_prioritization.py

# 6. Collapse to eQTL/sQTL-like per-SNP-per-gene tables
python 05_collapse_eqtl_sqtl.py
```

Steps 02 and 03 include **checkpoint/resume support** -- if interrupted, re-running the script will skip already-completed variant-tissue pairs.

## Output Files

| File | Produced By | Description |
|------|-------------|-------------|
| `region_info.txt` | 01 | Locus metadata (lead SNP, region boundaries, parameters) |
| `region_significant_snps.csv` | 01 | All genome-wide significant SNPs in the window |
| `index_snps_independent_signals.csv` | 01 | Independent signal index SNPs after LD clumping |
| `ld_proxies_union_all_signals.csv` | 01 | Union of LD proxies across all signals (hg19) |
| `ld_proxies_union_all_signals_GRCh38.csv` | 01B | Proxies with GRCh38 coordinates |
| `alphagenome_allGTEx_allOutputs_LONG.parquet` | 02 | Track-level AlphaGenome predictions (long format) |
| `alphagenome_allGTEx_allOutputs_COLLAPSED.csv` | 02 / 02B | Per-variant/tissue/modality summary |
| `variant_scores_tidy.csv` | 03 | Tidy scores from recommended variant scorers |
| `snp_prioritization_by_scorer.csv` | 04 | Per-SNP rankings within each scorer |
| `snp_prioritization_overall.csv` | 04 | Consensus SNP ranking across all scorers |
| `snp_prioritization_by_tissue.csv` | 04 | Tissue-specific SNP rankings |
| `variant_scores_per_snp_gene_eqtl.csv` | 05 | Per-SNP-per-gene eQTL-like scores |
| `variant_scores_per_snp_gene_sqtl.csv` | 05 | Per-SNP-per-gene sQTL-like scores |

## References

1. Cheng, J., Novati, G., Pan, J. *et al.* Advancing regulatory variant effect prediction with AlphaGenome. *Nature* (2025). https://doi.org/10.1038/s41586-025-10014-0
2. Morris, J.A., Kemp, J.P., Youlten, S.E. *et al.* An atlas of genetic influences on osteoporosis in humans and mice. *Nature Genetics* **51**, 258-266 (2019). https://doi.org/10.1038/s41588-018-0302-x

## License

This project is for research and educational purposes. AlphaGenome predictions are generated under [Google DeepMind's terms](https://deepmind.google/technologies/alphagenome/).
