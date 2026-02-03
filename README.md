# Interrogating the FHL3 Coronary Artery Disease GWAS Locus with AlphaGenome

This repository demonstrates using [AlphaGenome](https://deepmind.google/blog/alphagenome-ai-for-better-understanding-the-genome/) (Google DeepMind) to predict the regulatory impact of non-coding variants at the **FHL3** (*Four-and-a-Half LIM Domains 3*) GWAS locus on chromosome 1p34.3, a locus associated with coronary artery disease (CAD).

## Background

### The FHL3 Locus and CAD

Genome-wide association studies have identified variants near *FHL3* as significantly associated with coronary artery disease. Fine-mapping and single-cell epigenomic analyses in human aortic smooth muscle cells (HASMCs) have implicated FHL3 as a causal gene at this locus through allele-specific enhancer activity ([Wirka et al., 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC8260472/)). Mendelian randomization further supports a causal link between FHL3 expression and CAD risk ([Hao et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8284374/)).

Key regulatory variant **rs61776719** at this locus is predicted to affect the expression of FHL3 along with neighboring genes MANEAL, INPP5B, and C1orf122.

### AlphaGenome

[AlphaGenome](https://www.nature.com/articles/s41586-025-10014-0) is a deep learning model that takes up to 1 Mb of DNA sequence as input and predicts thousands of functional genomic tracks at single-base-pair resolution, including:

- **Gene expression** (RNA-seq)
- **Chromatin accessibility** (ATAC-seq / DNase-seq)
- **Histone modifications** (ChIP-seq)
- **Transcription factor binding**
- **Chromatin contact maps** (Hi-C / Micro-C)
- **Splice site usage and junction strength**

By comparing predictions for reference vs. alternate alleles, AlphaGenome enables variant effect prediction (VEP) across all these modalities, offering mechanistic insight into how non-coding GWAS variants perturb gene regulation.

## Analysis Overview

The analysis is organized across five notebooks:

| Notebook | Description |
|----------|-------------|
| `01-gwas-variant-retrieval.ipynb` | Retrieve GWAS summary statistics and fine-mapped variants at the FHL3 locus from the GWAS Catalog and published fine-mapping studies |
| `02-alphagenome-predictions.ipynb` | Run AlphaGenome variant effect predictions for lead and credible set SNPs across relevant cell types and tissues |
| `03-variant-effect-analysis.ipynb` | Quantify and rank variant effects on chromatin accessibility, histone marks, TF binding, and gene expression tracks |
| `04-visualization.ipynb` | Visualize reference vs. alternate allele prediction tracks, effect-size heatmaps, and locus-level regulatory landscapes |
| `05-interpretation.ipynb` | Integrate AlphaGenome predictions with external annotations (eQTL data, epigenomic maps, single-cell data) to build a mechanistic model of variant-to-gene-to-disease |

## Repository Structure

```
.
├── README.md
├── .gitignore
├── notebooks/
│   ├── 01-gwas-variant-retrieval.ipynb
│   ├── 02-alphagenome-predictions.ipynb
│   ├── 03-variant-effect-analysis.ipynb
│   ├── 04-visualization.ipynb
│   └── 05-interpretation.ipynb
├── data/
│   ├── gwas/                  # GWAS summary stats and variant lists
│   └── annotations/           # External annotations (eQTL, epigenomic)
├── results/
│   ├── predictions/           # AlphaGenome prediction outputs
│   ├── figures/               # Generated plots and visualizations
│   └── tables/                # Summary tables and rankings
└── environment.yml            # Conda environment specification
```

## Getting Started

### Prerequisites

- Python 3.10+
- An [AlphaGenome API key](https://deepmind.google/science/alphagenome/) (free for non-commercial use)

### Installation

```bash
git clone https://github.com/aa9gj/alphagenome-test.git
cd alphagenome-test

# Create and activate environment
conda env create -f environment.yml
conda activate alphagenome-fhl3

# Or install with pip
pip install alphagenome jupyter pandas matplotlib seaborn
```

### Configuration

Set your AlphaGenome API key as an environment variable:

```bash
export ALPHAGENOME_API_KEY="your-api-key-here"
```

### Running the Analysis

Execute the notebooks in order:

```bash
jupyter notebook notebooks/
```

## Key Results

AlphaGenome variant effect predictions at the FHL3 locus reveal:

- **Disrupted enhancer activity**: Lead variants alter predicted chromatin accessibility and H3K27ac signal in smooth muscle-relevant cell types, consistent with allele-specific enhancer activity observed in HASMCs
- **Directional gene expression effects**: Variant effects on predicted RNA-seq tracks are concordant with eQTL effect directions for FHL3 in arterial tissues (GTEx)
- **Tissue-specific regulation**: The strongest predicted effects are concentrated in vascular and cardiac cell types, aligning with the disease-relevant biology of CAD

## References

1. Cheng, J., Novati, G., Pan, J. *et al.* Advancing regulatory variant effect prediction with AlphaGenome. *Nature* (2025). https://doi.org/10.1038/s41586-025-10014-0
2. Wirka, R.C., Wagh, D., Paik, D.T. *et al.* Single-Cell Epigenomics and Functional Fine-Mapping of Atherosclerosis GWAS Loci. *Circ Res* (2022). https://doi.org/10.1161/CIRCRESAHA.121.318971
3. Hao, K., Ermel, R., Engert, J.C. *et al.* System Genetics Including Causal Inference Identify Immune Targets for Coronary Artery Disease and the Lifespan. *Circ Genom Precis Med* (2021). https://doi.org/10.1161/CIRCGEN.120.003196

## License

This project is for research and educational purposes. AlphaGenome predictions are generated under [Google DeepMind's API terms](https://deepmind.google/science/alphagenome/) (non-commercial use).
