#!/usr/bin/env python3
"""
Collapse AlphaGenome tidy scores to per-SNP-per-gene tables for:
  - eQTL-like signals (RNA-Seq gene-level scorers)
  - sQTL-like signals (splicing-related scorers/outputs)

Inputs:
  - variant_scores_tidy.csv

Outputs:
  - variant_scores_per_snp_gene_eqtl.csv
  - variant_scores_per_snp_gene_sqtl.csv
"""

import pandas as pd

IN_CSV = "variant_scores_tidy.csv"
OUT_EQTL = "variant_scores_per_snp_gene_eqtl.csv"
OUT_SQTL = "variant_scores_per_snp_gene_sqtl.csv"

USE_COLS = [
    "proxy_rsid",
    "variant_id",
    "scored_interval",
    "gene_id",
    "gene_name",
    "gene_type",
    "gene_strand",
    "variant_scorer",
    "output_type",
    "track_name",
    "gtex_tissue",
    "junction_Start",
    "junction_End",
    "raw_score",
    "quantile_score",
]


def collapse_per_snp_gene(df: pd.DataFrame) -> pd.DataFrame:
    # Use quantile_score when available, else raw_score
    score_col = "quantile_score" if df["quantile_score"].notna().any() else "raw_score"
    df = df.copy()
    df["abs_score"] = df[score_col].abs()

    keys = ["proxy_rsid", "gene_id", "gene_name"]

    agg = (
        df.groupby(keys)
        .agg(
            max_abs_score=("abs_score", "max"),
            mean_abs_score=("abs_score", "mean"),
            median_abs_score=("abs_score", "median"),
            n_rows=("abs_score", "size"),
            n_scorers=("variant_scorer", "nunique"),
            n_tissues=("gtex_tissue", "nunique"),
            n_tracks=("track_name", "nunique"),
            gene_type=("gene_type", "first"),
            gene_strand=("gene_strand", "first"),
        )
        .reset_index()
    )

    # Best row per SNP x gene (max abs score)
    idx = df.groupby(keys)["abs_score"].idxmax()
    best = df.loc[
        idx,
        [
            "proxy_rsid",
            "gene_id",
            "gene_name",
            score_col,
            "variant_scorer",
            "output_type",
            "track_name",
            "gtex_tissue",
            "junction_Start",
            "junction_End",
            "variant_id",
            "scored_interval",
        ],
    ].copy()

    best = best.rename(
        columns={
            score_col: "best_score",
            "variant_scorer": "best_scorer",
            "output_type": "best_output_type",
            "track_name": "best_track_name",
            "gtex_tissue": "best_gtex_tissue",
            "junction_Start": "best_junction_start",
            "junction_End": "best_junction_end",
            "variant_id": "best_variant_id",
            "scored_interval": "best_scored_interval",
        }
    )

    out = agg.merge(best, on=keys, how="left")
    out["rank_within_snp"] = (
        out.groupby("proxy_rsid")["max_abs_score"].rank(ascending=False, method="min").astype(int)
    )
    return out.sort_values(["proxy_rsid", "rank_within_snp"])


def main():
    df = pd.read_csv(IN_CSV, usecols=USE_COLS, low_memory=False)

    # Keep only gene-annotated rows (required for per-SNP-per-gene tables)
    df = df[df["gene_name"].notna() & df["gene_id"].notna()].copy()

    # eQTL-like: RNA_SEQ gene-level scorers (exclude PolyadenylationScorer)
    eqtl = df[
        (df["output_type"] == "RNA_SEQ")
        & (df["variant_scorer"].str.contains("GeneMask", na=False))
    ]
    eqtl_out = collapse_per_snp_gene(eqtl)
    eqtl_out.to_csv(OUT_EQTL, index=False)

    # sQTL-like: splicing-related outputs/scorers
    sqtl = df[
        df["output_type"].isin(["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"])
        | df["variant_scorer"].str.contains("Splicing", na=False)
        | df["variant_scorer"].str.contains("SpliceJunction", na=False)
    ]
    sqtl_out = collapse_per_snp_gene(sqtl)
    sqtl_out.to_csv(OUT_SQTL, index=False)

    print(f"[ok] wrote {OUT_EQTL} (rows={len(eqtl_out)})")
    print(f"[ok] wrote {OUT_SQTL} (rows={len(sqtl_out)})")


if __name__ == "__main__":
    main()
