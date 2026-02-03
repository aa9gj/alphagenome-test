#!/usr/bin/env python3
"""
Rank-based prioritization of SNPs from AlphaGenome variant scores.

Takes the tidy scores from 04_score_variant.py and produces per-SNP
priority rankings aggregated across scorers and tissues.

Inputs:
  - variant_scores_tidy.csv  (from 04_score_variant.py)

Outputs:
  - snp_prioritization_by_scorer.csv   (per SNP x scorer summary)
  - snp_prioritization_overall.csv     (single ranking across all scorers)
  - snp_prioritization_by_tissue.csv   (per SNP x tissue summary)

Usage:
  python3.11 05_rank_prioritization.py
"""

import numpy as np
import pandas as pd

IN_CSV = "variant_scores_tidy.csv"

OUT_BY_SCORER = "snp_prioritization_by_scorer.csv"
OUT_OVERALL = "snp_prioritization_overall.csv"
OUT_BY_TISSUE = "snp_prioritization_by_tissue.csv"


def main():
    df = pd.read_csv(IN_CSV)
    print(f"[input] {len(df)} rows, {df['proxy_rsid'].nunique()} SNPs")

    # Determine which score column to use
    if "quantile_score" in df.columns and df["quantile_score"].notna().any():
        score_col = "quantile_score"
        print("[info] Using quantile_score (model-calibrated)")
    else:
        score_col = "raw_score"
        print("[info] Using raw_score")

    df["abs_score"] = df[score_col].abs()

    # ------------------------------------------------------------------
    # 1) Per SNP x scorer: max absolute score across all tracks/genes
    # ------------------------------------------------------------------
    by_scorer = (
        df.groupby(["proxy_rsid", "variant_scorer"])
        .agg(
            max_abs_score=("abs_score", "max"),
            mean_abs_score=("abs_score", "mean"),
            n_tracks=("abs_score", "count"),
        )
        .reset_index()
    )

    # Rank within each scorer (highest effect = rank 1)
    by_scorer["rank_within_scorer"] = (
        by_scorer.groupby("variant_scorer")["max_abs_score"]
        .rank(ascending=False, method="min")
        .astype(int)
    )

    by_scorer = by_scorer.sort_values(["variant_scorer", "rank_within_scorer"])
    by_scorer.to_csv(OUT_BY_SCORER, index=False)
    print(f"[done] Wrote {OUT_BY_SCORER} ({len(by_scorer)} rows)")

    # ------------------------------------------------------------------
    # 2) Overall prioritization: aggregate ranks across scorers
    # ------------------------------------------------------------------
    n_snps = by_scorer["proxy_rsid"].nunique()
    n_scorers = by_scorer["variant_scorer"].nunique()
    print(f"[info] {n_snps} SNPs x {n_scorers} scorers")

    overall = (
        by_scorer.groupby("proxy_rsid")
        .agg(
            mean_rank=("rank_within_scorer", "mean"),
            median_rank=("rank_within_scorer", "median"),
            best_rank=("rank_within_scorer", "min"),
            n_scorers_top3=("rank_within_scorer", lambda x: (x <= 3).sum()),
            n_scorers_scored=("rank_within_scorer", "count"),
            best_scorer=("max_abs_score", lambda x: by_scorer.loc[x.idxmax(), "variant_scorer"]),
            best_scorer_score=("max_abs_score", "max"),
        )
        .reset_index()
    )

    # Final rank by mean rank (lower mean rank = higher priority)
    overall["overall_rank"] = (
        overall["mean_rank"]
        .rank(ascending=True, method="min")
        .astype(int)
    )
    overall = overall.sort_values("overall_rank")
    overall.to_csv(OUT_OVERALL, index=False)
    print(f"[done] Wrote {OUT_OVERALL} ({len(overall)} rows)")

    # Print top SNPs
    print("\n=== TOP PRIORITIZED SNPs ===")
    for _, row in overall.head(10).iterrows():
        print(
            f"  Rank {int(row['overall_rank']):>2}: {row['proxy_rsid']:<15} "
            f"mean_rank={row['mean_rank']:.1f}  "
            f"best_rank={int(row['best_rank'])}  "
            f"top3_in={int(row['n_scorers_top3'])}/{int(row['n_scorers_scored'])} scorers  "
            f"best={row['best_scorer']}({row['best_scorer_score']:.4f})"
        )

    # ------------------------------------------------------------------
    # 3) Per SNP x tissue: for tissue-specific prioritization
    # ------------------------------------------------------------------
    tissue_col = None
    for candidate in ["gtex_tissue", "biosample_name"]:
        if candidate in df.columns and df[candidate].notna().any():
            tissue_col = candidate
            break

    if tissue_col is not None:
        by_tissue = (
            df.groupby(["proxy_rsid", tissue_col])
            .agg(
                max_abs_score=("abs_score", "max"),
                mean_abs_score=("abs_score", "mean"),
                n_tracks=("abs_score", "count"),
            )
            .reset_index()
        )
        by_tissue["rank_within_tissue"] = (
            by_tissue.groupby(tissue_col)["max_abs_score"]
            .rank(ascending=False, method="min")
            .astype(int)
        )
        by_tissue = by_tissue.sort_values([tissue_col, "rank_within_tissue"])
        by_tissue.to_csv(OUT_BY_TISSUE, index=False)
        print(f"\n[done] Wrote {OUT_BY_TISSUE} ({len(by_tissue)} rows)")
    else:
        print("\n[info] No tissue column found in scores; skipping tissue-level prioritization.")


if __name__ == "__main__":
    main()
