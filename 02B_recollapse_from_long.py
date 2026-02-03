#!/usr/bin/env python3
"""
Recollapse the long AlphaGenome parquet output into the collapsed CSV
with max-based summary columns used by step 03.

This avoids re-calling the AlphaGenome API. It requires a parquet engine
such as pyarrow or fastparquet.

Inputs:
  - alphagenome_allGTEx_allOutputs_LONG.parquet

Outputs:
  - alphagenome_allGTEx_allOutputs_COLLAPSED.csv
"""

import os
import sys
import numpy as np
import pandas as pd

IN_PARQUET = "alphagenome_allGTEx_allOutputs_LONG.parquet"
OUT_CSV = "alphagenome_allGTEx_allOutputs_COLLAPSED.csv"

USE_COLS = [
    "query_variant_id",
    "query_tissue",
    "query_gtex_tissue",
    "query_modality",
    "local_max_abs_delta",
    "local_mean_abs_delta",
]


def main():
    if not os.path.exists(IN_PARQUET):
        raise FileNotFoundError(f"Missing input parquet: {IN_PARQUET}")

    try:
        df = pd.read_parquet(IN_PARQUET, columns=USE_COLS)
    except Exception as e:
        msg = (
            "Failed to read parquet. You likely need a parquet engine.\n"
            "Try: pip install pyarrow   (recommended)\n"
            "Or:  pip install fastparquet\n"
            f"Original error: {e}"
        )
        print(msg, file=sys.stderr)
        raise

    # Normalize tissue column name
    if "query_tissue" not in df.columns and "query_gtex_tissue" in df.columns:
        df["query_tissue"] = df["query_gtex_tissue"]

    required = ["query_variant_id", "query_tissue", "query_modality",
                "local_max_abs_delta", "local_mean_abs_delta"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise RuntimeError(f"Missing required columns: {missing}. Found: {list(df.columns)}")

    gcols = ["query_variant_id", "query_tissue", "query_modality"]
    max_col = "local_max_abs_delta"
    mean_col = "local_mean_abs_delta"

    def _agg(g):
        return pd.Series({
            "n_rows": len(g),
            "n_max_nonmissing": int(g[max_col].notna().sum()),
            "n_mean_nonmissing": int(g[mean_col].notna().sum()),
            "best_local_max_abs_delta": float(g[max_col].max()) if g[max_col].notna().any() else np.nan,
            "best_local_mean_abs_delta": float(g[mean_col].max()) if g[mean_col].notna().any() else np.nan,
            "mean_max": float(g[max_col].mean()) if g[max_col].notna().any() else np.nan,
            "sd_max": float(g[max_col].std()) if g[max_col].notna().any() else np.nan,
        })

    collapsed = df.groupby(gcols, dropna=False).apply(_agg).reset_index()
    collapsed.to_csv(OUT_CSV, index=False)
    print(f"[ok] wrote {OUT_CSV} ({len(collapsed)} rows)")


if __name__ == "__main__":
    main()
