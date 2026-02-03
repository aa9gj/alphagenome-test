#!/usr/bin/env python3
import numpy as np
import pandas as pd
import warnings

IN_CSV = "alphagenome_allGTEx_allOutputs_COLLAPSED.csv"
OUT_CSV = "alphagenome_allGTEx_allOutputs_COLLAPSED_zscores.csv"
MIN_N_PER_MODALITY = 20   # set to 5 for small test runs; 20 is safer for real scoring

# Backward/forward compatibility:
#  - Older collapsed outputs used best_local_max_abs_delta / best_local_mean_abs_delta
#  - Newer collapsed outputs (from 02_run_alphagenome_alltiss_allmodal.py) use mean_max / sd_max
MAX_COL_CANDIDATES = [
    "best_local_max_abs_delta",
    "best_local_max",
    "mean_max",
]
MEAN_COL_CANDIDATES = [
    "best_local_mean_abs_delta",
    "best_local_mean",
]

def robust_z(x: pd.Series) -> pd.Series:
    x = x.astype(float)
    ok = x.notna()
    if ok.sum() < MIN_N_PER_MODALITY:
        return pd.Series(np.nan, index=x.index)
    vals = x[ok].values
    med = np.median(vals)
    mad = np.median(np.abs(vals - med))
    scale = 1.4826 * mad
    if scale == 0 or not np.isfinite(scale):
        return pd.Series(np.nan, index=x.index)
    z = (x - med) / scale
    return z

def standard_z(x: pd.Series) -> pd.Series:
    x = x.astype(float)
    ok = x.notna()
    if ok.sum() < MIN_N_PER_MODALITY:
        return pd.Series(np.nan, index=x.index)
    vals = x[ok].values
    mu = vals.mean()
    sd = vals.std(ddof=0)
    if sd == 0 or not np.isfinite(sd):
        return pd.Series(np.nan, index=x.index)
    return (x - mu) / sd

def pick_first_existing(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for col in candidates:
        if col in df.columns:
            return col
    return None

def main():
    df = pd.read_csv(IN_CSV)

    max_col = pick_first_existing(df, MAX_COL_CANDIDATES)
    mean_col = pick_first_existing(df, MEAN_COL_CANDIDATES)

    if max_col is None:
        raise RuntimeError(
            "Could not find a max-score column to z-score. "
            f"Tried: {MAX_COL_CANDIDATES}. "
            f"Found columns: {list(df.columns)}"
        )

    for col in [max_col, mean_col]:
        if col is not None:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    print(f"[info] z-scoring max values from column: {max_col}")
    if mean_col is None:
        print("[warn] No mean-score column found; z_mean_modality / rz_mean_modality will be NaN.")
    else:
        print(f"[info] z-scoring mean values from column: {mean_col}")

    # QC table (counts per modality)
    qc = (df.groupby("query_modality")
            .agg(
                n_rows=("query_modality","size"),
                n_max_nonmissing=(max_col, lambda s: s.notna().sum()),
                mean_max=(max_col,"mean"),
                sd_max=(max_col,"std"),
            )
            .reset_index()
         )
    if mean_col is not None:
        mean_counts = (df.groupby("query_modality")[mean_col]
                         .apply(lambda s: s.notna().sum())
                         .reset_index(name="n_mean_nonmissing"))
        qc = qc.merge(mean_counts, on="query_modality", how="left")
    else:
        qc["n_mean_nonmissing"] = 0
    qc.to_csv("qc_modality_counts.csv", index=False)

    # Silence numpy warnings from edge cases
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)

        df["z_max_modality"]  = df.groupby("query_modality")[max_col].transform(standard_z)
        df["rz_max_modality"]  = df.groupby("query_modality")[max_col].transform(robust_z)

        if mean_col is not None:
            df["z_mean_modality"] = df.groupby("query_modality")[mean_col].transform(standard_z)
            df["rz_mean_modality"] = df.groupby("query_modality")[mean_col].transform(robust_z)
        else:
            df["z_mean_modality"] = np.nan
            df["rz_mean_modality"] = np.nan

    # Choose one normalized score per row for downstream comparison
    df["score_modality"] = df["rz_max_modality"]

    df.to_csv(OUT_CSV, index=False)
    print(f"[ok] wrote {OUT_CSV}")
    print(f"[ok] wrote qc_modality_counts.csv (MIN_N_PER_MODALITY={MIN_N_PER_MODALITY})")

if __name__ == "__main__":
    main()
