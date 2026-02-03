#!/usr/bin/env python3
"""
Run AlphaGenome for a set of GRCh38 variants across:
  - ALL available GTEx tissues (as exposed by ag.output_metadata().rna_seq where data_source == 'gtex')
  - ALL available OutputTypes (modalities) in your installed AlphaGenome client

Inputs:
  - ld_proxies_union_all_signals_GRCh38.csv  (must include rsID + chr/pos + ref/alt)

Outputs:
  - alphagenome_allGTEx_allOutputs_LONG.parquet  (track-level rows; can be big)
  - alphagenome_allGTEx_allOutputs_COLLAPSED.csv (per variant x tissue x modality summary)
  - progress_checkpoint.csv (resume support)

Usage:
  export ALPHA_GENOME_API_KEY="YOUR_KEY"
  python3.11 02_run_alphagenome_allGTEx_allOutputs.py
"""

import os
import time
import math
import numpy as np
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client


# =========================
# USER SETTINGS
# =========================
INPUT_CSV = "ld_proxies_union_all_signals_GRCh38.csv"

# Window AlphaGenome expects (their model uses fixed windows; 2**20 is used in your working scripts)
WINDOW_BP = 2**20
HALF_WINDOW = WINDOW_BP // 2

# Local summarization around variant (bp)
LOCAL_BP = 2000

# Throttling / robustness
SLEEP_BETWEEN_CALLS = 0.0
MAX_RETRIES = 6

# Debug caps (set to None for full run)
MAX_VARIANTS = None          # e.g. 50 for debug
MAX_TISSUES = None           # e.g. 10 for debug

DEBUG_PRINT_SHAPES_ONCE = False  # set True to print values_shape/interval for first rsid+tissue per modality


# Checkpointing
CHECKPOINT_EVERY_CALLS = 25
CHECKPOINT_FILE = "progress_checkpoint.csv"

# Outputs
OUT_LONG_PARQUET = "alphagenome_allGTEx_allOutputs_LONG.parquet"
OUT_COLLAPSED_CSV = "alphagenome_allGTEx_allOutputs_COLLAPSED.csv"


# =========================
# Client
# =========================
API_KEY = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not API_KEY:
    raise RuntimeError("Missing ALPHA_GENOME_API_KEY environment variable")

ag = dna_client.create(api_key=API_KEY)


# =========================
# Output types / modalities
# =========================
# Pull all OutputTypes available in this SDK build:
REQUESTED = list(dna_client.OutputType)
print(f"[init] OutputTypes available: {[ot.name for ot in REQUESTED]}")


# =========================
# Helpers
# =========================
def interval_centered(chrom_ucsc: str, pos_1based: int) -> genome.Interval:
    # Clamp start to >= 1 (Interval is 1-based-like semantics downstream)
    start = max(1, int(pos_1based) - HALF_WINDOW)
    end = int(pos_1based) + HALF_WINDOW
    return genome.Interval(chromosome=chrom_ucsc, start=start, end=end)


def normalize_variant_input(df: pd.DataFrame) -> pd.DataFrame:
    """
    Best-effort normalization for expected columns.
    Requires at least:
      - proxy_rsid
      - chrom_ucsc (e.g. 'chr1')
      - pos (1-based)
      - ref
      - alt
    """
    df = df.copy()

    # Common aliases
    colmap = {
        "rsid": "proxy_rsid",
        "rsID": "proxy_rsid",
        "chrom": "chrom_ucsc",
        "chr": "chrom_ucsc",
        "pos_1based": "pos",
        "position": "pos",
        "pos_hg38": "pos",
        "ref_allele": "ref",
        "alt_allele": "alt",
        "ref_hg38": "ref",
        "alt_hg38": "alt",
    }
    for a, b in colmap.items():
        if a in df.columns and b not in df.columns:
            df[b] = df[a]

    required = ["proxy_rsid", "chrom_ucsc", "pos", "ref", "alt"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"INPUT_CSV missing required columns: {missing}. Found: {list(df.columns)}")

    # Normalize chr prefix
    df["chrom_ucsc"] = df["chrom_ucsc"].astype(str)
    df.loc[~df["chrom_ucsc"].str.startswith("chr"), "chrom_ucsc"] = "chr" + df["chrom_ucsc"]

    # Clean types
    df["proxy_rsid"] = df["proxy_rsid"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    return df[required].copy()


def get_all_gtex_tissues_and_ontologies():
    """
    Query AlphaGenome output metadata and return:
      tissues: list[str]
      tissue_to_ontology: dict[str, str]
    using RNA-seq metadata rows where data_source == 'gtex'.
    """
    md = ag.output_metadata()
    rna = md.rna_seq
    if not isinstance(rna, pd.DataFrame):
        rna = pd.DataFrame(rna)

    # Try common column names
    # Expect something like: data_source, biosample_name, ontology_term_id
    # Different builds may use slightly different labels, so try robust fallbacks.
    def pick_col(candidates):
        for c in candidates:
            if c in rna.columns:
                return c
        return None

    ds_col = pick_col(["data_source", "dataSource", "source"])
    tissue_col = pick_col(["biosample_name", "biosample", "tissue", "context", "name"])
    ont_col = pick_col(["ontology_curie", "ontology_term_id", "ontologyTermId", "ontology", "ontology_id", "ont_id"])

    if ds_col is None or tissue_col is None or ont_col is None:
        raise RuntimeError(
            f"Could not find expected metadata columns in md.rna_seq. "
            f"Columns: {list(rna.columns)}"
        )

    rna_gtex = rna[rna[ds_col].astype(str).str.lower() == "gtex"].copy()
    rna_gtex = rna_gtex.dropna(subset=[tissue_col, ont_col])

    tissue_to_ont = dict(zip(rna_gtex[tissue_col].astype(str), rna_gtex[ont_col].astype(str)))
    tissues = sorted(tissue_to_ont.keys())

    return tissues, tissue_to_ont


def predict_with_retries(interval: genome.Interval, var: genome.Variant, ontologies, outputs):
    """
    Predict with retry + exponential-ish backoff.
    """
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            out = ag.predict_variant(
                interval=interval,
                variant=var,
                ontology_terms=ontologies,
                requested_outputs=outputs,
            )
            return out
        except Exception as e:
            if attempt == MAX_RETRIES:
                raise
            wait = min(2**attempt, 30)
            print(f"[warn] predict failed attempt {attempt}/{MAX_RETRIES}: {e} | sleeping {wait}s")
            time.sleep(wait)


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        ck = pd.read_csv(CHECKPOINT_FILE)
        # Expect columns: proxy_rsid, tissue
        if "proxy_rsid" in ck.columns and "tissue" in ck.columns:
            return set(zip(ck["proxy_rsid"].astype(str), ck["tissue"].astype(str)))
    return set()


def save_checkpoint(done_pairs):
    ck = pd.DataFrame(list(done_pairs), columns=["proxy_rsid", "tissue"])
    ck.to_csv(CHECKPOINT_FILE, index=False)


def extract_track(tref, talt):
    """
    Normalize an output track (if present) into:
      ref_vals: np.ndarray [L x K]
      alt_vals: np.ndarray [L x K]
      md_df: pandas DataFrame with K rows (track metadata)
      interval: genome.Interval
    Returns None if missing / non-track output.
    """
    if tref is None or talt is None:
        return None
    if not hasattr(tref, "values") or not hasattr(tref, "metadata") or not hasattr(tref, "interval"):
        return None
    md = tref.metadata
    if md is None:
        return None
    if not isinstance(md, pd.DataFrame):
        # Some versions may store metadata differently; try best-effort conversion
        try:
            md = pd.DataFrame(md)
        except Exception:
            return None
    return (
        np.asarray(tref.values),
        np.asarray(talt.values),
        md.copy(),
        tref.interval,
    )


def summarize_delta(ref_vals, alt_vals, interval, var_pos_1based: int):
    """
    Robust per-track summary near the variant, handling base-level and binned tracks.

    Returns:
      local_max_abs_delta: np.ndarray shape (K,)
      local_mean_abs_delta: np.ndarray shape (K,)

    Supports:
      - 1D arrays [L]  (single track)
      - 2D arrays [L, K] (K tracks)
    For other shapes (e.g. CONTACT_MAPS often 3D/4D), returns NaNs.
    """
    ref_vals = np.asarray(ref_vals)
    alt_vals = np.asarray(alt_vals)

    # Normalize to 2D [L, K]
    if ref_vals.ndim == 1:
        ref_vals = ref_vals[:, None]
        alt_vals = alt_vals[:, None]
    elif ref_vals.ndim != 2:
        # Unsupported (e.g., contact maps): return NaNs for best-effort downstream handling
        K = int(ref_vals.shape[-1]) if ref_vals.ndim >= 1 else 1
        return (np.full((K,), np.nan, dtype=float), np.full((K,), np.nan, dtype=float))

    delta = alt_vals - ref_vals
    L, K = delta.shape

    # Map genomic position (bp) -> track index (bin)
    interval_len_bp = int(interval.end) - int(interval.start)
    if interval_len_bp <= 0 or L <= 0:
        return (np.full((K,), np.nan, dtype=float), np.full((K,), np.nan, dtype=float))

    bp_per_bin = interval_len_bp / float(L)  # ~1 for base-resolution tracks; >1 for binned tracks
    i_bin = int((int(var_pos_1based) - int(interval.start)) / bp_per_bin)

    # Convert LOCAL_BP (bp) -> number of bins to include around the variant
    local_bins = max(1, int(math.ceil(LOCAL_BP / bp_per_bin)))

    sl = slice(max(i_bin - local_bins, 0), min(i_bin + local_bins + 1, L))
    if sl.stop <= sl.start:
        return (np.full((K,), np.nan, dtype=float), np.full((K,), np.nan, dtype=float))

    chunk = np.abs(delta[sl, :])
    if chunk.size == 0:
        return (np.full((K,), np.nan, dtype=float), np.full((K,), np.nan, dtype=float))

    local_max = np.max(chunk, axis=0)
    local_mean = np.mean(chunk, axis=0)
    return local_max, local_mean


def collapse_long(long_df: pd.DataFrame) -> pd.DataFrame:
    """
    Collapse per-track long rows to per-variant x tissue x modality summary.
    Produces:
      - n_rows
      - n_max_nonmissing
      - n_mean_nonmissing
      - best_local_max_abs_delta   (max across tracks)
      - best_local_mean_abs_delta  (max across tracks)
      - mean_max
      - sd_max
    """
    if long_df.empty:
        return pd.DataFrame()

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

    return long_df.groupby(gcols, dropna=False).apply(_agg).reset_index()


def main():
    # 1) Load variants
    raw = pd.read_csv(INPUT_CSV)
    cand = normalize_variant_input(raw)

    if MAX_VARIANTS is not None and len(cand) > MAX_VARIANTS:
        cand = cand.head(MAX_VARIANTS).copy()
        print(f"[debug] Capped variants to {MAX_VARIANTS}")

    # 2) Enumerate GTEx tissues
    tissues, tissue_to_ont = get_all_gtex_tissues_and_ontologies()
    print(f"[metadata] GTEx tissues available in this AlphaGenome build: {len(tissues)}")

    if MAX_TISSUES is not None and len(tissues) > MAX_TISSUES:
        tissues = tissues[:MAX_TISSUES]
        print(f"[debug] Capped GTEx tissues to {len(tissues)}")

    planned_calls = len(cand) * len(tissues)
    print(f"[plan] API calls (variants x tissues): {len(cand)} x {len(tissues)} = {planned_calls}")

    # 3) Resume support
    done = load_checkpoint()
    print(f"[resume] Completed (variant,tissue) pairs loaded: {len(done)}")

    # 4) Run
    long_chunks = []
    collapsed_rows = []
    done_to_append = []
    calls = 0
    t0 = time.time()

    debug_shapes_printed = False

    for _, v in cand.iterrows():
        rsid = str(v["proxy_rsid"])
        chrom = str(v["chrom_ucsc"])
        pos = int(v["pos"])
        ref = str(v["ref"])
        alt = str(v["alt"])

        interval = interval_centered(chrom, pos)
        var = genome.Variant(chromosome=chrom, position=pos, reference_bases=ref, alternate_bases=alt)

        for tissue in tissues:
            key = (rsid, tissue)
            if key in done:
                continue

            ont = tissue_to_ont[tissue]

            try:
                out = predict_with_retries(interval, var, [ont], REQUESTED)

                for ot in REQUESTED:
                    attr = ot.name.lower()

                    tref = getattr(out.reference, attr, None)
                    talt = getattr(out.alternate, attr, None)

                    parsed = extract_track(tref, talt)
                    if parsed is None:
                        continue

                    ref_vals, alt_vals, md_df, t_interval = parsed

                    if DEBUG_PRINT_SHAPES_ONCE and (not debug_shapes_printed):
                        try:
                            print(
                                f"[shape] rsid={rsid} tissue={tissue} modality={ot.name} "
                                f"values_shape={np.asarray(ref_vals).shape} interval=({t_interval.start},{t_interval.end})"
                            )
                        except Exception as e:
                            print(f"[shape] failed to print shapes for modality={ot.name}: {e}")
                        debug_shapes_printed = True

                    local_max, local_mean = summarize_delta(ref_vals, alt_vals, t_interval, pos)

                    md_df = md_df.copy()
                    md_df["local_max_abs_delta"] = local_max
                    md_df["local_mean_abs_delta"] = local_mean

                    # query context (prefix to avoid collisions)
                    md_df["query_variant_id"] = rsid
                    md_df["query_chrom"] = chrom
                    md_df["query_pos_1based"] = pos
                    md_df["query_ref"] = ref
                    md_df["query_alt"] = alt
                    md_df["query_tissue"] = tissue
                    md_df["query_ontology_term"] = ont
                    md_df["query_modality"] = ot.name

                    long_chunks.append(md_df)

                # Mark completed pair
                done.add(key)
                done_to_append.append(key)
                calls += 1

                if SLEEP_BETWEEN_CALLS > 0:
                    time.sleep(SLEEP_BETWEEN_CALLS)

                if calls % CHECKPOINT_EVERY_CALLS == 0:
                    save_checkpoint(done)
                    print(f"[checkpoint] saved after {calls} calls. elapsed={(time.time()-t0)/60:.1f} min")

            except Exception as e:
                print(f"[error] rsid={rsid} tissue={tissue} failed: {e}")
                # Do not mark as done; will retry on next run.
                continue

    # 5) Write outputs
    if long_chunks:
        long_df = pd.concat(long_chunks, ignore_index=True)
    else:
        long_df = pd.DataFrame()

    if not long_df.empty:
        long_df.to_parquet(OUT_LONG_PARQUET, index=False)
        collapsed = collapse_long(long_df)
        collapsed.to_csv(OUT_COLLAPSED_CSV, index=False)
        print(f"[done] wrote {OUT_LONG_PARQUET} and {OUT_COLLAPSED_CSV}")
    else:
        print("[done] No long rows produced; nothing to write.")

    # Final checkpoint
    save_checkpoint(done)
    print(f"[final] checkpoint saved to {CHECKPOINT_FILE}. total_done={len(done)}")


if __name__ == "__main__":
    main()
