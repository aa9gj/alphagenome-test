#!/usr/bin/env python3
"""
Score variants using AlphaGenome's built-in variant scorers.

Uses score_variants() with RECOMMENDED_VARIANT_SCORERS, which produce
properly calibrated per-gene and per-track effect scores that are
comparable across modalities without needing external normalization.

Inputs:
  - ld_proxies_union_all_signals_GRCh38.csv

Outputs:
  - variant_scores_tidy.csv   (long-format: one row per variant x scorer x track x gene)
  - variant_scores_checkpoint.csv  (resume support)

Usage:
  export ALPHA_GENOME_API_KEY="YOUR_KEY"
  python3.11 04_score_variant.py
"""

import os
import time
import pandas as pd

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers


# =========================
# USER SETTINGS
# =========================
INPUT_CSV = "ld_proxies_union_all_signals_GRCh38.csv"

SEQUENCE_LENGTH = dna_client.SEQUENCE_LENGTH_1MB

MAX_RETRIES = 6

# Debug cap (set to None for full run)
MAX_VARIANTS = None

# Checkpointing
CHECKPOINT_FILE = "variant_scores_checkpoint.csv"

# Output
OUT_CSV = "variant_scores_tidy.csv"


# =========================
# Client
# =========================
API_KEY = os.environ.get("ALPHA_GENOME_API_KEY", "")
if not API_KEY:
    raise RuntimeError("Missing ALPHA_GENOME_API_KEY environment variable")

ag = dna_client.create(api_key=API_KEY)


# =========================
# Scorers
# =========================
ALL_SCORERS = variant_scorers.RECOMMENDED_VARIANT_SCORERS
scorer_list = list(ALL_SCORERS.values())
print(f"[init] Using {len(scorer_list)} recommended scorers: {list(ALL_SCORERS.keys())}")


# =========================
# Helpers
# =========================
def normalize_variant_input(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
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

    df["chrom_ucsc"] = df["chrom_ucsc"].astype(str)
    df.loc[~df["chrom_ucsc"].str.startswith("chr"), "chrom_ucsc"] = "chr" + df["chrom_ucsc"]
    df["proxy_rsid"] = df["proxy_rsid"].astype(str)
    df["pos"] = df["pos"].astype(int)
    df["ref"] = df["ref"].astype(str)
    df["alt"] = df["alt"].astype(str)

    return df[required].copy()


def score_with_retries(interval, variant, scorers):
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            return ag.score_variant(
                interval=interval,
                variant=variant,
                variant_scorers=scorers,
            )
        except Exception as e:
            if attempt == MAX_RETRIES:
                raise
            wait = min(2**attempt, 30)
            print(f"[warn] score_variant failed attempt {attempt}/{MAX_RETRIES}: {e} | sleeping {wait}s")
            time.sleep(wait)


def load_checkpoint():
    if os.path.exists(CHECKPOINT_FILE):
        ck = pd.read_csv(CHECKPOINT_FILE)
        if "proxy_rsid" in ck.columns:
            return set(ck["proxy_rsid"].astype(str))
    return set()


def save_checkpoint(done_rsids):
    ck = pd.DataFrame({"proxy_rsid": sorted(done_rsids)})
    ck.to_csv(CHECKPOINT_FILE, index=False)


# =========================
# Main
# =========================
def main():
    raw = pd.read_csv(INPUT_CSV)
    cand = normalize_variant_input(raw)

    if MAX_VARIANTS is not None and len(cand) > MAX_VARIANTS:
        cand = cand.head(MAX_VARIANTS).copy()
        print(f"[debug] Capped variants to {MAX_VARIANTS}")

    print(f"[plan] Scoring {len(cand)} variants with {len(scorer_list)} scorers each")

    done = load_checkpoint()
    print(f"[resume] Already completed: {len(done)} variants")

    all_chunks = []
    t0 = time.time()

    for i, (_, v) in enumerate(cand.iterrows()):
        rsid = str(v["proxy_rsid"])
        if rsid in done:
            continue

        chrom = str(v["chrom_ucsc"])
        pos = int(v["pos"])
        ref = str(v["ref"])
        alt = str(v["alt"])

        variant = genome.Variant(
            chromosome=chrom,
            position=pos,
            reference_bases=ref,
            alternate_bases=alt,
            name=rsid,
        )
        interval = variant.reference_interval.resize(SEQUENCE_LENGTH)

        try:
            scores = score_with_retries(interval, variant, scorer_list)
            tidy = variant_scorers.tidy_scores(scores)
            if tidy is not None and not tidy.empty:
                tidy["proxy_rsid"] = rsid
                all_chunks.append(tidy)

            done.add(rsid)
            elapsed = (time.time() - t0) / 60
            print(f"[{i+1}/{len(cand)}] {rsid} done | elapsed={elapsed:.1f} min")

            if len(done) % 5 == 0:
                save_checkpoint(done)

        except Exception as e:
            print(f"[error] {rsid} failed: {e}")
            continue

    # Write outputs
    if all_chunks:
        df = pd.concat(all_chunks, ignore_index=True)
        df.to_csv(OUT_CSV, index=False)
        print(f"[done] Wrote {OUT_CSV} ({len(df)} rows)")
    else:
        print("[done] No scores produced; nothing to write.")

    save_checkpoint(done)
    print(f"[final] Checkpoint saved. total_done={len(done)}")


if __name__ == "__main__":
    main()
