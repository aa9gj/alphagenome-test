#!/usr/bin/env python3
import os
import pandas as pd
from pyliftover import LiftOver

# -----------------------------
# INPUTS
# -----------------------------
PROXIES_CSV = "ld_proxies_union_all_signals.csv"  # produced by your 01_prepare script
LIFTOVER_DIR = "/Users/aa9gj/Documents/fhl3_alphagenome_res/liftOver"
CHAIN_19_TO_38 = os.path.join(LIFTOVER_DIR, "hg19ToHg38.over.chain.gz")

OUT_CSV = "ld_proxies_union_all_signals_GRCh38.csv"

BASE_COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

def revcomp_allele(a: str) -> str:
    a = str(a).upper()
    if a not in BASE_COMP:
        return a
    return BASE_COMP[a]

def die(msg: str):
    raise SystemExit(f"[ERROR] {msg}")

def main():
    if not os.path.exists(PROXIES_CSV):
        die(f"Missing input: {PROXIES_CSV}")

    if not os.path.exists(CHAIN_19_TO_38):
        die(
            f"Missing chain file: {CHAIN_19_TO_38}\n"
            f"Expected it in: {LIFTOVER_DIR}\n"
            f"Make sure the filename is exactly hg19ToHg38.over.chain.gz"
        )

    # Load proxies
    df = pd.read_csv(PROXIES_CSV)
    required = ["index_rsid","proxy_rsid","chr","pos_hg19","ref_1000g","alt_1000g","maf_eur","r2_eur","distance_to_index_bp"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        die(f"Input CSV missing columns: {missing}\nFound: {df.columns.tolist()}")

    df["chr"] = df["chr"].astype(str).str.replace("^chr", "", regex=True)
    df["pos_hg19"] = df["pos_hg19"].astype(int)

    lo = LiftOver(CHAIN_19_TO_38)

    chr38_list = []
    pos38_list = []
    strand_list = []
    ref38_list = []
    alt38_list = []

    n_missed = 0
    for _, r in df.iterrows():
        chrom = f"chr{r['chr']}"
        pos1 = int(r["pos_hg19"])

        hits = lo.convert_coordinate(chrom, pos1 - 1)  # pyliftover uses 0-based
        if not hits:
            chr38_list.append(None)
            pos38_list.append(None)
            strand_list.append(None)
            ref38_list.append(None)
            alt38_list.append(None)
            n_missed += 1
            continue

        # Take best hit
        c2, p2_0, strand, _ = hits[0]
        chr38 = c2.replace("chr", "")
        pos38 = int(p2_0) + 1

        ref = str(r["ref_1000g"]).upper()
        alt = str(r["alt_1000g"]).upper()

        # If chain flips strand, alleles need reverse-complement
        if strand == "-":
            ref = revcomp_allele(ref)
            alt = revcomp_allele(alt)

        chr38_list.append(chr38)
        pos38_list.append(pos38)
        strand_list.append(strand)
        ref38_list.append(ref)
        alt38_list.append(alt)

    df["chr_hg38"] = chr38_list
    df["pos_hg38"] = pos38_list
    df["liftover_strand"] = strand_list
    df["ref_hg38"] = ref38_list
    df["alt_hg38"] = alt38_list

    out = df.dropna(subset=["chr_hg38","pos_hg38"]).copy()
    out.to_csv(OUT_CSV, index=False)

    print(f"[ok] wrote {OUT_CSV} with {len(out)} rows (missed {n_missed})")
    print("[note] If you see liftover_strand == '-', ref/alt were reverse-complemented.")

if __name__ == "__main__":
    main()
