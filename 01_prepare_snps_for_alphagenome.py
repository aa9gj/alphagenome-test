#!/usr/bin/env python3
"""
PI-grade Step 1: multi-signal aware LD proxies in a ±200kb window (GRCh37/hg19).

Fixes included:
- Normalize dotted sumstats cols: "P.NI" -> "P_NI"
- Handle both 1000G panel schemas:
    (sample, population)  OR  (sample, pop, super_pop, ...)
"""

import os
import math
import subprocess
import pandas as pd
import numpy as np
from cyvcf2 import VCF

# -----------------------
# USER SETTINGS
# -----------------------
SUMSTATS = "/Users/aa9gj/Downloads/Morrisetal2018.NatGen.SumStats.tar_0/Biobank2-British-Bmd-As-C-Gwas-SumStats.txt"
LEAD_RSID = "rs28428561"
WINDOW_BP = 200_000

GW_SIG = 5e-8
CLUMP_R2 = 0.10
PROXY_R2 = 0.80

MAF_MIN = 0.01
BIALLELIC_SNP_ONLY = True

# Prefer super_pop EUR if available; otherwise use these EUR pops.
EUR_POPS = {"CEU", "TSI", "FIN", "GBR", "IBS"}

PANEL_URL = "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"
PANEL_FILE = "integrated_call_samples_v3.20130502.ALL.panel"

UCSC_BASE = "https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3"
VCF_TEMPLATE = "ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

REGION_INFO = "region_info.txt"
REGION_VCF = "1000G_region.vcf.gz"
REGION_VCF_TBI = "1000G_region.vcf.gz.tbi"
OUT_SIG = "region_significant_snps.csv"
OUT_INDEX = "index_snps_independent_signals.csv"
OUT_PROXIES = "ld_proxies_union_all_signals.csv"


def run(cmd: list[str]):
    print("[cmd]", " ".join(cmd))
    subprocess.check_call(cmd)

def ensure_panel():
    if os.path.exists(PANEL_FILE):
        return
    run(["curl", "-L", "-o", PANEL_FILE, PANEL_URL])

def eur_samples_from_panel() -> set[str]:
    """
    Supports both common 1000G panel schemas:
      - columns: sample, population
      - columns: sample, pop, super_pop, gender, ...
    Preference:
      - If super_pop exists: take super_pop == 'EUR'
      - Else: use pop/population in EUR_POPS
    """
    panel = pd.read_csv(PANEL_FILE, sep="\t", dtype=str)

    if "sample" not in panel.columns:
        raise RuntimeError(f"Panel missing 'sample' column. columns={panel.columns.tolist()}")

    # Prefer super_pop if present
    if "super_pop" in panel.columns:
        eur = panel[panel["super_pop"] == "EUR"]["sample"].astype(str).tolist()
        return set(eur)

    # Otherwise fall back to pop/population
    pop_col = None
    if "population" in panel.columns:
        pop_col = "population"
    elif "pop" in panel.columns:
        pop_col = "pop"

    if pop_col is None:
        raise RuntimeError(f"Panel missing expected pop column. columns={panel.columns.tolist()}")

    eur = panel[panel[pop_col].isin(sorted(EUR_POPS))]["sample"].astype(str).tolist()
    return set(eur)

def bcftools_extract_region(chrom: str, start: int, end: int):
    vcf_url = f"{UCSC_BASE}/{VCF_TEMPLATE.format(chrom=chrom)}"
    region = f"{chrom}:{start}-{end}"
    if os.path.exists(REGION_VCF) and os.path.exists(REGION_VCF_TBI):
        return
    run(["bcftools", "view", "-r", region, "-Oz", "-o", REGION_VCF, vcf_url])
    run(["bcftools", "index", "-t", REGION_VCF])

def dosage_array(genotypes) -> np.ndarray:
    d = np.full((len(genotypes),), np.nan, dtype=float)
    for i, g in enumerate(genotypes):
        a1, a2 = g[0], g[1]
        if a1 < 0 or a2 < 0:
            continue
        d[i] = a1 + a2
    return d

def maf_from_dosage(d: np.ndarray) -> float:
    mask = ~np.isnan(d)
    if mask.sum() == 0:
        return np.nan
    p_alt = np.nanmean(d) / 2.0
    return float(min(p_alt, 1.0 - p_alt))

def r2(d1: np.ndarray, d2: np.ndarray) -> float:
    mask = ~np.isnan(d1) & ~np.isnan(d2)
    if mask.sum() < 20:
        return np.nan
    x = d1[mask].astype(float)
    y = d2[mask].astype(float)
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float((x * x).sum() * (y * y).sum()))
    if denom == 0.0:
        return np.nan
    r = float((x * y).sum() / denom)
    return r * r

def is_biallelic_snp(var) -> bool:
    if not BIALLELIC_SNP_ONLY:
        return True
    if len(var.ALT) != 1:
        return False
    ref = var.REF.upper()
    alt = var.ALT[0].upper()
    return (len(ref) == 1 and len(alt) == 1 and ref in "ACGT" and alt in "ACGT")

def allele_set_key(pos_1based: int, a1: str, a2: str):
    return (int(pos_1based), tuple(sorted([str(a1).upper(), str(a2).upper()])))

def normalize_sumstats_cols(df: pd.DataFrame) -> pd.DataFrame:
    rename_map = {}
    if "P.NI" in df.columns:
        rename_map["P.NI"] = "P_NI"
    if "P.I" in df.columns:
        rename_map["P.I"] = "P_I"
    if rename_map:
        df = df.rename(columns=rename_map)
    return df

def load_lead_row_and_window_signals(sumstats_path: str, lead_rsid: str, window_bp: int):
    usecols = ["SNPID", "RSID", "CHR", "BP", "EA", "NEA", "EAF", "INFO",
               "BETA", "SE", "P", "P.I", "P.NI", "N"]

    lead_row = None
    for chunk in pd.read_csv(sumstats_path, sep="\t", usecols=lambda c: c in usecols,
                             chunksize=500_000, dtype=str):
        chunk = normalize_sumstats_cols(chunk)
        hit = chunk[chunk["RSID"] == lead_rsid]
        if not hit.empty:
            lead_row = hit.iloc[0].to_dict()
            break
    if lead_row is None:
        raise RuntimeError(f"Lead RSID {lead_rsid} not found in sumstats.")

    chrom = str(lead_row["CHR"])
    lead_bp = int(float(lead_row["BP"]))
    start = max(1, lead_bp - window_bp)
    end = lead_bp + window_bp

    sig_rows = []
    for chunk in pd.read_csv(sumstats_path, sep="\t", usecols=lambda c: c in usecols,
                             chunksize=500_000, dtype=str):
        chunk = normalize_sumstats_cols(chunk)
        chunk = chunk[chunk["CHR"] == chrom]
        if chunk.empty:
            continue
        bp = pd.to_numeric(chunk["BP"], errors="coerce")
        chunk = chunk[(bp >= start) & (bp <= end)]
        if chunk.empty:
            continue
        pni = pd.to_numeric(chunk["P_NI"], errors="coerce")
        chunk = chunk[pni < GW_SIG]
        if chunk.empty:
            continue
        sig_rows.append(chunk)

    sig_df = pd.concat(sig_rows, ignore_index=True) if sig_rows else pd.DataFrame(columns=normalize_sumstats_cols(pd.DataFrame(columns=usecols)).columns)
    if not sig_df.empty:
        sig_df["BP"] = pd.to_numeric(sig_df["BP"], errors="coerce").astype("Int64")
        sig_df["P_NI"] = pd.to_numeric(sig_df["P_NI"], errors="coerce")
    return lead_row, sig_df, (chrom, start, end)

def main():
    lead, sig_df, (chrom, start, end) = load_lead_row_and_window_signals(SUMSTATS, LEAD_RSID, WINDOW_BP)
    lead_bp = int(float(lead["BP"]))
    lead_ea = str(lead["EA"]).upper()
    lead_nea = str(lead["NEA"]).upper()
    lead_pni = float(pd.to_numeric(lead.get("P_NI", np.nan), errors="coerce"))

    print(f"[region] {chrom}:{start}-{end} (lead {LEAD_RSID} at {lead_bp})")
    print(f"[sumstats] significant SNPs in window (P_NI<{GW_SIG}): {len(sig_df)}")
    sig_df.to_csv(OUT_SIG, index=False)
    print(f"Wrote {OUT_SIG}")

    with open(REGION_INFO, "w") as f:
        f.write(f"lead_rsid\t{LEAD_RSID}\nchrom\t{chrom}\nlead_bp_hg19\t{lead_bp}\n")
        f.write(f"window_bp\t{WINDOW_BP}\nregion_hg19\t{chrom}:{start}-{end}\n")
        f.write(f"GW_SIG(P_NI)\t{GW_SIG}\nCLUMP_R2\t{CLUMP_R2}\nPROXY_R2\t{PROXY_R2}\nMAF_MIN\t{MAF_MIN}\n")

    # seeds = lead + all GWS in window
    seeds = [{
        "RSID": LEAD_RSID, "CHR": chrom, "BP": lead_bp, "EA": lead_ea, "NEA": lead_nea, "P_NI": lead_pni, "is_lead": True
    }]
    for _, row in sig_df.iterrows():
        rsid = str(row["RSID"])
        if rsid == LEAD_RSID:
            continue
        seeds.append({
            "RSID": rsid,
            "CHR": chrom,
            "BP": int(row["BP"]),
            "EA": str(row["EA"]).upper(),
            "NEA": str(row["NEA"]).upper(),
            "P_NI": float(row["P_NI"]),
            "is_lead": False,
        })
    seeds_df = pd.DataFrame(seeds).sort_values(["is_lead", "P_NI"], ascending=[False, True]).reset_index(drop=True)

    ensure_panel()
    eur_samples = eur_samples_from_panel()
    print(f"[1000G] EUR sample IDs loaded: {len(eur_samples)}")

    bcftools_extract_region(chrom, start, end)

    vcf = VCF(REGION_VCF)
    samples = np.array(vcf.samples)
    eur_mask = np.array([s in eur_samples for s in samples], dtype=bool)
    n_eur = int(eur_mask.sum())
    if n_eur < 50:
        raise RuntimeError(f"Too few EUR samples selected ({n_eur}). Check panel parsing.")
    print(f"[1000G] EUR samples in VCF: {n_eur}")

    needed_by_key = {allele_set_key(int(r.BP), r.EA, r.NEA): r.RSID for r in seeds_df.itertuples(index=False)}
    needed_by_rsid = {str(r.RSID): str(r.RSID) for r in seeds_df.itertuples(index=False)}

    seed_dosages = {}
    seed_match_info = {}

    for var in vcf:
        if not is_biallelic_snp(var):
            continue
        pos = int(var.POS)
        ref = var.REF.upper()
        alt = var.ALT[0].upper()
        key = allele_set_key(pos, ref, alt)

        matched = None
        if key in needed_by_key:
            matched = needed_by_key[key]
        elif var.ID in needed_by_rsid:
            matched = needed_by_rsid[var.ID]

        if matched is not None and matched not in seed_dosages:
            d = dosage_array(var.genotypes)[eur_mask]
            seed_dosages[matched] = d
            seed_match_info[matched] = {
                "seed_rsid": matched,
                "vcf_id": var.ID,
                "pos_hg19": pos,
                "ref_1000g": ref,
                "alt_1000g": alt,
            }

    if LEAD_RSID not in seed_dosages:
        raise RuntimeError(f"Lead {LEAD_RSID} not found in 1000G region VCF by allele-set or rsID.")

    found = sorted(seed_dosages.keys())
    missing = [rs for rs in seeds_df["RSID"].tolist() if rs not in seed_dosages]
    print(f"[1000G] seed SNPs requested: {len(seeds_df)}; found: {len(found)}; missing: {len(missing)}")

    seeds_df = seeds_df[seeds_df["RSID"].isin(found)].copy().reset_index(drop=True)
    if seeds_df.empty:
        raise RuntimeError("None of the seed SNPs were found in 1000G; cannot proceed.")

    index_list = []
    index_rows = []

    for r in seeds_df.itertuples(index=False):
        rsid = str(r.RSID)
        d = seed_dosages[rsid]

        same_signal = False
        for idx in index_list:
            r2_val = r2(seed_dosages[idx], d)
            if not np.isnan(r2_val) and r2_val >= CLUMP_R2:
                same_signal = True
                break

        if not same_signal:
            index_list.append(rsid)
            mi = seed_match_info.get(rsid, {})
            index_rows.append({
                "index_rsid": rsid,
                "pos_hg19": int(r.BP),
                "EA": str(r.EA),
                "NEA": str(r.NEA),
                "P_NI": float(r.P_NI),
                "is_lead": bool(r.is_lead),
                "matched_1000g_id": mi.get("vcf_id", ""),
                "ref_1000g": mi.get("ref_1000g", ""),
                "alt_1000g": mi.get("alt_1000g", ""),
            })

    index_df = pd.DataFrame(index_rows).sort_values(["is_lead", "P_NI"], ascending=[False, True])
    index_df.to_csv(OUT_INDEX, index=False)
    print(f"[signals] independent index SNPs (same signal if r2>={CLUMP_R2}): {len(index_df)}")
    print(f"Wrote {OUT_INDEX}")

    vcf = VCF(REGION_VCF)

    index_rsids = index_df["index_rsid"].astype(str).tolist()
    index_dos = [seed_dosages[rs] for rs in index_rsids]
    index_pos = {row.index_rsid: int(row.pos_hg19) for row in index_df.itertuples(index=False)}

    proxy_rows = []
    for var in vcf:
        if not is_biallelic_snp(var):
            continue
        d = dosage_array(var.genotypes)[eur_mask]
        maf = maf_from_dosage(d)
        if np.isnan(maf) or maf < MAF_MIN:
            continue

        pos = int(var.POS)
        ref = var.REF.upper()
        alt = var.ALT[0].upper()

        for idx_rsid, idx_d in zip(index_rsids, index_dos):
            r2_val = r2(idx_d, d)
            if np.isnan(r2_val) or r2_val < PROXY_R2:
                continue
            proxy_rows.append({
                "index_rsid": idx_rsid,
                "proxy_rsid": var.ID,
                "chr": str(var.CHROM),
                "pos_hg19": pos,
                "ref_1000g": ref,
                "alt_1000g": alt,
                "maf_eur": maf,
                "r2_eur": r2_val,
                "distance_to_index_bp": pos - index_pos[idx_rsid],
            })

    proxies_df = pd.DataFrame(proxy_rows).drop_duplicates(
        subset=["index_rsid", "proxy_rsid", "pos_hg19", "ref_1000g", "alt_1000g"]
    ).sort_values(["index_rsid", "r2_eur", "maf_eur"], ascending=[True, False, False])

    proxies_df.to_csv(OUT_PROXIES, index=False)
    print(f"[proxies] total proxies across all signals: {len(proxies_df)}")
    print(f"Wrote {OUT_PROXIES}")

if __name__ == "__main__":
    main()
