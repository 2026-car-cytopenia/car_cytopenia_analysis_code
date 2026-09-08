"""Shared table handling for the standalone analysis scripts."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ALLELE = ["CHR", "POS", "TUMOR", "REF"]
ANNOTATION = ALLELE + ["GENE", "TYPE"]
CH_GENES = ["ASXL1", "CBL", "CHEK2", "DNMT3A", "GNAS", "GNB1", "IDH1", "IDH2",
            "JAK2", "MYD88", "PPM1D", "SRSF2", "STAT3", "TET2", "U2AF1", "TP53", "SF3B1"]
NONFUNCTIONAL = ["ncRNA_exonic", "ncRNA_intronic", "intergenic", "intronic", "synonymous",
                 "upstream", "upstream;downstream", "UTR5", "UTR3", "downstream"]


def read_table(path, required=()):
    path = Path(path)
    df = pd.read_csv(path, sep="\t" if path.suffix in {".tsv", ".bed"} else ",", low_memory=False)
    df = df.loc[:, ~df.columns.str.startswith("Unnamed:")].copy()
    if "Progression" in df and "progression" not in df:
        df = df.rename(columns={"Progression": "progression"})
    missing = set(required) - set(df)
    if missing:
        raise ValueError(f"{path.name}: missing columns {sorted(missing)}")
    if "POS" in df:
        df["POS"] = pd.to_numeric(df.POS, errors="raise").astype("int64")
    for col in ["AF", "Depth", "sample_day", "progression"]:
        if col in df:
            df[col] = pd.to_numeric(df[col], errors="raise")
    if "AF" in df and ((df.AF.dropna() < 0) | (df.AF.dropna() > 100)).any():
        raise ValueError(f"{path.name}: AF must be a percentage between 0 and 100")
    return df


def write_table(df, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def write_json(value, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2) + "\n")


def variant_id(df, patient=False):
    # Preserve the original exported identifier convention; joins use columns.
    cols = ALLELE + (["patient_id"] if patient else [])
    return df[cols].astype(str).agg("".join, axis=1)


def attach_metadata(df, metadata=None):
    if metadata is None:
        return df
    meta = read_table(metadata, ["sample_id"])
    if meta.sample_id.duplicated().any():
        raise ValueError("Metadata must have one row per sample_id")
    if not df.sample_id.isin(meta.sample_id).all():
        raise ValueError("Metadata is missing one or more monitoring samples")
    # An explicit metadata table supersedes fields embedded in monitoring files.
    return df.drop(columns=[c for c in meta if c != "sample_id" and c in df]).merge(
        meta, on="sample_id", how="left", validate="many_to_one")


def restrict_regions(df, bed, position_base=0):
    """Intersect POS with half-open BED intervals; original notebooks use POS as-is."""
    intervals = pd.read_csv(bed, sep="\t", comment="#", header=None, usecols=[0, 1, 2])
    intervals.columns = ["CHR", "start", "end"]
    keep = np.zeros(len(df), dtype=bool)
    for chrom, rows in intervals.groupby("CHR"):
        idx = np.flatnonzero(df.CHR.to_numpy() == chrom)
        pos = df.iloc[idx].POS.to_numpy() - position_base
        for start, end in rows[["start", "end"]].itertuples(index=False, name=None):
            keep[idx] |= (pos >= start) & (pos < end)
    return df.loc[keep].copy()
