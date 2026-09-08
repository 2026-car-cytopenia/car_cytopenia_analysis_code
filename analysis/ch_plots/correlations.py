#!/usr/bin/env python3
"""Exact-day longitudinal VAF correlations from the corrected correlation notebook."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "preprocessing"))
import argparse
import numpy as np
import pandas as pd
from common import ALLELE, read_table, write_table, attach_metadata, variant_id

GENES = ["ASXL1", "CHEK2", "DNMT3A", "PPM1D", "SF3B1", "TET2", "TP53"]


def correlations(df, genes=GENES, min_days=3):
    if "detailed_sample_type" not in df:
        raise ValueError("Provide detailed_sample_type (cfDNA/PBMC/PDWB) in sample metadata; IDs are not parsed")
    d = df.dropna(subset=["AF"]).drop_duplicates(ALLELE + ["GENE", "sample_id"]).copy()
    if "var_id" not in d:
        d["var_id"] = variant_id(d, patient=True)
    output = {k: [] for k in ["cfDNA_vs_PBL", "cfDNA_vs_PBMC", "cfDNA_vs_PDWB", "PBMC_vs_PDWB"]}
    for patient, frame in d.groupby("patient_id", sort=False):
        plasma = frame.loc[frame.detailed_sample_type == "cfDNA"]
        pbmc = frame.loc[frame.detailed_sample_type == "PBMC"]
        pdwb = frame.loc[frame.detailed_sample_type == "PDWB"]
        normal = pdwb if len(pdwb) else pbmc
        for name, a, b in [("cfDNA_vs_PBL", plasma, normal), ("cfDNA_vs_PBMC", plasma, pbmc),
                           ("cfDNA_vs_PDWB", plasma, pdwb), ("PBMC_vs_PDWB", pbmc, pdwb)]:
            for v in np.intersect1d(a.var_id, b.var_id):
                aa = a.loc[a.var_id == v].drop_duplicates("sample_day").set_index("sample_day")
                bb = b.loc[b.var_id == v].drop_duplicates("sample_day").set_index("sample_day")
                days = np.intersect1d(aa.index, bb.index)
                gene = aa.GENE.iloc[0]
                if len(days) < min_days or gene not in genes:
                    continue
                x, y = aa.loc[days].AF.to_numpy(), bb.loc[days].AF.to_numpy()
                if np.std(x) == 0 or np.std(y) == 0:
                    continue
                output[name].append(dict(patient=patient, gene=gene, corr=float(np.corrcoef(x, y)[0, 1]),
                                         var_id=v, n_days=len(days)))
    return {k: pd.DataFrame(v, columns=["patient", "gene", "corr", "var_id", "n_days"]) for k, v in output.items()}


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--monitoring", required=True)
    p.add_argument("--metadata", help="Optional if detailed_sample_type is already in monitoring")
    p.add_argument("--out", required=True)
    a = p.parse_args()
    for name, frame in correlations(attach_metadata(read_table(a.monitoring), a.metadata)).items():
        write_table(frame, Path(a.out) / f"{name}.csv")
