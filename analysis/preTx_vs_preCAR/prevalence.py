#!/usr/bin/env python3
"""Summarize baseline CH using the study's low/high VAF bins and cohort denominator."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "preprocessing"))
import argparse
import pandas as pd
from common import CH_GENES, read_table, write_table


def prevalence(ch, members, source, genes=CH_GENES, min_reads=5, high_af=2):
    if not {"sample_id", "patient_id"}.issubset(members):
        raise ValueError("Cohort members must include sample_id and patient_id, including CH-negative participants")
    total = members.patient_id.nunique()
    if total == 0:
        raise ValueError("The prevalence cohort is empty")
    d = ch.loc[ch.sample_id.isin(members.sample_id) & ch.sample_day.le(0) & ch.sample_type.isin(["PBMC", "Normal", "PDWB"])].copy()
    if len(d):
        d = d.loc[d.sample_day.eq(d.groupby("patient_id").sample_day.transform("min"))]
    d = d.loc[d.Depth.ge(min_reads) & d.GENE.isin(genes)].copy()
    d["low_af"] = d.AF.lt(high_af).astype(int)
    grouped = d.groupby(["patient_id", "GENE", "low_af"]).Depth.sum().unstack(fill_value=0).reindex(columns=[0, 1], fill_value=0)
    grouped.columns = ["depth_high_af", "depth_low_af"]
    rows = []
    combined = grouped.groupby(level="patient_id").sum()
    high = int(combined.depth_high_af.ge(min_reads).sum())
    low = int((combined.depth_high_af.lt(min_reads) & combined.depth_low_af.ge(min_reads)).sum())
    rows.append(dict(gene="Any", low_af=low, high_af=high, wildtype=total-low-high))
    for gene in genes:
        g = grouped.loc[grouped.index.get_level_values("GENE") == gene]
        # Preserve notebook behavior: per-gene low/high bins are counted separately.
        low = int(g.depth_low_af.ge(min_reads).sum())
        high = int(g.depth_high_af.ge(min_reads).sum())
        rows.append(dict(gene=gene, low_af=low, high_af=high, wildtype=total-low-high))
    result = pd.DataFrame(rows)
    result["source"] = source
    result["n_total"] = total
    counts = d.groupby("patient_id").var_id.nunique().reindex(members.patient_id.unique(), fill_value=0).rename("unique_mutations").reset_index()
    counts["source"] = source
    return result, grouped.reset_index(), counts


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for k in ["ch", "members", "source", "out"]:
        p.add_argument("--" + k, required=True)
    a = p.parse_args()
    summary, grouped, counts = prevalence(read_table(a.ch), read_table(a.members), a.source)
    for name, d in [("prevalence", summary), ("patient_gene_depth", grouped), ("mutation_counts", counts)]:
        write_table(d, Path(a.out) / f"{name}.csv")

