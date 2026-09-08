#!/usr/bin/env python3
"""Join baseline PBL samples to supplied coverage percentiles (original depth notebook)."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "preprocessing"))
from common import read_table, write_table


def sample_depths(monitoring, coverage, members=None):
    required = {"sample_id", "type", "percentile", "depth"}
    if not required.issubset(coverage):
        raise ValueError(f"Coverage requires {sorted(required)}; map sample IDs explicitly before running")
    samples = monitoring.loc[monitoring.sample_type.isin(["PBMC", "Normal", "PDWB"]) & monitoring.sample_day.le(0),
                             ["sample_id", "patient_id", "sample_day"]].drop_duplicates()
    if members is not None:
        samples = samples.loc[samples.patient_id.isin(members.patient_id)]
    depth = coverage.loc[coverage.type.eq("all") & coverage.percentile.eq(50),
                         ["sample_id", "type", "percentile", "depth"]]
    return samples.merge(depth, on="sample_id", how="left")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["monitoring", "coverage", "out"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--members")
    a = p.parse_args()
    write_table(sample_depths(read_table(a.monitoring), read_table(a.coverage),
                             read_table(a.members) if a.members else None), a.out)
