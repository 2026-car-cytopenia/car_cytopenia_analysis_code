#!/usr/bin/env python3
"""Create de novo CH calls, longitudinal monitoring and baseline CH tables."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "preprocessing"))
import argparse
import pandas as pd
from common import ALLELE, CH_GENES, NONFUNCTIONAL, read_table, write_table, variant_id, restrict_regions


def prepare_ch(monitoring, calls, snv_depth, indel_depth, *, cohort="precar", members=None, bed=None,
               position_base=0, genes=CH_GENES, min_reads=5):
    keys = ALLELE + ["sample_id"]
    # Calls are caller-positive alleles, not every allele in the monitoring table.
    eligibility = monitoring.sample_type.isin(["PBMC", "Normal", "PDWB"])
    if cohort != "pretx":
        if "progression" not in monitoring:
            raise ValueError("Post-treatment CH discovery requires progression")
        eligibility |= monitoring.sample_type.eq("cfDNA") & monitoring.sample_day.ge(28) & monitoring.progression.eq(0)
    raw = monitoring.loc[eligibility].merge(calls[keys].drop_duplicates(), on=keys, how="inner", validate="many_to_one")
    raw = raw.drop_duplicates()
    longitudinal = monitoring.merge(raw[ALLELE + ["patient_id"]].drop_duplicates(),
                                    on=ALLELE + ["patient_id"], how="inner", validate="many_to_one")
    longitudinal["var_id"] = variant_id(longitudinal, patient=True)
    bases = ["A", "C", "G", "T"]
    is_snv = raw.TUMOR.isin(bases) & raw.REF.isin(bases)
    enriched = []
    for data, depth, snv in [(raw.loc[is_snv], snv_depth, True), (raw.loc[~is_snv], indel_depth, False)]:
        # Preserve original join multiplicity. Existing annotation-specific depth
        # rows are part of the reference workflow; do not silently change weighting.
        depth = depth[keys + ["Depth"]].copy()
        depth.Depth = depth.Depth.fillna(0).astype(int)
        d = data.drop(columns="Depth", errors="ignore").merge(depth, on=keys, how="left")
        if snv:
            d.Depth = d.Depth * d.AF / 100  # SNV table reports total depth.
        # Indel Depth is already alternate-supporting depth in the source pipeline.
        enriched.append(d)
    ch = pd.concat(enriched, ignore_index=True)
    ch.sample_type = ch.sample_type.replace({"Normal": "PBMC", "PDWB": "PBMC"})
    ch = ch.loc[ch.Depth.ge(min_reads) & ~ch.TYPE.isin(NONFUNCTIONAL) & ch.GENE.isin(genes) & ch.sample_day.le(0)].copy()
    ch["var_id"] = variant_id(ch)
    if members is not None:
        ch = ch.loc[ch.patient_id.isin(members.patient_id)]
    before_regions = ch.copy()
    if bed:
        ch = restrict_regions(ch, bed, position_base)
        ch["valid"] = True
    return raw, longitudinal, before_regions, ch


def write_outputs(result, output):
    names = ["raw_calls", "ch_monitoring", "ch_baseline_all_regions", "ch_baseline"]
    for name, df in zip(names, result):
        write_table(df, Path(output) / f"{name}.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for key in ["monitoring", "calls", "snv-depth", "indel-depth", "out"]:
        p.add_argument("--" + key, required=True)
    p.add_argument("--cohort", choices=["precar", "pretx"], default="precar")
    p.add_argument("--members")
    p.add_argument("--bed")
    p.add_argument("--position-base", type=int, choices=[0, 1], default=0)
    a = p.parse_args()
    write_outputs(prepare_ch(read_table(a.monitoring), read_table(a.calls), read_table(a.snv_depth),
                  read_table(a.indel_depth), cohort=a.cohort,
                  members=read_table(a.members) if a.members else None, bed=a.bed, position_base=a.position_base), a.out)

