#!/usr/bin/env python3
"""Filter concatenated monitoring files, following the final filtering notebook."""
import argparse
from pathlib import Path
import pandas as pd
from common import ALLELE, ANNOTATION, read_table, write_table, write_json, attach_metadata


def filter_monitoring(snv, indel, genes):
    excluded = ["ncRNA_exonic", "ncRNA_intronic", "intergenic"]
    snv = snv.loc[snv.GENE.isin(genes) & ~snv.TYPE.isin(excluded)].copy()
    indel = indel.loc[indel.GENE.isin(genes) & ~indel.TYPE.isin(excluded)].copy()
    # Original indel screen: remove an annotation if AF > 0.5 in >80% of patients.
    n = indel.patient_id.nunique()
    counts = indel.loc[indel.AF > .5].groupby(ANNOTATION, dropna=False).patient_id.nunique()
    bad = counts[counts > .8 * n].reset_index()[ANNOTATION]
    if len(bad):
        indel = indel.merge(bad.assign(_exclude=True), on=ANNOTATION, how="left")
        indel = indel.loc[indel._exclude.isna()].drop(columns="_exclude")
    return snv.drop_duplicates(), indel, bad


def run(snv_path, indel_path, snv_depth_path, indel_depth_path, gene_path, outdir, metadata=None):
    out = Path(outdir)
    required = ANNOTATION + ["SAMPLE", "sample_id", "patient_id", "AF", "sample_day", "sample_type"]
    snv = attach_metadata(read_table(snv_path, required), metadata)
    indel = attach_metadata(read_table(indel_path, required), metadata)
    genes = pd.read_csv(gene_path, sep="\t", usecols=["symbol"])["symbol"].dropna().tolist()
    n_before = [len(snv), len(indel)]
    snv, indel, bad = filter_monitoring(snv, indel, genes)
    write_table(snv, out / "snv_filtered.csv")
    write_table(indel, out / "indel_filtered.csv")
    write_table(bad, out / "excluded_common_indels.csv")
    for name, path, variants in [("snv", snv_depth_path, snv), ("indel", indel_depth_path, indel)]:
        depth = attach_metadata(read_table(path, ALLELE + ["sample_id", "Depth"]), metadata)
        keys = ALLELE + ["sample_id"]
        selected = depth.merge(variants[keys].drop_duplicates(), on=keys, how="inner", validate="many_to_one")
        write_table(selected, out / f"{name}_depth_filtered.csv")
    write_json({"input_rows": dict(zip(["snv", "indel"], n_before)),
                "filtered_rows": {"snv": len(snv), "indel": len(indel)},
                "excluded_indel_annotations": len(bad)}, out / "filter_counts.json")
    return snv, indel, genes


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["snv", "indel", "snv-depth", "indel-depth", "genes", "out"]:
        p.add_argument("--" + name, required=True)
    p.add_argument("--metadata")
    a = p.parse_args()
    run(a.snv, a.indel, a.snv_depth, a.indel_depth, a.genes, a.out, a.metadata)
