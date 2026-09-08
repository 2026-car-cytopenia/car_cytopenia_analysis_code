#!/usr/bin/env python3
"""Combine filtered monitoring and prepare the corrected enrichment input."""
import argparse
from pathlib import Path
import pandas as pd
from common import ALLELE, read_table, write_table


def prepare_enrichment(snv, indel, genes):
    d = pd.concat([snv.assign(variant_class="snv"), indel.assign(variant_class="indel")], ignore_index=True)
    d.AF = d.AF.fillna(0)
    d = d.loc[~d.selector.fillna("").str.contains("MRD", regex=False)].copy()
    d.GENE = d.GENE.str.split(";").str[0]
    # Corrected study script removes synonymous annotations. Intronic/UTR variants
    # are not removed by that particular filter; CH classification is separate.
    d = d.loc[~d.TYPE.fillna("").str.contains("synonymous", regex=False)]
    d = d.drop_duplicates(ALLELE + ["SAMPLE"], keep="first")
    return d.loc[d.GENE.isin([g for g in genes if "LINC" not in g])].copy()



def run(filtered, genes_path, output):
    filtered, output = Path(filtered), Path(output)
    snv = read_table(filtered / "snv_filtered.csv")
    indel = read_table(filtered / "indel_filtered.csv")
    genes = pd.read_csv(genes_path, sep="\t", usecols=["symbol"])["symbol"].dropna().tolist()
    write_table(pd.concat([snv, indel], ignore_index=True), output / "monitoring.csv")
    write_table(prepare_enrichment(snv, indel, genes), output / "enrichment_input.csv")
    for kind in ["snv", "indel"]:
        write_table(read_table(filtered / f"{kind}_depth_filtered.csv"), output / f"{kind}_depth.csv")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    for name in ["filtered", "genes", "out"]:
        p.add_argument("--" + name, required=True)
    a = p.parse_args()
    run(a.filtered, a.genes, a.out)
