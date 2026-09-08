#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${PYTHON_BIN:-python3}"
if [[ $# -eq 1 && "$1" == "--example" ]]; then
  input_dir="$repo_dir/preprocessing/example_inputs"
  set -- "$input_dir/snv_monitoring.csv" "$input_dir/indel_monitoring.csv" \
    "$input_dir/snv_depth.csv" "$input_dir/indel_depth.csv" \
    "$input_dir/protein_coding_genes.tsv" "$repo_dir/outputs/example"
fi
if [[ $# -ne 6 ]]; then
  echo "Usage: ./run_preprocessing.sh SNV.csv INDEL.csv SNV_DEPTH.csv INDEL_DEPTH.csv GENES.tsv OUTPUT_DIR" >&2
  echo "   or: ./run_preprocessing.sh --example" >&2
  exit 2
fi
"$python_bin" "$repo_dir/preprocessing/01_filter_variants.py" \
  --snv "$1" --indel "$2" --snv-depth "$3" --indel-depth "$4" --genes "$5" --out "$6/filtered"
"$python_bin" "$repo_dir/preprocessing/02_prepare_analysis_tables.py" \
  --filtered "$6/filtered" --genes "$5" --out "$6/processed"
echo "Analysis inputs written to $6/processed"
