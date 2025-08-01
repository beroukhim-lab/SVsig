"""
find_overlap_pairs.py

Given two files containing annotated top 5 gene_i and gene_j columns,
this script finds overlapping gene_i - gene_j pairs between the two datasets.

Usage:
    python find_overlap_pairs.py --sig_driv sig_driv.csv --driv driv.csv [--sep ,]

The input files must contain columns named:
    gene_i_annot_1 ... gene_i_annot_4
    gene_j_annot_1 ... gene_j_annot_4

Output is printed to stdout.
"""

import argparse
import pandas as pd
import sys

def load_pairs(df, prefix_i="gene_i_annot_", prefix_j="gene_j_annot_", num=4): 
    pairs = set()
    for i in range(1, num+1):
        for j in range(1, num+1):
            col_i = f"{prefix_i}{i}"
            col_j = f"{prefix_j}{j}"
            if col_i not in df.columns or col_j not in df.columns:
                continue  # skip missing columns silently
            sub = df[[col_i, col_j]].dropna()
            # Convert each row to tuple (gene_i, gene_j)
            pairs.update(tuple(x) for x in sub.values)
    return pairs

def main():
    parser = argparse.ArgumentParser(description="Find overlapping gene_i - gene_j pairs between two tables.")
    parser.add_argument("--pcawg", required=True, help="Path to pcawg file (CSV or TSV).")
    parser.add_argument("--ccle", required=True, help="Path to ccle file (CSV or TSV).")
    parser.add_argument("--sep", default=None,
                        help="Field separator for input files. If not provided, will infer from extension: "
                             "'.tsv' -> tab, else comma.")
    parser.add_argument("--top", type=int, default=None,
                        help="If provided, limit output to the first N overlapping pairs (after sorting).")
    args = parser.parse_args()

    # Determine separator
    def infer_sep(path, override):
        if override is not None:
            return override
        if path.lower().endswith(".tsv") or path.lower().endswith(".txt"):
            return "\t"
        return ","

    sep_sig = infer_sep(args.pcawg, args.sep)
    sep_driv = infer_sep(args.ccle, args.sep)

    try:
        sig_driv_df = pd.read_csv(args.pcawg, sep=sep_sig, dtype=str)
    except Exception as e:
        print(f"Error reading sig_driv file '{args.pcawg}': {e}", file=sys.stderr)
        sys.exit(1)

    try:
        driv_df = pd.read_csv(args.ccle, sep=sep_driv, dtype=str)
    except Exception as e:
        print(f"Error reading driv file '{args.ccle}': {e}", file=sys.stderr)
        sys.exit(1)

    sig_gene_pairs = load_pairs(sig_driv_df)
    driv_gene_pairs = load_pairs(driv_df)

    overlap_pairs = sig_gene_pairs.intersection(driv_gene_pairs)

    print(f"Number of overlapping gene_i - gene_j pairs: {len(overlap_pairs)}")
    if len(overlap_pairs) == 0:
        return

    print("Overlapping gene_i - gene_j pairs:")
    for idx, pair in enumerate(sorted(overlap_pairs)):
        if args.top is not None and idx >= args.top:
            break
        print(f"{pair[0]} - {pair[1]}")

if __name__ == "__main__":
    main()
