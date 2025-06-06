#!/usr/bin/env python3
"""
merge_maf.py: Merge all MAF files in an input directory, extracting specific columns,
creating sample-specific frequency/depth columns, deduplicating rows on key fields,
and producing a two-line header indicating sample names.

Usage:
    python merge_maf.py -i /path/to/input_folder
"""
import argparse
import pandas as pd
import sys
from pathlib import Path

def parse_args():
    parser = argparse.ArgumentParser(
        description="Merge MAFs into a matrix with sample-specific freq/depth and a two-line header."
    )
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='Directory containing the .maf files'
    )
    return parser.parse_args()

def detect_sample_name(path: Path) -> str:
    """Detect sample by keyword, case-insensitive."""
    name = path.name.lower()
    if 'ffpe' in name:
        return 'FFPE'
    if 'leuco' in name:
        return 'Leuco'
    if 'ctdna' in name or 'cfdna' in name:
        return 'ctDNA'
    # fallback: prefix before “.hard”
    return path.name.split('.hard')[0]

def read_maf(path: Path) -> pd.DataFrame:
    # Keep exactly these 9 columns: 7 identity/annotation + 2 measures
    cols = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'HGVSc',
        'HGVSp',
        'Consequence',
        'Variant_Frequencies',
        'Total_Depth'
    ]
    df = pd.read_csv(path, sep='\t', comment='#', dtype=str, low_memory=False)
    missing = set(cols) - set(df.columns)
    if missing:
        print(f"Error: {path.name} missing columns: {', '.join(sorted(missing))}", file=sys.stderr)
        sys.exit(1)
    return df[cols].copy()

def main():
    args = parse_args()
    input_dir = Path(args.input_dir)
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory.", file=sys.stderr)
        sys.exit(1)

    maf_paths = sorted(input_dir.glob("*.maf"))
    if not maf_paths:
        print(f"Error: no .maf files found in {input_dir}.", file=sys.stderr)
        sys.exit(1)

    # Columns to keep
    orig_cols = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'HGVSc',
        'HGVSp',
        'Consequence'
    ]
    key_cols = [
        'Hugo_Symbol',
        'Chromosome',
        'Start_Position',
        'End_Position',
        'HGVSc',
        'HGVSp'
    ]
    extra_cols = ['Variant_Frequencies', 'Total_Depth']

    # Read each file, detect sample, preserve order
    sample_order = []
    sample_dfs = {}
    for path in maf_paths:
        sample = detect_sample_name(path)
        if sample not in sample_order:
            sample_order.append(sample)
        sample_dfs[sample] = read_maf(path)

    # Build master list of unique variants by the key columns
    all_keys = pd.concat(
        [df[key_cols] for df in sample_dfs.values()],
        ignore_index=True
    ).drop_duplicates().reset_index(drop=True)

    # Retrieve the annotation columns for each variant
    ann = pd.concat(
        [df[orig_cols] for df in sample_dfs.values()],
        ignore_index=True
    ).drop_duplicates(subset=key_cols).reset_index(drop=True)

    master = all_keys.merge(ann, on=key_cols, how='left')

    # Merge in each sample's frequency/depth
    for sample in sample_order:
        df = sample_dfs[sample]
        rename_map = {
            'Variant_Frequencies': f"{sample}_Variant_Frequencies",
            'Total_Depth':          f"{sample}_Total_Depth"
        }
        df_renamed = df[key_cols + extra_cols].rename(columns=rename_map)
        master = master.merge(df_renamed, on=key_cols, how='left')

    # Prepare output directory and file
    out_dir = input_dir / 'Matrice_Merge'
    out_dir.mkdir(exist_ok=True)
    out_file = out_dir / 'merged.maf'

    # Build the two-line header
    header1 = [''] * len(orig_cols)
    for sample in sample_order:
        header1 += [sample, sample]
    header2 = orig_cols.copy()
    for _ in sample_order:
        header2 += extra_cols

    # Write output
    with open(out_file, 'w', newline='') as fh:
        fh.write('\t'.join(header1) + '\n')
        fh.write('\t'.join(header2) + '\n')
        master.to_csv(fh, sep='\t', header=False, index=False)

    print(f"Wrote merged matrix to {out_file} "
          f"({len(master)} variants x {len(header1)} columns).")

if __name__ == '__main__':
    main()

