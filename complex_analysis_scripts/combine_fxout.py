#!/usr/bin/env python3
"""
FoldX5 output concat

Requirements:
- Pandas

Usage:
    python combine_fxout.py path/to/folder [--output-path ./output.tsv]
"""

import argparse
import os
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine multiple fxout files"
    )
    parser.add_argument(
        'folder',
        help='Folder containing fxout files'
    )
    parser.add_argument(
        '--output-path',
        default='./output.tsv',
        help='Output path for combined file (default: ./output.tsv)'
    )
    return parser.parse_args()

def get_complex_files(file_paths):
    return [file for file in file_paths if file.startswith('Interaction')]

def main():
    args = parse_args()

    all_files = os.listdir(args.folder)
    matched_files = get_complex_files(all_files)
    full_paths = [os.path.join(args.folder, file) for file in matched_files]
    dfs = [pd.read_csv(file, sep='\t', skiprows=8) for file in full_paths]
    combined_df = pd.concat(dfs)

    combined_df.to_csv(args.output_path, sep="\t", index=False)

if __name__ == "__main__":
    main()