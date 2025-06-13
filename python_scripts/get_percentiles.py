#!/usr/bin/env python3
"""
Author: Adam, Carolina de Lima
Date: 2025
Purpose:
    Calculate summary statistics (mean, variance, 5th and 95th percentiles)
    for allele metrics per locus from a tab-separated input file.

Input format:
    A tab-separated file where:
    - First column locus identifier
    - Subsequent columns are numeric allele statistics (can contain '.' for missing data)

    *Lines with fewer than two numeric values after filtering missing data are skipped.
    
Output:
    A tab-separated file with columns:
    locus_id, mean, variance, percentile_5, percentile_95
"""

import argparse
import statistics

parser = argparse.ArgumentParser()
parser.add_argument("-i", required=True, help="Input file")
args = parser.parse_args()

outfile = args.i.rsplit(".", 1)[0] + "_percentiles.bed"

def percentile(data, percent):
    k = (len(data) - 1) * (percent / 100)
    f = int(k)
    c = f + 1
    if c >= len(data):
        return data[f]
    return data[f] + (data[c] - data[f]) * (k - f)

with open(args.i) as fin, open(outfile, "w") as fout:
    fout.write("locus_id\tmean\tvariance\tpercentile_5\tpercentile_95\n")
    for line in fin:
        parts = line.strip().split('\t')
        if len(parts) < 2:
            continue
        locus_id = parts[0]
        # Filter out missing values ('.') and convert rest to float
        try:
            values = [float(x) for x in parts[1:] if x != "."]
        except ValueError:
            print(f"Warning: Could not parse numeric values in line: {line.strip()}")
            continue
        if len(values) < 2:
            print(f"Warning: Not enough data for locus {locus_id}, skipping")
            continue
        values.sort()
        mean = statistics.mean(values)
        variance = statistics.variance(values)
        p5 = percentile(values, 5)
        p95 = percentile(values, 95)
        fout.write(f"{locus_id}\t{mean:.3f}\t{variance:.3f}\t{p5:.3f}\t{p95:.3f}\n")

print(f"Output written to {outfile}")
