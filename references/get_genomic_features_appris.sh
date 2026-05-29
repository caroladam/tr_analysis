#!/bin/bash

# ---
# Title: get_genomic_features_appris.sh
# Date: 2025
# Author: Adam, Carolina L.
# Purpose: Determine CHM13 genomic features from annotation file (GFF)
# Inputs:
    # A GFF annotation file
    # A text file with chromosome lengths
# Outputs:
    # BED files containing genomic coordinates for "CDS", "5' UTR", "3' UTR", "promoter", "intron", and "intergenic" regions
# ---

usage() {
    echo "Usage: $0 <annotation_file> <chr_size_file>"
    exit 1
}

INPUT_GFF=$1
CHROM_SIZES=$2

# Check input
if [[ ! -f $INPUT_GFF ]]; then
    echo "Error: Annotation file not found."
    usage
fi

if [[ ! -f $CHROM_SIZES ]]; then
    echo "Error: Chromosome size file not found."
    usage
fi

# Create output directory
OUTDIR="feature_beds"
mkdir -p "$OUTDIR"

# Filter for appris_principal annotations
grep "appris_principal" "$INPUT_GFF" > appris.gff3

# Extract full exons
awk '$3 == "exon" { $4 = $4 - 1; print $1, $4, $5, $2, $3, $9, ".", $7}' OFS='\t' appris.gff3 > "$OUTDIR/exons.bed"

# Extract coding sequences
awk '$3 == "CDS" { $4 = $4 - 1; print $1, $4, $5, $2, $3, $9, ".", $7}' OFS='\t' appris.gff3 > "$OUTDIR/cds.bed"

# Extract introns
awk '$3 == "intron" { $4 = $4 - 1; print $1, $4, $5, $2, $3, $9, ".", $7}' OFS='\t' appris.gff3 > "$OUTDIR/introns.bed"

# Extract TSS from transcript entries
awk '$3 == "transcript"' appris.gff3 | awk 'BEGIN{OFS="\t"} {
  if ($7 == "+") print $1, $4-1, $4, $9, ".", $7;
  else if ($7 == "-") print $1, $5-1, $5, $9, ".", $7;
}' > "$OUTDIR/tss.bed"

# Define promoters 1Kb upstream of TSS
bedtools slop -i "$OUTDIR/tss.bed" -g "$CHROM_SIZES" -l 1000 -r 0 -s > "$OUTDIR/promoters.bed"

# Infer 3' UTRs

EXONS_TMP=$(mktemp)
STOP_TMP=$(mktemp)
STOP_INFO=$(mktemp)

awk '$3 == "exon"' appris.gff3 > "$EXONS_TMP"
awk '$3 == "stop_codon"' appris.gff3 > "$STOP_TMP"

awk 'BEGIN{OFS="\t"} {
    match($9, /transcript_id=([^;]+)/, arr);
    transcript=arr[1];
    print transcript, $1, $4, $5, $7;
}' "$STOP_TMP" > "$STOP_INFO"

awk -v OFS='\t' '
FNR==NR {
    stop[$1] = $2"\t"$3"\t"$4"\t"$5;
    next
}
{
    match($9, /transcript_id=([^;]+)/, a);
    tid = a[1];
    if (tid in stop) {
        split(stop[tid], s, "\t");
        chr = s[1]; sc = s[2]; ec = s[3]; strand = s[4];

        exon_start = $4;
        exon_end = $5;

        if (strand == "+") {
            if (exon_end > ec) {
                utr_start = (exon_start > ec) ? exon_start : ec + 1;
                utr_end = exon_end;
                if (utr_start <= utr_end) {
                    $4 = utr_start - 1;
                    $5 = utr_end;
                    print $0;
                }
            }
        } else if (strand == "-") {
            if (exon_start < sc) {
                utr_start = exon_start;
                utr_end = (exon_end < sc) ? exon_end : sc - 1;
                if (utr_start <= utr_end) {
                    $4 = utr_start - 1;
                    $5 = utr_end;
                    print $0;
                }
            }
        }
    }
}' "$STOP_INFO" "$EXONS_TMP" > "$OUTDIR/3utr_inferred.bed"

# Infer 5' UTRs

START_TMP=$(mktemp)
START_INFO=$(mktemp)

awk '$3 == "start_codon"' appris.gff3 > "$START_TMP"

awk 'BEGIN{OFS="\t"} {
    match($9, /transcript_id=([^;]+)/, arr);
    transcript=arr[1];
    print transcript, $1, $4, $5, $7;
}' "$START_TMP" > "$START_INFO"

awk -v OFS='\t' '
FNR==NR {
    start[$1] = $2"\t"$3"\t"$4"\t"$5;
    next
}
{
    match($9, /transcript_id=([^;]+)/, a);
    tid = a[1];
    if (tid in start) {
        split(start[tid], s, "\t");
        chr = s[1]; sc = s[2]; ec = s[3]; strand = s[4];

        exon_start = $4;
        exon_end = $5;

        if (strand == "+") {
            if (exon_start < sc) {
                utr_start = exon_start;
                utr_end = (exon_end < sc) ? exon_end : sc - 1;
                if (utr_start <= utr_end) {
                    $4 = utr_start - 1;
                    $5 = utr_end;
                    print $0;
                }
            }
        } else if (strand == "-") {
            if (exon_end > ec) {
                utr_start = (exon_start > ec) ? exon_start : ec + 1;
                utr_end = exon_end;
                if (utr_start <= utr_end) {
                    $4 = utr_start - 1;
                    $5 = utr_end;
                    print $0;
                }
            }
        }
    }
}' "$START_INFO" "$EXONS_TMP" > "$OUTDIR/5utr_inferred.bed"

# Cleanup
rm "$EXONS_TMP" "$STOP_TMP" "$STOP_INFO" "$START_TMP" "$START_INFO"

# Subtract UTRs from exons
bedtools intersect -v -a "$OUTDIR/exons.bed" -b "$OUTDIR/3utr_inferred.bed" | bedtools intersect -v -a - -b "$OUTDIR/5utr_inferred.bed" > "$OUTDIR/exons_no_utr.bed"

echo "Processing done."
echo "Features saved in: $OUTDIR/"
