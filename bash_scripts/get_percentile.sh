#!/bin/bash

# Usage: ./get_percentile.sh -i input_file.txt -p 90

while getopts "i:p:" opt; do
  case $opt in
    i) infile="$OPTARG";;
    p) perc="$OPTARG";;
    *) echo "Usage: $0 -i input_file -p percentile_range (e.g., 90)"; exit 1;;
  esac
done

if [ -z "$infile" ] || [ -z "$perc" ]; then
  echo "Input file with allele lengths per locus and percentile range are required"
  echo "Usage: $0 -i input_file -p percentile_range (e.g., 90)"
  exit 1
fi

# Compute lower and upper percentiles
lower=$(echo "scale=3; (100 - $perc) / 2" | bc)
upper=$(echo "scale=3; 100 - $lower" | bc)

outfile="${infile%.*}_percentiles_P${perc}.bed"

{
  echo -e "locus_id\tmean\tpercentile_${lower}\tpercentile_${upper}"
  while IFS=$'\t' read -r locus_id rest; do
    values=($rest)
    sorted=($(printf "%s\n" "${values[@]}" | sort -n))
    count=${#sorted[@]}

    # Calculate index positions
    lower_index=$(printf "%.0f\n" "$(echo "($count - 1) * $lower / 100" | bc -l)")
    upper_index=$(printf "%.0f\n" "$(echo "($count - 1) * $upper / 100" | bc -l)")

    p_lower=${sorted[$lower_index]}
    p_upper=${sorted[$upper_index]}

    sum=0
    for val in "${sorted[@]}"; do
      sum=$((sum + val))
    done
    mean=$(echo "scale=3; $sum / $count" | bc)

    echo -e "${locus_id}\t${mean}\t${p_lower}\t${p_upper}"
  done < "$infile"
} > "$outfile"

echo "Output written to $outfile"
