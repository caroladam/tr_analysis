#!/bin/bash

# Usage: ./calculate_variance_sd.sh -i input_file.txt

while getopts "i:" opt; do
  case $opt in
    i) infile="$OPTARG";;
    *) echo "Usage: $0 -i input_file"; exit 1;;
  esac
done

if [ -z "$infile" ]; then
  echo "Input file with allele length per locus is required"
  echo "Usage: $0 -i input_file"
  exit 1
fi

outfile="${infile%.*}_variance_sd.bed"

{
  echo -e "locus_id\tvariance\tstd_dev"
  while IFS=$'\t' read -r locus_id rest; do
    values=($rest)

    # Calculate the mean
    count=${#values[@]}
    sum=0
    for val in "${values[@]}"; do
      sum=$((sum + val))
    done
    mean=$(echo "scale=3; $sum / $count" | bc)

    # Calculate variance
    variance_sum=0
    for val in "${values[@]}"; do
      variance_sum=$(echo "$variance_sum + ($val - $mean)^2" | bc)
    done
    variance=$(echo "scale=3; $variance_sum / $count" | bc)

    # Calculate standard deviation
    std_dev=$(echo "scale=3; sqrt($variance)" | bc)

    # Output the result
    echo -e "${locus_id}\t${variance}\t${std_dev}"
  done < "$infile"
} > "$outfile"

echo "Output written to $outfile"

