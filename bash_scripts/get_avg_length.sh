#!/bin/bash

input_file=$1
filename=$(basename -- "$input_file")
prefix="${filename%.*}"
output_file="${prefix}_stats.bed"

# Loop through each line of the input file
awk 'NR > 1 {
    # Extract chr, start, end, motif, and allele length columns
    chr = $1
    start = $2
    end = $3
    motif = $4
    n = 0
    sum = 0
    sum_sq = 0
    lengths = ""
    
    # Read allele lengths per individual
    for (i = 5; i <= NF; i++) {
        # Split allele pair
        split($i, alleles, ",")
        # Calculate allele length average per individual
        avg_length = (alleles[1] + alleles[2]) / 2
        lengths = lengths avg_length ","
        
        # Sum and sum of squares for variance calculation
        sum += avg_length
        sum_sq += avg_length^2
        n++
    }
    
    # Calculate mean, median, variance, and standard deviation
    mean = sum / n
    # Sort the allele lengths and calculate the median
    n_length = split(lengths, len_arr, ",")
    asort(len_arr)
    if (n_length % 2 == 1) {
        median = len_arr[(n_length + 1) / 2]
    } else {
        median = (len_arr[n_length / 2] + len_arr[(n_length / 2) + 1]) / 2
    }
    
    variance = (sum_sq / n) - (mean^2)
    stdev = sqrt(variance)
    
    print chr "\t" start "\t" end "\t" motif "\t" mean "\t" median "\t" variance "\t" stdev
}' "$input_file" > "$output_file"

echo "BED file containing average, median, variance and st. deviation of allele lengths per locus are saved to "${output_file}"
