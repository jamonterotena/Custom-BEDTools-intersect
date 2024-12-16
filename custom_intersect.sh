#!/bin/bash

# bedtools intersect adapted for class type
## binary: add 0 / 1 if no overlap / overlap
## integer: 0, 1, 2, 3, ..., n
## ratio: 0 - 1
## numeric.length: 0.23413, 4324.123, 99.234 (average of basepair lengths)
## numeric.value: 0.23413, 4324.123, 99.234 (average of anything but basepair lengths)

# Define class
class=$3

# Check if more than 3 arguments are passed
if [ "$#" -gt 3 ]; then
    echo "Warning: More than 3 arguments passed." >&2
fi

# Check if the third argument is either "binary", "ratio", "ratio.value", "integer", "numeric" or "numeric.sum"
if [ -n "$class" ]; then
    if [[ "$class" != "binary" && "$class" != "ratio" && "$class" != "integer" && "$class" != "numeric.length" && "$class" != "numeric.value" && "$class" != "ratio.value" && "${class}" != "numeric.sum" ]]; then
        echo "Error: non-accepted class"
        # Print usage message
        echo "Usage: custom_intersect.sh crossovers feature class" >&2
        exit 1  # Exit with error status
    fi
else
    # Print usage message if the third argument is missing
    echo "Warning: Missing third argument 'class'." >&2
    echo "Usage: custom_intersect.sh crossovers feature class" >&2
    exit 1  # Exit with error status
fi

# Return the count of intersected intervals and replace counts >1 by 1 (binary: 0/1)
if [[ "${class}" == "binary" ]]; then
    # Convert counts >1 to 1
    bedtools intersect -a "$1" -b "$2" -c | awk 'BEGIN {OFS="\t"} { if ($NF > 1) $NF = 1; print }' | sponge "$1"
    head "$1"
fi

# Return the count of intersected intervals
if [[ "${class}" == "integer" ]]; then
    # Simply return count
    bedtools intersect -a "$1" -b "$2" -c | sponge "$1"
    head "$1"
fi

# Return the ratio or proportion of the genomic interval covered by the feature overlap
if [[ "${class}" == "ratio" ]]; then
    # Output overlaps in file 1, overlaps in file 2 plus the base pair length of the overlaps
    bedtools intersect -wao -a "$1" -b "$2" | awk -F'\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF-4; i++) {printf "%s\t", $i;} print $NF}' | sponge "$1"
    head "$1"
    # Divide last column by the interval length in file 1
    cat "$1" |  awk -v OFS="\t" '{ overlap_length = $NF; interval_length = $3 - $2; perc_overlap_length = overlap_length / interval_length; for (i=1; i<=NF-1; i++) {printf "%s\t", $i;} print perc_overlap_length}' | sponge "$1"
    head "$1"
    # Extract the number of columns of the file in order to know which is last column
    ncol=$(head -n 1 "$1" | tr "\t" "\n" | wc -l)
    ncol1=$((ncol-1))
    # Sum last column grouped by coordinates
    cat "$1" | bedtools groupby -g 1-"$ncol1" -c "$ncol" -o sum | sponge "$1"
fi

# Return the average length of the feature overlap within the genomic interval
if [[ "${class}" == "numeric.length" ]]; then
    # Output overlaps in file 1, overlaps in file 2 plus the base pair length of the overlaps
    bedtools intersect -wao -a "$1" -b "$2" | awk -F'\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF-4; i++) {printf "%s\t", $i;} print $NF}' | sponge "$1"
    head "$1"
    # Extract the number of columns of the file in order to know which is last column
    ncol=$(head -n 1 "$1" | tr "\t" "\n" | wc -l)
    ncol1=$((ncol-1))
    # Calculate mean of last column
    cat "$1" | bedtools groupby -g 1-"$ncol1" -c "$ncol" -o mean | sponge "$1"
fi

# Return the average value of the 4th column of file 2 (not the average length)
if [[ "${class}" == "numeric.value" || "${class}" == "ratio.value" ]]; then
    bedtools intersect -a "$1" -b "$2" -wa -wb -loj | head
    # Output overlaps in file 1 and overlaps in file 2
    bedtools intersect -a "$1" -b "$2" -wa -wb -loj | awk -F'\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF-4; i++) {printf "%s\t", $i;} printf "%.5f\n", $NF}' | sed "s/.$/0/g" | sponge "$1"
    head "$1"
    # Extract the number of columns of the file in order to know which is last column
    ncol=$(head -n 1 "$1" | tr "\t" "\n" | wc -l)
    ncol1=$((ncol-1))
    # Calculate mean of last column
    cat "$1" | bedtools groupby -g 1-"$ncol1" -c "$ncol" -o mean | sponge "$1"
fi

# Return the summed value of the 4th column of file 2
if [[ "${class}" == "numeric.sum" ]]; then
    # Output overlaps in file 1 and overlaps in file 2
    bedtools intersect -a "$1" -b "$2" -wa -wb -loj | awk -F'\t' 'BEGIN {OFS="\t"} { for (i=1; i<=NF-4; i++) {printf "%s\t", $i;} printf "%.5f\n", $NF}' | sed "s/.$/0/g" | sponge "$1"
    head "$1"
    # Extract the number of columns of the file in order to know which is last column
    ncol=$(head -n 1 "$1" | tr "\t" "\n" | wc -l)
    ncol1=$((ncol-1))
    # Calculate sum of last column
    cat "$1" | bedtools groupby -g 1-"$ncol1" -c "$ncol" -o sum | sponge "$1"
fi
