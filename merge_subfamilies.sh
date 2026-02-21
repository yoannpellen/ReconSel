#!/bin/bash

# Run as: bash merge_subfamilies.sh sf1 sf2 [sf3 ...]

if [ "$#" -lt 2 ]; then
    echo "Please give at least 2 subfamilies to merge: $0 sf1 sf2 [sf3 ...]"
	echo "	Subfamilies must be given simply as 9E, M, X, etc..."
    exit 1
fi

# Store given subfamilies into a list
sf_list=("$@")

# Build new subfamily name by joining each subfamily with '-'
new_sf=$(IFS=-; echo "${sf_list[*]}")

# Create output directory and initialize output fasta
new_sf_dir="sf_${new_sf}"
new_sf_fasta="${new_sf_dir}/${new_sf}.nucl.trimmed.fna"

mkdir -p "${new_sf_dir}"
> "${new_sf_fasta}"

# Loop through the subfamilies given as input and concatenate files
for sf in "${sf_list[@]}"; do
    ref_fasta="sf_${sf}/${sf}.nucl.trimmed.fna"

    if [ ! -f "${ref_fasta}" ]; then
        echo "Error: File not found: ${ref_fasta}"
        exit 1
    fi

    cat "${ref_fasta}" >> "${new_sf_fasta}"
done

# Rename original folders so we can get it back in case of mistake, otherwise it can just be deleted
for sf in "${sf_list[@]}"; do
    mv "sf_${sf}" "OLD.sf_${sf}"
done

echo "Merged ${#sf_list[@]} files into:"
echo "	${new_sf_fasta}"
