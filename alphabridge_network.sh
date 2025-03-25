#!/usr/bin/env bash

#check which is the current work directory --> variable che andra aggiunta a + /src/network_int.py 

# Check if all arguments are provided
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <input_directory> <threshold> <output_directory_project> <binary_or_complex>"
    exit 1
fi

# Initialize variables
input_directory="$1"
threshold=$(echo "$2" | xargs) # Remove unwanted spaces
output_directory_project="$3"
binary_or_complex="$4"

# Ensure the output directory exists
mkdir -p "$output_directory_project"

# Check if input directory contains only files (no subdirectories)
contains_only_files=true
for item in "$input_directory"/*; do
    if [ -d "$item" ]; then
        contains_only_files=false
        break
    fi
done
##just complex for noww
#running alphabridge 
if [ "$contains_only_files" = true ]; then
    # Process the single directory
    echo "Running: python3 define_interfaces.py -i $input_directory -t $threshold"
    python3 /local_data/carmen/zata/AlphaBridge/define_interfaces.py -i "$input_directory" -t "$threshold" || exit 1
else
    # Iterate over subdirectories
    for dir in "$input_directory"/*; do
        [ -d "$dir" ] || continue
        echo "Processing subdirectory: $dir"
        python3 /local_data/carmen/zata/AlphaBridge/define_interfaces.py -i "$dir" -t "$threshold" || exit 1
    done
fi

# Run the appropriate main script based on binary_or_complex argument
if [ "$binary_or_complex" = "bin" ]; then
    echo "Running: python3 binary_interactions.py -i $input_directory -o $output_directory_project -t $threshold"
    python3 /local_data/carmen/code/interactive_network/src/binary_interactions.py -i "$input_directory" -o "$output_directory_project" -t "$threshold" || exit 1
elif [ "$binary_or_complex" = "complex" ]; then
    echo "Running: python3 network_int.py -i $input_directory -o $output_directory_project -t $threshold"
    python3 /local_data/carmen/code/interactive_network/src/network_int.py -i "$input_directory" -o "$output_directory_project" -t "$threshold" || exit 1
else
    echo "Error: Invalid value for binary_or_complex. Use 'bin' or 'complex'."
    exit 1
fi

