#!/usr/bin/env bash

# Check if at least 3 arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <input_directory> <threshold1> [<threshold2> ...] <output_directory_net>"
    exit 1
fi

# Initialize variables
input_directory="$1"
output_directory_net="${@: -1}"  # Last argument
threshold_list=("${@:2:$#-2}")   # All arguments except first and last

# Ensure the output directory exists
mkdir -p "$output_directory_net"

# Check if input directory exists
if [ ! -d "$input_directory" ]; then
    echo "Error: Input directory $input_directory does not exist."
    exit 1
fi

# Check which is the current working directory
cwd=$(pwd)
echo "Current working directory: $cwd"

# Parse the config file content as a variable
# Read the path from the .ini file
alphabridge_path=$(grep '^alphabridge_path' config_alphabridge_path.ini | cut -d '=' -f 2 | xargs)

# PRINT the alphabridge path to verify it's being extracted correctly
echo "Alphabridge path extracted: $alphabridge_path"

# Check if the path is valid
if [ -z "$alphabridge_path" ]; then
    echo "Error: Alphabridge path not found in config_alphabridge_path.ini."
    exit 1
fi

# Check if the alphabridge script exists at the path
if [ ! -f "$alphabridge_path" ]; then
    echo "Error: Alphabridge script not found at the specified path: $alphabridge_path."
    exit 1
fi

# Step 1: Running alphabridge, NO THRESHOLD HERE
echo "Running: python3 $alphabridge_path -i $input_directory"
python3 "$alphabridge_path" -i "$input_directory"

echo "Step 1 is done!"

# Step 2: Run the network script
echo "Running: python3 $cwd/src/network_int.py -i $input_directory -o $output_directory_net -t ${threshold_list[*]}"
python3 "$cwd/src/network_int.py" -i "$input_directory" -o "$output_directory_net" -t "${threshold_list[@]}"

echo "Done! Check your input folder for Alphabridge output and your output directory for the network :)"




