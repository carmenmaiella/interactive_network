#!/usr/bin/env bash

# Check if all arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_directory> <threshold> <output_directory_net>"
    exit 1
fi

# Initialize variables
input_directory="$1"
threshold=$(echo "$2" | xargs) # Remove unwanted spaces
output_directory_net="$3"

# Ensure the output directory exists
mkdir -p "$output_directory_net"

#check which is the current work directory 
cwd=$(pwd)
echo "current working directory: $cwd"

#parsing the confg file content as a variable 
# path to the conig 
config_file="$cwd/config_alphabridge_path.ini"
#grep -vE '^#|^$' "$config_file" --> ensures only non-empty, non-comment lines are passed through
#head -n 1: Takes only the first valid line from the filtered output
#path=$( ... ) -> captures the command’s output and stores it in the path variable
#creating the varibale
path=$(grep -vE '^#|^$' "$config_file" | head -n 1)

echo "current used path for the define_interface.py $path"


#STEP 1: running alphabridge 
echo "Running: python3 define_interfaces.py -i $input_directory -t $threshold"
#this is the path that the user should be able to change
python3 "$cwd"/config_file -i "$input_directory" -t "$threshold"

echo "step 1 is done!"

#STEP 2: run the network script 
echo "Running: python3 network_int.py -i $input_directory -o $output_directory_net -t $threshold"
#relative path
python3 "$cwd"/src/network_int.py -i "$input_directory" -o "$output_directory_net" -t "$threshold"

echo "done! check your input folder for alhabridge output and your output directory for the netork :)"


