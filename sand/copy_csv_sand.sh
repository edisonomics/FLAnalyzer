#!/bin/bash

# Define the output directory where files will be copied
output_dir=""

# Loop through each initial directory
for initial_dir in /Change_to_path_destination_directory/*; do
    if [ -d "$initial_dir" ]; then
        # Get the name of the initial directory
        dir_name=$(basename "$initial_dir")

        # Find the 0001.ft1 file inside each initial directory
        file_to_copy="$initial_dir/csv/0001.csv"

        # Check if the file exists
        if [ -f "$file_to_copy" ]; then
            # Copy the file to the output directory with the desired name
            cp "$file_to_copy" "$output_dir/$dir_name.csv"
            echo "Copied $file_to_copy to $output_dir/$dir_name.csv"
        else
            echo "File $file_to_copy not found."
        fi
    fi
done

