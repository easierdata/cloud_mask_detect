#!/bin/bash

# Check if a directory path was provided as an argument
if [ $# -eq 0 ]
then
    echo "Please provide a directory path as an argument."
    exit 1
fi

# Set the input directory path
input_dir=$1

# Set the output directory path
output_dir=inputs/landsatSceneMosaics

# Check if the output directory exists, if not create it
if [ ! -d "$output_dir" ]
then
    mkdir $output_dir
fi

# Set the output file name
output_file="$output_dir/$(basename $input_dir)_mosaic.tif"

# Initialize an empty list to hold the input files
input_files=""

# Find the bands we want (4, 3, and 2)
for band in 4 3 2
do
    # Build the file name for this band
    file_name="$input_dir/*B$band.TIF"

    # Add the file name to the list of input files
    input_files="$input_files $file_name"
done

# Create the mosaic
echo $input_files
gdal_merge.py -v -o $output_file $input_files

echo "Mosaic created successfully and saved to $output_file"
