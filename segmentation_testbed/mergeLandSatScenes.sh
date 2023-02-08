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

# Set the input files to be used in the mosaic
input_files="$input_dir/*.TIF"

# Set the output file name
output_file="$output_dir/$(basename $input_dir)_mosaic.tif"

# Activate the Python environment
source venv/bin/activate

# Create the mosaic
gdal_merge.py -o $output_file $input_files

echo "Mosaic created successfully and saved to $output_file"
