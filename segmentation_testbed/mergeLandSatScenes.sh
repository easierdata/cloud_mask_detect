#!/bin/bash

# Check if a directory path argument is provided
if [ $# -eq 0 ]; then
  echo "Error: No directory path provided."
  echo "Usage: sh mergeLandsatScene.sh <directory path>"
  exit 1
fi

# Set the input directory path
input_dir=$1

# Set the output file name
output_file="merged_landsat_scene.tif"

# Get all the TIF files in the input directory
files=$(ls $input_dir | grep ".TIF")

# Build the gdal_merge.py command
cmd="gdal_merge.py -o $output_file"
for file in $files; do
  cmd="$cmd $input_dir/$file"
done

# Run the gdal_merge.py command
eval $cmd

echo "Merged Landsat scene created: $output_file"
