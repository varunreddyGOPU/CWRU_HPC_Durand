#!/bin/bash

#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH --mem-per-cpu=8G 
#SBATCH -o slurm.%N.%J.%u.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.err # STDERR
#SBATCH --mail-type=ALL 
#SBATCH --job-name="format_db" 

module load rclone

dataDir="SPARC_Project_Copy/Recording_Data/Chronic_Implants/VNKA/sub-VN091818"

# Get an array of partial paths to each .rhd file in dataDir
dataFiles=$(rclone lsf -R --absolute --files-only --include "*.rhd" box:"$dataDir")

# Loop through the list of files
for file in $dataFiles
do
    # Get the filename without the path and the file extension
    filename="${file##*/}"
    basename="${filename%.*}"
    parentdir="$(basename "$(dirname "$file")")"
	
    # Check if the file is already in a folder of the same name
	if [ $basename = $parentdir ]; then
        echo "$file already in $basename folder, skipping"
    else
        # Move the file to the new folder
        echo "moving $dataDir$file to $dataDir${file%/*}/$basename"
        rclone move --delete-empty-src-dirs box:"$dataDir$file" box:"$dataDir${file%/*}/$basename"
    fi
done