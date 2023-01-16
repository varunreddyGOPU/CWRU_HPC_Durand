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

# Get an array of partial paths to each directory in the current parent directory that has prefix "sam-"
directories=$(rclone lsf -R --absolute --dirs-only box:"$dataDir" | grep -E '.*\/sam-.*\/.*\/')
for dir in $directories
do
	# Get the directory name
	dirname=$(basename "$dir")
	# Move the directory up one level
	echo "moving $dataDir$dir to $dataDir""/perf-$dirname"
	rclone move --delete-empty-src-dirs box:"$dataDir$dir" box:"$dataDir""/perf-$dirname"
done