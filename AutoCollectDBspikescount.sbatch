#!/bin/bash

#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH --mem-per-cpu=8G 
#SBATCH -o slurm.%N.%J.%u.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.err # STDERR
#SBATCH --mail-type=ALL 
#SBATCH --job-name="intan_job" 

dataDir="SPARC_Project_Copy/Recording_Data/Chronic_Implants/VNKA/VN091818"
key="dataB_filt_spikes_count"

module load rclone
dataFiles=$(rclone lsf -R --absolute --files-only --include "*$key.txt" box:"$dataDir")
dataFilesClean=$(tr -s 'txt\n' 'txt;' <<< "$dataFiles")
IFS=';' read -a arr <<< $dataFilesClean

mkdir "/scratch/pbsjobs/$USER/$dataDir/$key/"

for((j=0; j<${#arr[@]}; j++)); do

	filePath="$dataDir${arr[j]}"
	y=$(basename "$filePath")
	#z=${y%.*}
	
	rclone copy box:"$filePath" "/scratch/pbsjobs/$USER/$dataDir/$key/"
	cd "/scratch/pbsjobs/$USER/$dataDir/$key/"
	value=$(cat "$y")
	printf '"%s",%s\n' "$filePath" "$value" >> $key.csv

done

# Write to csv, copy csv, delete rest in scratch
rclone copy "/scratch/pbsjobs/$USER/$dataDir/$key/$key.csv" box:"$dataDir"
