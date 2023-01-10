#!/bin/bash

#SBATCH -N 1 
#SBATCH -c 1 
#SBATCH --mem-per-cpu=8G 
#SBATCH -o slurm.%N.%J.%u.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.err # STDERR
#SBATCH --mail-type=ALL 
#SBATCH --job-name="intan_job" 

dataDir="SPARC_Project_Copy/Recording_Data/Chronic_Implants/VNKA/sub-VN091818"

curDir=$(pwd)
req="$curDir/Requirements_SPARC"
mkdir "/scratch/pbsjobs/$USER"

module load rclone
dataFiles=$(rclone lsf -R --absolute --files-only --include "*.rhd" box:"$dataDir")
dataFilesClean=$(tr -s 'rhd\n' 'rhd;' <<< "$dataFiles")
IFS=';' read -a arr <<< $dataFilesClean

for((j=0; j<${#arr[@]}; j++)); do
	filePath="$dataDir${arr[j]}"
	z="${arr[j]%.*}"
	
	mkdir -p "/scratch/pbsjobs/$USER/$dataDir/$z/"
	cp "$curDir/progress.sbatch" "/scratch/pbsjobs/$USER/$dataDir/$z/"
	
	cd "/scratch/pbsjobs/$USER/$dataDir/$z"
	printf "%s" "$filePath" > filePath.txt
	printf "%s" "$req" > reqPath.txt
	
	sbatch progress.sbatch
done