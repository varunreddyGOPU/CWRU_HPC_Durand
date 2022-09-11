#!/bin/bash

#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem-per-cpu=8G
#SBATCH --time=0-00:30:00 # 30 minutes 
#SBATCH --mail-type=ALL 
#SBATCH --job-name="intan_processing" 
#SBATCH -o slurm.%N.%J.%u.out # STDOUT
#SBATCH -e slurm.%N.%J.%u.err # STDERR


file="filePath.txt"
filePath_=$(cat "$file")
name=$(basename "$filePath_")
genName=${name%.*}

file="reqPath.txt"
reqPath=$(cat "$file")

module load rclone 
rclone copy box:"$filePath_" "$(pwd)"

module load matlab

matlab -nosplash -noFigureWindows -r "addpath('"$reqPath"'); fxn_intan_script_CWRU4_miniboard_2CH_Diff('$name'); clear all; close all; functionSpikeCounting('$name',4.2312,1.2156); clear all; close all; quit"

rm progress.sbatch
rm *.rhd
rm *.mat
rm filePath.txt
rm reqPath.txt

dataDir=$(dirname "$filePath_")

rclone copy "$(pwd)" box:"$dataDir/${name%.*}"

rm *