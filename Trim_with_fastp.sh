#!/bin/bash

## Trim files with fastp (https://github.com/OpenGene/fastp)
#######

# Number of cores
ncores=12

# Experimental info (update for single or paired reads)
paired="yes"
homedir="/share/ScratchGeneral/scoyou/sarah_projects"
data="raw_data"
project="fusion_erihey"
tool="fastp"
results="project_results"
QCDir="logs"
inFileExt=".fastq.gz"

# Path for log files
logDir=$homedir/$project/$QCDir
mkdir -p $logDir
echo "logDir $logDir"

# Path to folder containing sample folders
sample_Path=$homedir/$project/$data
echo "sample_Path "$sample_Path

# Make an array containing names of directories containing samples
sample_arr=( $(ls $sample_Path) )
echo "# in samples array ${#sample_arr[@]}"
echo "names in samples array ${sample_arr[@]}"

# Submit command for each sample in array
for sample in ${sample_arr[@]}; do

# Runs loop for only the first sample in array (used for development)
# for sample in ${sample_arr[0]}; do

# Define input directory, define and make output and log directories
inPath=$sample_Path/$sample
echo "inPath $inPath"

outDir=$homedir/$project/$results/$tool/$sample
  mkdir -p $outDir
echo "outDir $outDir"

# Make an array containing names of fastq files within each sample directory
fastq_arr=( $(ls $inPath/*$inFileExt) )
echo "# in fastq array ${#fastq_arr[@]}"
echo "names in fastq array ${fastq_arr[@]}"

# Define input files
inFile1=${fastq_arr[0]}
inFile2=${fastq_arr[1]}
echo "inFile1 $inFile1 inFile2 $inFile2"

outFile1=$outDir/$sample"_trimmed_R1.fastq.gz"
outFile2=$outDir/$sample"_trimmed_R2.fastq.gz"
echo "outFile1 $outFile1 outFile2 $outFile2"

# Command to be executed
CommandSE="fastp --thread $ncores -q 10 -i $inFile1 -o $outFile1"
CommandPE="fastp --thread $ncores -q 10 -i $inFile1 -I $inFile2 -o $outFile1 -O $outFile1"
CommandUMI="fastp --thread $ncores -f12 -t 12 -F 12 -T 12 -q 10 -i $inFile1 -I $inFile2 -o $outFile1 -O $outFile1"

# Submit to queue
if [ $paired = "yes" ]
  then
    if [[ $sample == "SP-"* ]]
    then
    echo "CommandUMI "$CommandUMI
    # $CommandUMI
    qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandUMI
    elif [[ $sample == "K19" ]]
    then
    echo "CommandUMI "$CommandUMI
    # $CommandUMI
    qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandUMI
    else
    echo "CommandPE "$CommandPE
    # $CommandPE
    qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandPE
    fi
    
  else
  echo "CommandSE "$CommandSE
  # $CommandSE
  qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandSE
  fi
  
done