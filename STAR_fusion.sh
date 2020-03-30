#!/bin/bash

## Align using STAR with encode parameters
#######

# Number of cores
ncores=20

# Experimental info (update for single or paired reads, and readlength). Set paired to "yes" if appropriate. NOTE: reads are assumed to be stranded.
project="fusion_erihey"

homedir="/share/ScratchGeneral/scoyou/sarah_projects"
data="fastp"
tool="STAR_fusion"
results="project_results"
QCDir="logs"
inFileExt=".fastq.gz"

# Path for log files
logDir=$homedir/$project/$QCDir
mkdir -p $logDir
echo "logDir $logDir"

# Path to folder containing sample folders
sample_Path=$homedir/$project/$results/$data
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

# Define input files
inFile1=$inPath/$sample"_trimmed_R1.fastq.gz"
inFile2=$inPath/$sample"_trimmed_R2.fastq.gz"
echo "inFile1 $inFile1 inFile2 $inFile2"

# Command to be executed.
Command="/share/ScratchGeneral/scoyou/local/bin/STAR-Fusion-v1.8.0_FULL/STAR-Fusion --genome_lib_dir /share/ScratchGeneral/scoyou/local/resources/GRCh38_gencode_v32_CTAT_lib_Dec062019.plug-n-play/ctat_genome_lib_build_dir \
             --left_fq $inFile1 \
             --right_fq $inFile2 \
             --CPU $ncores \
             --min_junction_reads 5 \
             --min_sum_frags 10 \
             --examine_coding_effect \
             --extract_fusion_reads \
             --FusionInspector validate \
             --denovo_reconstruct \
             --output_dir $outDir"

# Submit to queue
echo "Command "$Command
# $CommandPE
qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -hold_jid 'fastp'$sample -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $Command

done