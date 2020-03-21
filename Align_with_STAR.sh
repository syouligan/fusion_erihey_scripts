#!/bin/bash

## Align using STAR with parameters tuned for STAR-fusion
#######

# Number of cores
ncores=20

# Experimental info (update for single or paired reads, and readlength). Set paired to "yes" if appropriate. NOTE: reads are assumed to be stranded.
paired="yes"
readLength=125
project="fusion_erihey"
genome="GRCh38_w_ercc"

homedir="/share/ScratchGeneral/scoyou/sarah_projects"
data="fastp"
tool="STAR"
results="project_results"
QCDir="logs"
inFileExt=".fastq.gz"

genomeDir="/share/ClusterShare/biodata/contrib/scoyou"
species="genomes/human"
STARDir=$genomeDir/$species/$genome/$tool"_"$readLength

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

# Command to be executed
CommandPE="STAR --genomeDir $STARDir \
--readFilesIn $inFile1 $inFile2 \
--outFileNamePrefix $outDir/$sample. \
--readFilesCommand "gunzip -c" \
--outReadsUnmapped None \
--twopassMode Basic \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 12 \
--chimOutJunctionFormat 1 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--runThreadN $ncores \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1"

CommandSE="STAR --genomeDir $STARDir \
--readFilesIn $inFile1 \
--outReadsUnmapped None \
--twopassMode Basic \
--outFileNamePrefix $outDir/$sample. \
--readFilesCommand "gunzip -c" \
--outSAMstrandField intronMotif \
--outSAMunmapped Within \
--chimSegmentMin 12 \
--chimJunctionOverhangMin 12 \
--chimOutJunctionFormat 1 \
--alignSJDBoverhangMin 10 \
--alignMatesGapMax 100000 \
--alignIntronMax 100000 \
--alignSJstitchMismatchNmax 5 -1 5 5 \
--outSAMattrRGline ID:GRPundef \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimMultimapNmax 20 \
--chimNonchimScoreDropMin 10 \
--runThreadN $ncores \
--peOverlapNbasesMin 12 \
--peOverlapMMp 0.1"

# Submit to queue
if [ $paired = "yes" ]
then
echo "CommandPE "$CommandPE
# $CommandPE
qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -hold_jid 'fastp'$sample -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandPE
else
  echo "CommandSE "$CommandSE
# $CommandSE
qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -hold_jid 'fastp'$sample -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $CommandSE
fi

done