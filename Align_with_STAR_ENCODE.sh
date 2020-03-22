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
parameters="ENCODE"
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

outDir=$homedir/$project/$results/$tool"_"$parameters/$sample
mkdir -p $outDir
echo "outDir $outDir"

# Define input files
inFile1=$inPath/$sample"_trimmed_R1.fastq.gz"
inFile2=$inPath/$sample"_trimmed_R2.fastq.gz"
echo "inFile1 $inFile1 inFile2 $inFile2"

# Command to be executed
# Command to be executed
CommandPE="STAR --genomeDir $STARDir \
--outFileNamePrefix $outDir/$sample. \
--readFilesIn $inFile1 $inFile2 \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--readFilesCommand zcat \
--runThreadN $ncores \
--limitBAMsortRAM 10000000000 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outWigStrand Stranded"

CommandSE="STAR --genomeDir $STARDir \
--outFileNamePrefix $outDir/$sample. \
--readFilesIn $inFile1 \
--outSAMunmapped Within \
--outFilterType BySJout \
--outSAMattributes NH HI AS NM MD \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--sjdbScore 1 \
--readFilesCommand zcat \
--runThreadN $ncores \
--limitBAMsortRAM 10000000000 \
--outSAMtype BAM SortedByCoordinate \
--quantMode TranscriptomeSAM \
--outWigStrand Stranded"

# Submit to queue each sample
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