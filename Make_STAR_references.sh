#!/bin/bash

## Make a STAR genome reference
#---------------------------------------------------------------

# Number of cores
ncores=32

#Make indexes for all standard lengths
for i in 125
do

# Experimental info (update for genome)
tool="STAR"
readLength=$i
overHang=$(expr $readLength - 1)
QCDir="logs"
species="genomes/human"
genome="GRCh38_w_ercc"
sequence="GRCh38.p13.genome_ERCC92.fa"
annotation="gencode.v32.chr_patch_hapl_scaff.annotation_ERCC92.gtf"

# Make directory structure
homedir="/share/ClusterShare/biodata/contrib/scoyou"

inputDir=$homedir/$species/$genome
echo "inPath "$inputDir

outDir=$homedir/$species/$genome/$tool"_"$readLength
mkdir -p $outDir
echo "outPath "$outDir

logDir=$outDir/$QCDir
mkdir -p $logDir
echo "logDir "$logDir

# Build star index for genome
starCommand="STAR --runMode genomeGenerate --genomeDir $outDir --genomeFastaFiles $inputDir/$sequence --outFileNamePrefix $outDir/ --runThreadN $ncores --sjdbGTFfile $inputDir/$annotation --sjdbOverhang $overHang"
echo "starCommand "$starCommand

# Submit to queue
qsub -P OsteoporosisandTranslationalResearch -N $tool$genome$i -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $starCommand

done