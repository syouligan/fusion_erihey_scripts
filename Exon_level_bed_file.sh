#!/bin/bash

# Makes exon level bed for all genes in GENCODE v32 (used for alignment).
#######

wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.chr_patch_hapl_scaff.annotation.gff3.gz \
| gunzip --stdout \
| awk '$3 == "exon"' - \
| convert2bed -i gff - \
> /share/ClusterShare/biodata/contrib/scoyou/genomes/human/GRCh38_w_ercc/GRCh38v2_exons.bed
