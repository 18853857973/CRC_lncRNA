#! /bin/sh
genome_fa=/disk/database/human/hg38/Gencode/genome.fa

rsem-prepare-reference -p 25 --star --gtf lncRNA/protein_coding.final.gtf  $genome_fa lncRNA/pc_RSEM_index/protein_coding.final.RSEM.star
