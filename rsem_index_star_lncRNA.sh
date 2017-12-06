#! /bin/sh
genome_fa=/disk/database/human/hg38/Gencode/genome.fa

rsem-prepare-reference -p 25 --star --gtf lncRNA/lncRNA.final.v2.gtf  $genome_fa lncRNA/lncRNA_RSEM_index/lncRNA.final.v2.RSEM.star
