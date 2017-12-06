#!/usr/bin/perl -w
use strict;


open FH,"CRC_datainfo.txt" or die;
open SH,">run_RSEM_STAR.sh" or die;
while(<FH>){
	chomp;
	my ($id,$seq1,$seq2)=split "\t";
	
	print SH "rsem-calculate-expression --paired-end --no-bam-output -p 4 --forward-prob 0 --star --gzipped-read-file /disk/zhw/CRClncRNA/fastq_data/$id\_1.clean.fq.gz /disk/zhw/CRClncRNA/fastq_data/$id\_2.clean.fq.gz lncRNA/lncRNA_RSEM_index/lncRNA.final.v2.RSEM.star RSEM/$id\_star_lncRNA\n";
	print SH "rsem-calculate-expression --paired-end --no-bam-output -p 4 --forward-prob 0 --star --gzipped-read-file /disk/zhw/CRClncRNA/fastq_data/$id\_1.clean.fq.gz /disk/zhw/CRClncRNA/fastq_data/$id\_2.clean.fq.gz lncRNA/pc_RSEM_index/protein_coding.final.RSEM.star RSEM/$id\_star_pc\n";
	
}

`cat run_RSEM_STAR.sh|xargs -iCMD -P4 bash -c CMD`;
