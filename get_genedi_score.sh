
nk=$1
Amp_Del=$2
can=$3
echo $nk,$Amp_Del,$can
perl /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/uniq_mean.pl ${can}.rep_all_${nk}.${Amp_Del}.geneid_precent.gistic >rep_all_${nk}.${Amp_Del}.geneid_precent_gene_num.gistic
python /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/addzero_lncRNA2.py ${can}.rep_all_${nk}.${Amp_Del}.geneid_precent.gistic ${nk}.${Amp_Del}_no_bedtools_lncRNA.bed ${nk} ${Amp_Del}
cat rep_all_${nk}.${Amp_Del}.geneid_precent_gene_num.gistic ${nk}.${Amp_Del}_no_bedtools_lncRNA.bed > all_${nk}.${Amp_Del}.geneid_precent2.gistic
sort -k1,1V -k2,2n all_${nk}.${Amp_Del}.geneid_precent2.gistic>${can}.res_${nk}.${Amp_Del}.geneid_precent_sorted.gistic

rm rep_all_${nk}.${Amp_Del}.geneid_precent_gene_num.gistic
rm ${nk}.${Amp_Del}_no_bedtools_lncRNA.bed
rm all_${nk}.${Amp_Del}.geneid_precent2.gistic
