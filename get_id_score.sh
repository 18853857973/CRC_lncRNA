can=$1

cat /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/gistic/${can}_scores.gistic |awk '{print "chr"$2,$3,$4,$1,$8}'|sed 's/ /\t/g'>${can}.scores.gistic.bed
/disk/soft/CrossMap-0.2.7/bin/CrossMap.py bed /disk/database/human/convert_chain/hg19ToHg38.over.chain.gz ${can}.scores.gistic.bed ${can}.19_38scores.gistic
sort -k 1V,1 -k 2n,2  ${can}.19_38scores.gistic >${can}.19_38scores.sorted.gistic

bedtools intersect -a ${can}.19_38scores.sorted.gistic -b /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.sorted.bed -wb >${can}.novel_scores.gistic.bed
bedtools intersect -a ${can}.19_38scores.sorted.gistic -b /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.sorted.bed -wb >${can}.known_scores.gistic.bed

cat ${can}.novel_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>${can}.rep_all_novel.Amp.geneid_precent.gistic
cat ${can}.novel_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>${can}.rep_all_novel.Del.geneid_precent.gistic

cat ${can}.known_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>${can}.rep_all_known.Amp.geneid_precent.gistic
cat ${can}.known_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>${can}.rep_all_known.Del.geneid_precent.gistic


sh /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/get_genedi_score.sh novel Del $1
sh /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/get_genedi_score.sh novel Amp $1
sh /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/get_genedi_score.sh known Del $1
sh /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/get_genedi_score.sh known Amp $1

rm ${can}.scores.gistic.bed
rm ${can}.19_38scores.gistic
rm ${can}.novel_scores.gistic.bed
rm ${can}.known_scores.gistic.bed
rm ${can}.19_38scores.sorted.gistic
rm ${can}.rep_all_novel.Amp.geneid_precent.gistic
rm ${can}.rep_all_novel.Del.geneid_precent.gistic
rm ${can}.rep_all_known.Amp.geneid_precent.gistic
rm ${can}.rep_all_known.Del.geneid_precent.gistic
rm ${can}.19_38scores.gistic.unmap