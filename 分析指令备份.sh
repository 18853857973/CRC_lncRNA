#分别筛选出known和novel,集群
cd /disk/zhw/CRClncRNA/filter/lncRNA
cat lncRNA.final.v2.gtf|grep 'known' >lncRNA.final.v2.known.gtf
cat lncRNA.final.v2.gtf|grep 'novel'>lncRNA.final.v2.novel.gtf

#转为bed文件，zuo
cd /data1/users/zzuo/xulab/CRC_lncRNA/hongwan
gtf2bed < lncRNA.final.v2.known.gtf|sort-bed - > lncRNA.final.v2.known.bed
gtf2bed < lncRNA.final.v2.novel.gtf|sort-bed - > lncRNA.final.v2.novel.bed

#RSEM 集群
cd /disk/zhw/CRClncRNA/fastq_data/first
ls *_1.clean.fq.gz|sed 's/_1.clean.fq.gz//g'>SRR.txt

/disk/soft/CPAT-1.2.3/bin/cpat.py -r /disk/database/human/hg38/Gencode/genome.fa -g /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.bed
#fastp clean data
cd /disk/zhw/CRClncRNA/fastq_data/sec
fastp -i D6_HKNNLCCXX_L4_1.fq.gz -o D6_HKNNLCCXX_L4_1.clean.fq.gz -I D6_HKNNLCCXX_L4_2.fq.gz -O D6_HKNNLCCXX_L4_2.clean.fq.gz -w 8

#集群
#转为三列bed
cd /disk/zhw/CRClncRNA/filter/lncRNA
cat lncRNA.final.v2.known.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.known2.bed
cat lncRNA.final.v2.novel.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.novel2.bed
cat protein_coding.final.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>protein_coding.final2.bed

#34fuwuq
#novel
cd /data2/zhw/CRC_lncRNA
less lncRNA.final.v2.novel2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >lncRNA.final.v2.novel3.bed
/data/software/bwtool-master/bwtool extract bed lncRNA.final.v2.novel3.bed /data/database/hg38/hg38.phastCons100way.bw novel_bwtools
less -S novel_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >novel.mean.txt

#known
less lncRNA.final.v2.known2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >lncRNA.final.v2.known3.bed
/data/software/bwtool-master/bwtool extract bed lncRNA.final.v2.known3.bed /data/database/hg38/hg38.phastCons100way.bw known_bwtools
less -S known_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >knonw.mean.txt

#protein
less protein_coding.final2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >protein_coding.final3.bed
/data/software/bwtool-master/bwtool extract bed protein_coding.final3.bed /data/database/hg38/hg38.phastCons100way.bw protein_bwtools
less -S protein_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >protein.mean.txt


#合并（集群）
cd /disk/zhw/CRClncRNA/filter/RSEM

python pcRNA_RSEM_merge.py 5 pcRNA.rsem.count
python pcRNA_RSEM_merge.py 6 pcRNA.rsem.TPM
python pcRNA_RSEM_merge.py 7 pcRNA.rsem.FPKM

python lncRNA_RSEM_merge.py 5 lncRNA.rsem.count
python lncRNA_RSEM_merge.py 6 lncRNA.rsem.TPM
python lncRNA_RSEM_merge.py 7 lncRNA.rsem.FPKM

cd /disk/zhw/CRClncRNA/filter/lncRNA
cat lncRNA.final.v2.gtf|awk '{print $1,$2,$4,$5,$10}'|sed 's/"//g'|sed 's/;//g'|sed 's/ /\t/g' >lncRNA.final.v2.cp.2.txt



#画circos图
#分类
#运行D:\CRC_lncRNA目录的distinguish_known_novel.R
#node1 root

#test
head -1000 hs_exp_data/hs_normal_know_gtf_data.txt| awk {'print $1,$2,$3,$4'}>hs_normal_know_gtf_data3.txt
head -1000 hs_exp_data/hs_normal_novel_gtf_data.txt| awk {'print $1,$2,$3,$4'}>hs_normal_novel_gtf_data3.txt
head -1000 hs_exp_data/hs_unrecurrence_know_gtf_data.txt| awk {'print $1,$2,$3,$4'}>hs_unrecurrence_know_gtf_data3.txt
head -1000 hs_exp_data/hs_unrecurrence_novel_gtf_data.txt | awk {'print $1,$2,$3,$4'}>hs_unrecurrence_novel_gtf_data3.txt
head -1000 hs_exp_data/hs_recurrence_novel_gtf_data.txt| awk {'print $1,$2,$3,$4'}>hs_recurrence_novel_gtf_data3.txt
head -1000 hs_exp_data/hs_recurrence_know_gtf_data.txt| awk {'print $1,$2,$3,$4'}>hs_recurrence_know_gtf_data3.txt

/disk/soft/circos-0.69-6/bin/circos -conf circoshhj_test.conf

/disk/soft/circos-0.69-6/bin/circos -conf circoshhj.conf





#编码能力分析
#集群节点2
#把gtf转为fa文件
cd /disk/zhw/CRClncRNA/filter/lncRNA

gffread lncRNA.final.v2.known.gtf -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w lncRNA.final.v2.known.fa -W
gffread lncRNA.final.v2.novel.gtf -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w lncRNA.final.v2.novel.fa -W
gffread protein_coding.final.gtf  -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w protein_coding.final.fa -W

#CPAT
cd /disk/zhw/CRClncRNA/CPAT
cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/protein_coding.final.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_protein_coding.final
cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_lncRNA.final.v2.novel
cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_lncRNA.final.v2.known

#CNCI
cd /disk/zhw/CRClncRNA/CNCI
python /disk/soft/CNCI-master/CNCI.py -f /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.fa -o CNCI_lncRNA.final.v2.novel -m ve -p 4
python /disk/soft/CNCI-master/CNCI.py -f /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.fa -o CNCI_lncRNA.final.v2.known -m ve -p 4
python /disk/soft/CNCI-master/CNCI.py -f /disk/zhw/CRClncRNA/filter/lncRNA/protein_coding.final.fa -o CNCI_protein_coding.final -m ve -p 4

#bwtool
#去除长度为0的
cd /disk/zhw/CRClncRNA/bwtools
less /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' > lncRNA.final.v2.known_filt.bed
less /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' > lncRNA.final.v2.novel_filt.bed
less /disk/zhw/CRClncRNA/filter/lncRNA/protein_coding.final2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' > protein_coding.final_filt.bed

/disk/soft/bwtool/bwtool extract bed  lncRNA.final.v2.known_filt.bed /disk/database/human/hg38/hg38.phastCons100way.bw  lncRNA.final.v2.known
/disk/soft/bwtool/bwtool extract bed  lncRNA.final.v2.novel_filt.bed /disk/database/human/hg38/hg38.phastCons100way.bw  lncRNA.final.v2.novel
/disk/soft/bwtool/bwtool extract bed  protein_coding.final_filt.bed /disk/database/human/hg38/hg38.phastCons100way.bw  protein_coding.final

#求平均值
less -S protein_coding.final | perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' > bwtool_protein_coding.final_mean.txt
less -S lncRNA.final.v2.known | perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' > bwtool_lncRNA.final.v2.known_mean.txt
less -S lncRNA.final.v2.novel | perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' > bwtool_lncRNA.final.v2.novel_mean.txt


#找出lncRNA附件的蛋白质基因名字
cd /disk/zhw/CRClncRNA/filter/lncRN
cat lncRNA.final.v2.gtf |awk '{print $10}'|sed 's/"//g'|sed 's/;//g'|sed 's/ /\t/g'|uniq > lncRNA.final.v2.gtf.genename.txt


#画ecdf图
#运行 D:\CRC_lncRNA\coding_potential\coding_potential_ecdf.R

#画密度曲线图
D:\CRC_lncRNA\filter\lncRNA\\exon_length.py
D:\CRC_lncRNA\filter\lncRNA\\length_desitny.R


#差异表达


#表达箱线图
D:\CRC_lncRNA\filter\RSEM_expression
boxplot.R

#kegg富集分析
#excell 分割出基因symbol,david分析，excell画图





####CNV

cd /disk/zhw/CRClncRNA/filter/lncRNA
sort -k 1V,1 -k 2n,2 lncRNA.final.v2.novel.bed >lncRNA.final.v2.novel.sorted.bed
sort -k 1V,1 -k 2n,2 lncRNA.final.v2.known.bed >lncRNA.final.v2.known.sorted.bed


#所有novel/known lncRNA 基因id
less lncRNA.final.v2.novel.sorted.bed | perl -wanle'next if $F[0]=~/_/;print $_' >lncRNA.final.v2.novel.sorted2.bed
less lncRNA.final.v2.known.sorted.bed | perl -wanle'next if $F[0]=~/_/;print $_' >lncRNA.final.v2.known.sorted2.bed
mv lncRNA.final.v2.novel.sorted2.bed lncRNA.final.v2.novel.sorted.bed
mv lncRNA.final.v2.known.sorted2.bed lncRNA.final.v2.known.sorted.bed

cat lncRNA.final.v2.novel.sorted.bed|awk '{print $11}'|sed 's/";//g'|sed 's/"//g'|uniq >lncRNA.final.v2.novel.geneid2.txt  
cat lncRNA.final.v2.known.sorted.bed|awk '{print $11}'|sed 's/";//g'|sed 's/"//g'|uniq >lncRNA.final.v2.known.geneid2.txt
sort lncRNA.final.v2.novel.geneid2.txt |uniq >lncRNA.final.v2.novel.geneid.txt
sort lncRNA.final.v2.known.geneid2.txt |uniq >lncRNA.final.v2.known.geneid.txt

#所有novel/known lncRNA 转录本id
cat lncRNA.final.v2.novel.sorted2.bed|awk '{print $13}'|sed 's/";//g'|sed 's/"//g'|uniq >lncRNA.final.v2.novel.transcriptid.txt
cat lncRNA.final.v2.known.sorted2.bed|awk '{print $13}'|sed 's/";//g'|sed 's/"//g'|uniq >lncRNA.final.v2.known.transcriptid.txt


cd /disk/zhw/CRClncRNA/cnv
cat scores.gistic |awk '{print "chr"$2,$3,$4,$1,$8}'|sed 's/ /\t/g'>scores.gistic.bed
/disk/soft/CrossMap-0.2.7/bin/CrossMap.py bed /disk/database/human/convert_chain/hg19ToHg38.over.chain.gz scores.gistic.bed 19_38scores.gistic
sort -k 1V,1 -k 2n,2  19_38scores.gistic >19_38scores.sorted.gistic

#分为Amp Del
cat 19_38scores.sorted.gistic|grep "Amp"|awk '{print $1,$2,$3,$5}' |sed 's/chr/hs/g'>Amp_19_38scores.sorted.gistic
cat 19_38scores.sorted.gistic|grep "Del"|awk '{print $1,$2,$3,$5}' |sed 's/chr/hs/g'>Del_19_38scores.sorted.gistic



bedtools intersect -a 19_38scores.sorted.gistic -b /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.sorted.bed -wb >novel_scores.gistic.bed
bedtools intersect -a 19_38scores.sorted.gistic -b /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.sorted.bed -wb >known_scores.gistic.bed

#circos图使用，不要了
#cat novel_scores.gistic.bed |grep "Amp"|awk {'print $1,$7,$8,$5*100'}|sed 's/chr/hs/g'>Amp_novel_scores.gistic.bed
#cat novel_scores.gistic.bed |grep "Del"|awk {'print $1,$7,$8,$5*100'}|sed 's/chr/hs/g'>Del_novel_scores.gistic.bed
#cat known_scores.gistic.bed |grep "Amp"|awk {'print $1,$7,$8,$5*100'}|sed 's/chr/hs/g'>Amp_known_scores.gistic.bed
#cat known_scores.gistic.bed |grep "Del"|awk {'print $1,$7,$8,$5*100'}|sed 's/chr/hs/g'>Del_known_scores.gistic.bed

#取位点中点和+1
#le novel_scores.gistic.bed |grep "Amp"| perl -wanle'$mean = int($F[6]+$F[7]);$mean2 = $mean+1;$haha=$F[4]*100;$F[0]=~s/chr/hs/;print "$F[0] $mean $mean2 $haha"'|uniq>Amp_novel_scores.gistic.bed
#le novel_scores.gistic.bed |grep "Del"| perl -wanle'$mean = int($F[6]+$F[7]);$mean2 = $mean+1;$haha=$F[4]*100;$F[0]=~s/chr/hs/;print "$F[0] $mean $mean2 $haha"'|uniq>Del_novel_scores.gistic.bed
#le known_scores.gistic.bed |grep "Amp"| perl -wanle'$mean = int($F[6]+$F[7]);$mean2 = $mean+1;$haha=$F[4]*100;$F[0]=~s/chr/hs/;print "$F[0] $mean $mean2 $haha"'|uniq>Amp_known_scores.gistic.bed
#le known_scores.gistic.bed |grep "Del"| perl -wanle'$mean = int($F[6]+$F[7]);$mean2 = $mean+1;$haha=$F[4]*100;$F[0]=~s/chr/hs/;print "$F[0] $mean $mean2 $haha"'|uniq>Del_known_scores.gistic.bed



#绘制circos图
cd /disk/zhw/CRClncRNA/cnv/circos
/disk/soft/circos-0.69-6/bin/circos -conf circoscnv.conf

#test:
#/disk/soft/circos-0.69-6/bin/circos -conf circos2.conf 
#head -1000 novel_scores.gistic.bed |grep "Amp"|awk {'print $1,$2,$3,$5*100'}|sed 's/chr/hs/g'>Amp_novel_scores.gistic2.bed
#head -1000 known_scores.gistic.bed |grep "Amp"|awk {'print $1,$2,$3,$5*100'}|sed 's/chr/hs/g'>Amp_known_scores.gistic2.bed
#head -1000 novel_scores.gistic.bed |grep "Del"|awk {'print $1,$2,$3,$5*100'}|sed 's/chr/hs/g'>Del_novel_scores.gistic2.bed
#head -1000 known_scores.gistic.bed |grep "Del"|awk {'print $1,$2,$3,$5*100'}|sed 's/chr/hs/g'>Del_known_scores.gistic2.bed


#筛选percent大于25的、绘制柱状图、做overlap
cat novel_scores.gistic.bed|awk -F"\t" '$5>0.25{print $0}'>percentages25novel_scores.gistic.bed
cat known_scores.gistic.bed|awk -F"\t" '$5>0.25{print $0}'>percentages25known_scores.gistic.bed

#novel/known lncRNA 转录本id
cat percentages25novel_scores.gistic.bed |awk '{print $18}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.transcriptid.txt
cat percentages25known_scores.gistic.bed |awk '{print $18}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.transcriptid.txt
#novel/known lncRNA 基因id
cat percentages25novel_scores.gistic.bed |awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.geneid.txt
cat percentages25known_scores.gistic.bed |awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.geneid.txt


#分别提取出novel和known的基因id和对应的CNV频率（Amp和Del分开）
cat novel_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_novel.Amp.geneid_precent.gistic
cat novel_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_novel.Del.geneid_precent.gistic

cat known_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_known.Amp.geneid_precent.gistic
cat known_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_known.Del.geneid_precent.gistic


#python获取到没有比对部分，合并两个文件#获取到有cnv的 取到没有CNV变异的复制为0 合并 #排序
#cd /disk/zhw/CRClncRNA/cnv/heatmap
#perl uniq_mean.pl ../rep_all_novel.Amp.geneid_precent.gistic >rep_all_novel.Amp.geneid_precent_gene_num.gistic
#python addzero_lncRNA2.py /disk/zhw/CRClncRNA/cnv/rep_all_novel.Amp.geneid_precent.gistic novel_Amp_no_bedtools_lncRNA.bed
#cat rep_all_novel.Amp.geneid_precent_gene_num.gistic novel_Amp_no_bedtools_lncRNA.bed > all_novel.Amp.geneid_precent2.gistic
#sort -k1,1V -k2,2n all_novel.Amp.geneid_precent2.gistic>all_novel.Amp.geneid_precent2_dorted.gistic3

#Del
#perl uniq_mean.pl ../rep_all_novel.Del.geneid_precent.gistic >rep_all_novel.Del.geneid_precent_gene_num.gistic
#python addzero_lncRNA2.py /disk/zhw/CRClncRNA/cnv/rep_all_novel.Del.geneid_precent.gistic novel_Del_no_bedtools_lncRNA.bed
#cat rep_all_novel.Del.geneid_precent_gene_num.gistic novel_Del_no_bedtools_lncRNA.bed > all_novel.Del.geneid_precent2.gistic
#sort -k1,1V -k2,2n all_novel.Del.geneid_precent2.gistic>all_novel.Del.geneid_precent2_dorted.gistic3

#perl uniq_mean.pl all_novel.Del.geneid_precent2.gistic >all_novel.Del.geneid_precent_gene_num.gistic
#perl uniq_mean.pl all_known.Amp.geneid_precent2.gistic >all_known.Amp.geneid_precent_gene_num.gistic
#perl uniq_mean.pl all_known.Del.geneid_precent2.gistic >all_known.Del.geneid_precent_gene_num.gistic

ls | awk -F "_" '{print $1}' >../cancer.txt

#生成矩阵
sh foralllcancer.sh
#把Del改为负数
sh Del_negetive_num.sh

#合并known\novel
sh bind_knownnovel.sh
运行R文件D:\CRC_lncRNA\cnv\percentCNV\allcancer_chro_heatmap.R

cat percentages25novel_scores.gistic.bed|grep "Del"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.Del.geneid.txt
cat percentages25novel_scores.gistic.bed|grep "Amp"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.Amp.geneid.txt

cat percentages25known_scores.gistic.bed|grep "Del"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.Del.geneid.txt
cat percentages25known_scores.gistic.bed|grep "Amp"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.Amp.geneid.txt

cat percentages25novel.Del.geneid.txt percentages25known.Del.geneid.txt>percentages25.Del.geneid.txt
cat percentages25novel.Amp.geneid.txt percentages25known.Amp.geneid.txt>percentages25.Amp.geneid.txt

#运行D:\CRC_lncRNA\diffexp\gtf_bind_same_lncNRA.py
#先将D:\CRC_lncRNA\filter\lncRNA\lncRNA.final.v2.gtf文件only_min_max_position_lncRNA.final.v2.gtf,只有最大最小的位置，将一个lncRNA多个位置合并为一个
#运行：D:\CRC_lncRNA\diffexp\different_gene_position_log2FoldChange_for_circos.py

#cat rec_DESeq2_edgeR_up_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_rec_DESeq2_edgeR_up_gene_for_circos.txt
#cat rec_DESeq2_edgeR_down_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_rec_DESeq2_edgeR_down_gene_for_circos.txt
#cat normal_DESeq2_edgeR_up_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_normal_DESeq2_edgeR_up_gene_for_circos.txt
#cat normal_DESeq2_edgeR_down_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_normal_DESeq2_edgeR_down_gene_for_circos.txt

/disk/soft/circos-0.69-6/bin/circos -conf circoscnvright.conf 

#test
head -100 uniq_rec_DESeq2_edgeR_up_gene_for_circos.txt>uniq_rec_DESeq2_edgeR_up_gene_for_circos2.txt
head -100 uniq_rec_DESeq2_edgeR_down_gene_for_circos.txt> uniq_rec_DESeq2_edgeR_down_gene_for_circos2.txt
head -100 uniq_normal_DESeq2_edgeR_up_gene_for_circos.txt>uniq_normal_DESeq2_edgeR_up_gene_for_circos2.txt
head -100 uniq_normal_DESeq2_edgeR_down_gene_for_circos.txt >uniq_normal_DESeq2_edgeR_down_gene_for_circos2.txt
head -10000 ../Del_19_38scores.sorted.gistic >../Del_19_38scores.sorted.gistic2
head -10000 ../Amp_19_38scores.sorted.gistic >../Amp_19_38scores.sorted.gistic2
/disk/soft/circos-0.69-6/bin/circos -conf circoscnvtest.conf


#差异lncRNA在13种癌症中的热图，上为上调，下为下调，logFC排序
#将bed 文件和gist文件比对获取差异lncRNA的score值
bedtools intersect -a /disk/zhw/CRClncRNA/cnv/19_38scores.sorted.gistic -b intersect_normal_cnv_up_down_gtf.bed -wb | awk '{print $1,$2,$3,$4,$5,$9}'|sed 's/\t/ /g'>intersect_normal_cnv_up_down_scores.gistic.bed
dos2unix intersect_normal_cnv_up_down_scores.gistic.bed
#对应多个去平均值合并
perl /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/uniq_mean.pl  /disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap/intersect_normal_cnv_up_down_scores.gistic.bed >uniq_intersect_normal_cnv_up_down_scores.gistic.bed


cat /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.sorted.bed /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.sorted.bed >/disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel_knoiwn.sorted.bed
cat /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.geneid.txt /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.geneid.txt > /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel_known.geneid.txt

#############################
#单个流程
cat /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/gistic/BRCA_scores.gistic|awk '{print "chr"$2,$3,$4,$1,$8}'|sed 's/ /\t/g'>BRCA_scores.gistic.bed
/disk/soft/CrossMap-0.2.7/bin/CrossMap.py bed /disk/database/human/convert_chain/hg19ToHg38.over.chain.gz BRCA_scores.gistic.bed BRCA.19_38scores.gistic
sort -k 1V,1 -k 2n,2  BRCA.19_38scores.gistic >BRCA.19_38scores.sorted.gistic
bedtools intersect -a  BRCA.19_38scores.sorted.gistic  -b  /disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap/intersect_normal_cnv_up_down_gtf.bed -wb | awk '{print $1,$2,$3,$4,$5,$9}'|sed 's/\t/ /g'>BRCA_scores.gistic.bed
dos2unix BRCA_scores.gistic.bed
cat BRCA_scores.gistic.bed |grep "Amp"|sed 's/";//g'|sed 's/"//g'|sort|uniq>BRCA.rep_all.Amp.geneid_precent.gistic
cat BRCA_scores.gistic.bed |grep "Del"|sed 's/";//g'|sed 's/"//g'|sort|uniq>BRCA.rep_all.Del.geneid_precent.gistic

perl /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/uniq_mean.pl  BRCA.rep_all.Amp.geneid_precent.gistic >BRCA.rep_all_Amp.geneid_precent_gene_num.gistic
perl /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/uniq_mean.pl  BRCA.rep_all.Del.geneid_precent.gistic >BRCA.rep_all_Del.geneid_precent_gene_num.gistic

python addzero_lncRNA2.py BRCA.rep_all_Amp.geneid_precent_gene_num.gistic  BRCA.Amp_no_bedtools_lncRNA.bed
python addzero_lncRNA2.py BRCA.rep_all_Del.geneid_precent_gene_num.gistic  BRCA.Del_no_bedtools_lncRNA.bed

cat UVM.rep_all_Amp.geneid_precent_gene_num.gistic UVM.Amp_no_bedtools_lncRNA.bed>UVM.res_novel_known.Amp.geneid_precent_sorted.gistic
cat UVM.rep_all_Del.geneid_precent_gene_num.gistic UVM.Del_no_bedtools_lncRNA.bed>UVM.res_novel_known.Del.geneid_precent_sorted.gistic

cat UCEC.res_novel_known.Del.geneid_precent_sorted.gistic | awk '{print $1"\t"$2"\t"$3"\t-"$4}'> UCEC.res_novel_known.Del.geneid_precent_sorted2.gistic

##############################
#单个
sh pineline.sh
#批量： 
cd  /disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap
sh foralllcancer13.sh 
sh Del_negetive_num.sh 
#Del改为负数

#chip-seq
fastq-dump SRR577511 SRR577511.fastq

bowtie /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome SRR577511.fastq SRR577511.sam




















