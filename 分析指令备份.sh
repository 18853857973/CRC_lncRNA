#分别筛选出known和novel,集群
cd /disk/zhw/CRClncRNA/filter/lncRNA
cat lncRNA.final.v2.gtf|grep 'known' >lncRNA.final.v2.known.gtf
cat lncRNA.final.v2.gtf|grep 'novel'>lncRNA.final.v2.novel.gtf

#转为bed文件
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


#第一步筛选：表达量大于2的
#DEsewq_zhw.R
#得到lncRNA gene name,exprmorethan2_lncRNA_genename.txt


#chip-seq流程获取到组蛋白结合位置的read丰度
fastq-dump SRR577511 SRR577511.fastq
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504923.fastq SRR504923.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504924.fastq SRR504924.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504927.fastq SRR504927.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504928.fastq SRR504928.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR577511.fastq SRR577511.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504934.cleaned.fastq SRR504934.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR504935.cleaned.fastq SRR504935.sam&
bowtie -S /disk/database/human/hg38/Gencode/bowtie_index/bowtie_genome /disk/zhw/CRClncRNA/chipseq/rawdata/SRR577543.fastq SRR577543.sam&
samtools view -bS SRR577511.sam >SRR577511.bam
samtools merge SRR504934merge.bam SRR504934.sam SRR504935.sam&
samtools merge SRR504927merge.bam SRR504927.sam SRR504928.sam
samtools merge SRR504923merge.bam SRR504924.sam SRR504923.sam&
samtools merge SRR57754merge.bam SRR577544.sam SRR577543.sam&

mkdir SRR504923
cd SRR504923
macs2 callpeak -t ../SRR504923merge.bam ../SRR504927merge.bam -f BAM -g hs -n test -B -q 0.01
mkdir SRR504934
cd SRR504934
macs2 callpeak -t ../SRR504934merge.bam -c ../SRR504927merge.bam -f BAM -g hs -n test -B -q 0.01
mkdir SRR57754
cd SRR57754
macs2 callpeak -t ../SRR57754merge.bam -c ../SRR577511.bam -f BAM -g hs -n test -B -q 0.01

#运行D:\CRC_lncRNA\chipseq\CRClnaRNA.py，获取到73个lncRNA的转录因子前后2kb位置，正链为start，负链为end，将相同geneid的位置取最大最小，并有链方向
dos2unix exp_more_than2_lncRNA.bed

#lncRNA对应的peak
bedtools intersect -a exp_more_than2_lncRNA.bed -b /disk/zhw/CRClncRNA/chipseq/SRR504923/test_summits.bed  -wb> H3K27ac.bed
bedtools intersect -a exp_more_than2_lncRNA.bed -b /disk/zhw/CRClncRNA/chipseq/SRR57754/test_summits.bed  -wb> H3K4me3.bed
bedtools intersect -a exp_more_than2_lncRNA.bed -b /disk/zhw/CRClncRNA/chipseq/SRR504934/test_summits.bed  -wb> H3K4me1.bed
#得到有转录活性的lncRNA,提取出对应的lncRNA

#差异表达D:\CRC_lncRNA\allcode\TF_diffanay_heatmap_venn.R
#绘制热图、venn、柱状图

#kegg富集分析
#excell 分割出基因symbol,david分析，excell画图

#集群
#转为三列bed
cd /disk/zhw/CRClncRNA/filter/lncRNA
cat num2_lncRNA.final.v2.known.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.known2.bed
cat num2_lncRNA.final.v2.novel.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.novel2.bed
cat protein_coding.final.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>protein_coding.final2.bed

#cat TF_lncRNA.final.v2.known.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.known2.bed
#cat TF_lncRNA.final.v2.novel.gtf|awk '{print $1,$4,$5}'|sort -u |sed 's/ /\t/g'>lncRNA.final.v2.novel2.bed



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
#运行coding_potential_ecdf.R
#node1 root
cd /disk/zhw/CRClncRNA/circos
/disk/soft/circos-0.69-6/bin/circos -conf circoshhj.conf

#有转录活性的基因D:\CRC_lncRNA\diffexp\lncRNA_TF.txt

#运行chip-seq.py 从gtf中筛去除没有转录活性的lncRNA，得到TF_lncRNA.final.v2.novel.gtf、TF_lncRNA.final.v2.known.gtf
#编码能力分析
#把gtf转为fa文件
cd /disk/zhw/CRClncRNA/filter/lncRNA

#dos2unix TF_lncRNA.final.v2.novel.gtf
#dos2unix TF_lncRNA.final.v2.known.gtf
dos2unix num2_lncRNA.final.v2.novel.gtf
dos2unix num2_lncRNA.final.v2.known.gtf


gffread TF_lncRNA.final.v2.known.gtf -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w lncRNA.final.v2.known.fa -W
gffread TF_lncRNA.final.v2.novel.gtf -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w lncRNA.final.v2.novel.fa -W
gffread protein_coding.final.gtf  -o-|gffread -  -g /disk/database/human/hg38/Gencode/genome.fa -w protein_coding.final.fa -W

#CPAT#集群节点2
ssh node2
cd /disk/zhw/CRClncRNA/CPAT
/disk/soft/CPAT-1.2.3/bin/cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/protein_coding.final.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_protein_coding.final
/disk/soft/CPAT-1.2.3/bin/cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_lncRNA.final.v2.novel
/disk/soft/CPAT-1.2.3/bin/cpat.py -g /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.fa -d Human_logitModel.RData -x Human_Hexamer.tsv -o CPAT_lncRNA.final.v2.known

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

#求出外显子密度长度exon_length()
#画ecdf图\CPAT、CNCI、bwtool、画密度曲线图、表达箱线图
#运行 D:\CRC_lncRNA\allcode\CPAT_CNCI_bwtool_lengthdesitiny_boxplot.R


####CNV

#转为bed文件
cd /disk/zhw/CRClncRNA/filter/lncRNA
/disk/soft/bedops/gtf2bed < num2_lncRNA.final.v2.known.gtf|sort-bed - > lncRNA.final.v2.known.bed
/disk/soft/bedops/gtf2bed < num2_lncRNA.final.v2.novel.gtf|sort-bed - > lncRNA.final.v2.novel.bed

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

#运行gtf_bind_same_lncRNA():获取复发样本和未复发样本间以及 肿瘤和正常样本间 差异lncRNA对应的Deseq2 logFoldchange值
cat rec_DESeq2_edgeR_up_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_rec_DESeq2_edgeR_up_gene_for_circos.txt
cat rec_DESeq2_edgeR_down_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_rec_DESeq2_edgeR_down_gene_for_circos.txt
cat normal_DESeq2_edgeR_up_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_normal_DESeq2_edgeR_up_gene_for_circos.txt
cat normal_DESeq2_edgeR_down_gene_for_circos.txt |sort -k 1V,1 -k 2n,2 |uniq >uniq_normal_DESeq2_edgeR_down_gene_for_circos.txt

cat normal_DESeq2_edgeR_up_gene_for_circos.txt | awk {'print $1,$2,$3,-$4'}>uniq_normal_DESeq2_edgeR_up_gene_for_circos.txt
cat normal_DESeq2_edgeR_down_gene_for_circos.txt | awk {'print $1,$2,$3,-$4'}>uniq_normal_DESeq2_edgeR_down_gene_for_circos.txt

#绘制cnv circos图
cd /disk/zhw/CRClncRNA/cnv/circos
/disk/soft/circos-0.69-6/bin/circos -conf circoscnvright.conf


#筛选percent大于25的、绘制柱状图、做overlap
cat novel_scores.gistic.bed|awk -F"\t" '$5>0.25{print $0}'>percentages25novel_scores.gistic.bed
cat known_scores.gistic.bed|awk -F"\t" '$5>0.25{print $0}'>percentages25known_scores.gistic.bed

#novel/known lncRNA 转录本id
cat percentages25novel_scores.gistic.bed |awk '{print $18}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.transcriptid.txt
cat percentages25known_scores.gistic.bed |awk '{print $18}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.transcriptid.txt
#novel/known lncRNA 基因id
cat percentages25novel_scores.gistic.bed |awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.geneid.txt
cat percentages25known_scores.gistic.bed |awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.geneid.txt

cat percentages25novel_scores.gistic.bed|grep "Del"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.Del.geneid.txt
cat percentages25novel_scores.gistic.bed|grep "Amp"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25novel.Amp.geneid.txt
cat percentages25known_scores.gistic.bed|grep "Del"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.Del.geneid.txt
cat percentages25known_scores.gistic.bed|grep "Amp"|awk '{print $16}'|sed 's/";//g'|sed 's/"//g'|uniq >percentages25known.Amp.geneid.txt

cat percentages25novel.Del.geneid.txt percentages25known.Del.geneid.txt>percentages25.Del.geneid.txt
cat percentages25novel.Amp.geneid.txt percentages25known.Amp.geneid.txt>percentages25.Amp.geneid.txt
cat percentages25.Del.geneid.txt percentages25.Amp.geneid.txt >percentages25.geneid.txt

#分别提取出novel和known的基因id和对应的CNV频率（Amp和Del分开）
cat novel_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_novel.Amp.geneid_precent.gistic
cat novel_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_novel.Del.geneid_precent.gistic

cat known_scores.gistic.bed|grep "Amp"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_known.Amp.geneid_precent.gistic
cat known_scores.gistic.bed|grep "Del"|awk '{print $1,$2,$3,$4,$5,$16}'|sed 's/";//g'|sed 's/"//g'|sort|uniq>rep_all_known.Del.geneid_precent.gistic

#所有癌症
ls | awk -F "_" '{print $1}' >../cancer.txt

#单个
#sh pineline.sh
#批量：
#node1 
cd  /disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap
sh foralllcancer13.sh 
sh Del_negetive_num.sh 
#Del改为负数
#known和novel未分开
#运行R文件D:\CRC_lncRNA\cnv\percentCNV\cnv_histone_heatmap.R


#
#cd /disk/zhw/CRClncRNA/chipseq
samtools sort SRR504923merge.bam -o H3K27ac.srt.bam -O BAM
samtools index H3K27ac.srt.bam
samtools sort SRR504934merge.bam -o H3K4me1.srt.bam -O BAM
samtools index H3K4me1.srt.bam

samtools sort SRR57754merge.bam -o H3K4me3.srt.bam -O BAM
samtools index H3K4me3.srt.bam

#每5bp bed文件，比对得到对应的read覆盖度CNClncRNA\getlncRNApositionper5bp>D:\\CRC_lncRNA\\diffexp\\5151_lncRNA_TF_per5bp.bed
dos2unix 5151_lncRNA_TF_per5bp.bed 

bedtools multicov -bams H3K27ac.srt.bam -bed 5151_lncRNA_TF_per5bp.bed >H3K27ac_5151_lncRNA_TF_per5bp.bed
bedtools multicov -bams H3K4me1.srt.bam -bed 5151_lncRNA_TF_per5bp.bed >H3K4me1_5151_lncRNA_TF_per5bp.bed
bedtools multicov -bams H3K4me3.srt.bam -bed 5151_lncRNA_TF_per5bp.bed >H3K4me3_5151_lncRNA_TF_per5bp.bed

samtools view -c  H3K4me1.srt.bam #62073675
samtools view -c  H3K4me3.srt.bam #55858490
samtools view -c  H3K27ac.srt.bam #42343594

#overlap
bedtools intersect -a TF_normal_rec_0.25CNV_lncRNA_position.bed -b H3K27ac.ucsc.bedGraph -wb>H3K27ac.ucsc.bedGraph2.bed
bedtools intersect -a TF_normal_rec_0.25CNV_lncRNA_position.bed -b H3K4me1.ucsc.bedGraph >H3K4me1.ucsc.bedGraph.bed
bedtools intersect -a TF_normal_rec_0.25CNV_lncRNA_position.bed -b H3K4me3.ucsc.bedGraph >H3K4me3.ucsc.bedGraph.bed

#up_down_all.r 补0
#heatmap画热图











#差异lncRNA在13种癌症中的热图，上为上调，下为下调，logFC排序
#将bed文件和gist文件比对获取差异lncRNA的score值
bedtools intersect -a /disk/zhw/CRClncRNA/cnv/19_38scores.sorted.gistic -b intersect_normal_cnv_up_down_gtf.bed -wb | awk '{print $1,$2,$3,$4,$5,$9}'|sed 's/\t/ /g'>intersect_normal_cnv_up_down_scores.gistic.bed
dos2unix intersect_normal_cnv_up_down_scores.gistic.bed
#对应多个去平均值合并
perl /disk/zhw/CRClncRNA/cnv/CNV_gistic_data/uniq_mean.pl  /disk/zhw/CRClncRNA/cnv/differentgene_updown_heatmap/intersect_normal_cnv_up_down_scores.gistic.bed >uniq_intersect_normal_cnv_up_down_scores.gistic.bed

cat /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.sorted.bed /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.sorted.bed >/disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel_knoiwn.sorted.bed
cat /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel.geneid.txt /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.known.geneid.txt > /disk/zhw/CRClncRNA/filter/lncRNA/lncRNA.final.v2.novel_known.geneid.txt
####################################################################################################################
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


#############################################################################################################################



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






#novel
cd /data2/zhw/CRC_lncRNA
less lncRNA.final.v2.novel2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >lncRNA.final.v2.novel3.bed
#/data/software/bwtool-master/bwtool extract bed lncRNA.final.v2.novel3.bed /data/database/hg38/hg38.phastCons100way.bw novel_bwtools
/disk/soft/bwtool/bwtool extract bed lncRNA.final.v2.novel3.bed /disk/database/human/hg38/hg38.phastCons100way.bw novel_bwtools
less -S novel_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >novel.mean.txt

#known
less lncRNA.final.v2.known2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >lncRNA.final.v2.known3.bed
#/data/software/bwtool-master/bwtool extract bed lncRNA.final.v2.known3.bed /data/database/hg38/hg38.phastCons100way.bw known_bwtools
/disk/soft/bwtool/bwtool extract bed lncRNA.final.v2.known3.bed /disk/database/human/hg38/hg38.phastCons100way.bw known_bwtools
less -S known_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >knonw.mean.txt

#protein
less protein_coding.final2.bed | perl -wanle'$b++;$a=$F[2]-$F[1];print $_ unless $a==0' >protein_coding.final3.bed
/data/software/bwtool-master/bwtool extract bed protein_coding.final3.bed /data/database/hg38/hg38.phastCons100way.bw protein_bwtools
less -S protein_bwtools| perl -wanle'my @G=split /,/,$F[4];my $a=0;my $b=0;for(@G){next if /NA/;$a+=$_;$b++}if($b==0){$c="NA"}else{$c=$a/$b}print join "\t", @F[0..3], $c' >protein.mean.txt


