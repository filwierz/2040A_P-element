diagnostic\_SNPs
================
Filip Wierzbicki
4/5/2023

``` bash
cd /Volumes/Temp2/filip/2040A/data/strains/Dsim/trimmed
nohup sh -c 'for i in *.fq.gz;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dsim_3scg.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/${n}.sort.bam;done' &
cd /Volumes/Temp2/filip/2040A/data/strains/Dmel/trimmed
nohup sh -c 'for i in *.fq.gz;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dmel_3scg.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/${n}.sort.bam;done' &

##cross-map
cd /Volumes/Temp2/filip/2040A/data/strains/Dsim/trimmed
nohup sh -c 'for i in *.fq.gz;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dmel_3scg.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map/${n}.sort.bam;done' &
cd /Volumes/Temp2/filip/2040A/data/strains/Dmel/trimmed
nohup sh -c 'for i in *.fq.gz;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dsim_3scg.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map/${n}.sort.bam;done' &


##need to include population data because to many segregating SNPs:
###from https://doi.org/10.1371/journal.pgen.1005406
#ERR557065: http://ftp.sra.ebi.ac.uk/vol1/run/ERR557/ERR557065/DsimVPG-PCRfreec-rep2.sort.bam
#ERR557050: http://ftp.sra.ebi.ac.uk/vol1/run/ERR557/ERR557050/Dmel-PCRfree-rep1.sort.bam   
cd /Volumes/Temp2/filip/2040A/data/populations
samtools view Dmel-PCRfree-rep1.sort.bam|awk '{print "@"$1;print $10;print "+"$1;print $11}' > fastq/Dmel.fastq
samtools view DsimVPG-PCRfreec-rep2.sort.bam|awk '{print "@"$1;print $10;print "+"$1;print $11}' > fastq/Dsim.fastq
cd fastq
bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dsim_3scg.fasta Dsim.fastq|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/Dsim.sort.bam
bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dmel_3scg.fasta Dmel.fastq|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/Dmel.sort.bam

##cross-map
bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dmel_3scg.fasta Dsim.fastq|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map/Dsim.sort.bam
bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/scg-only/Dsim_3scg.fasta Dmel.fastq|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map/Dmel.sort.bam

#dmel-ref:
cd /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map
nohup samtools mpileup -q 20 ../Canton-S.sort.bam ../Iso1.sort.bam ../Oregon-R.sort.bam ../Dmel.sort.bam Mod6.sort.bam w501.sort.bam wXD1.sort.bam Dsim.sort.bam -o /Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/withPop_Dmel-scg.mpileup &
cd /Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/
java -ea -Xmx7g -jar /Volumes/Temp3/filip/programs/popoolation2_1201/mpileup2sync.jar --input withPop_Dmel-scg.mpileup --output withPop_Dmel-scg.sync --fastq-type sanger --min-qual 20 --threads 8

python /Volumes/Temp2/filip/2040A/programs/2040A_P-element/helper-scripts/diagnosticSNPs-finder.py --sync withPop_Dmel-scg.sync --selfIDs 1,2,3,4 --nonselfIDs 5,6,7,8 --min-cov 5 > withPop_Dmel-scg.SNPs

#dsim-ref:
cd /Volumes/Temp2/filip/2040A/map/3scg-only/cross-map
nohup samtools mpileup -q 20 ../Mod6.sort.bam ../w501.sort.bam ../wXD1.sort.bam ../Dsim.sort.bam Canton-S.sort.bam Iso1.sort.bam Oregon-R.sort.bam Dmel.sort.bam -o /Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/withPop_Dsim-scg.mpileup &
cd /Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/
java -ea -Xmx7g -jar /Volumes/Temp3/filip/programs/popoolation2_1201/mpileup2sync.jar --input withPop_Dsim-scg.mpileup --output withPop_Dsim-scg.sync --fastq-type sanger --min-qual 20 --threads 8

python /Volumes/Temp2/filip/2040A/programs/2040A_P-element/helper-scripts/diagnosticSNPs-finder.py --sync withPop_Dsim-scg.sync --selfIDs 1,2,3,4 --nonselfIDs 5,6,7,8 --min-cov 5 > withPop_Dsim-scg.SNPs
```

Dmel genes

``` r
m<-read.table("/Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/withPop_Dmel-scg.SNPs")
names(m)<-c("gene","position","Dmel-variant","Dsim-variant")
knitr::kable(m)
```

| gene        | position | Dmel-variant | Dsim-variant |
| :---------- | -------: | :----------- | :----------- |
| Dmel\_rhi   |     1252 | A            | T            |
| Dmel\_rhi   |     1560 | G            | C            |
| Dmel\_rhi   |     1986 | C            | A            |
| Dmel\_rhi   |     1999 | T            | A            |
| Dmel\_rhi   |     2249 | T            | A            |
| Dmel\_rhi   |     2388 | T            | C            |
| Dmel\_rhi   |     2456 | A            | G            |
| Dmel\_rhi   |     2460 | A            | C            |
| Dmel\_rhi   |     2484 | G            | T            |
| Dmel\_rhi   |     2515 | G            | A            |
| Dmel\_rhi   |     2558 | T            | C            |
| Dmel\_rhi   |     2560 | C            | T            |
| Dmel\_rhi   |     2746 | G            | T            |
| Dmel\_rhi   |     3482 | A            | G            |
| Dmel\_rhi   |     3509 | T            | C            |
| Dmel\_rhi   |     3897 | T            | C            |
| Dmel\_rhi   |     3978 | A            | G            |
| Dmel\_rhi   |     4031 | A            | G            |
| Dmel\_rhi   |     4129 | G            | C            |
| Dmel\_rhi   |     4131 | C            | T            |
| Dmel\_rhi   |     4231 | A            | G            |
| Dmel\_rhi   |     4239 | C            | A            |
| Dmel\_rhi   |     4391 | T            | A            |
| Dmel\_rhi   |     4509 | T            | C            |
| Dmel\_rhi   |     4510 | T            | G            |
| Dmel\_rhi   |     4550 | T            | G            |
| Dmel\_rhi   |     4595 | C            | T            |
| Dmel\_rhi   |     4599 | A            | T            |
| Dmel\_rhi   |     4603 | C            | A            |
| Dmel\_rhi   |     4667 | T            | C            |
| Dmel\_rhi   |     4732 | T            | C            |
| Dmel\_rhi   |     4742 | T            | C            |
| Dmel\_rhi   |     4760 | T            | G            |
| Dmel\_rhi   |     4987 | A            | G            |
| Dmel\_rhi   |     4989 | A            | T            |
| Dmel\_rhi   |     5186 | A            | T            |
| Dmel\_rhi   |     5193 | C            | T            |
| Dmel\_rhi   |     5380 | G            | A            |
| Dmel\_rhi   |     5487 | G            | A            |
| Dmel\_rhi   |     5496 | A            | G            |
| Dmel\_rhi   |     5760 | C            | A            |
| Dmel\_rhi   |     5788 | G            | A            |
| Dmel\_rhi   |     5807 | T            | C            |
| Dmel\_rhi   |     5811 | A            | G            |
| Dmel\_rhi   |     5822 | T            | G            |
| Dmel\_rhi   |     5839 | C            | T            |
| Dmel\_rhi   |     5850 | G            | A            |
| Dmel\_rhi   |     6161 | T            | C            |
| Dmel\_rhi   |     6191 | G            | A            |
| Dmel\_rhi   |     6196 | A            | T            |
| Dmel\_rhi   |     6476 | C            | T            |
| Dmel\_rhi   |     6490 | G            | A            |
| Dmel\_rhi   |     6492 | A            | G            |
| Dmel\_rhi   |     6493 | T            | C            |
| Dmel\_rhi   |     6678 | C            | T            |
| Dmel\_rhi   |     6685 | A            | T            |
| Dmel\_rhi   |     6688 | A            | G            |
| Dmel\_rhi   |     6713 | A            | G            |
| Dmel\_rhi   |     6731 | A            | G            |
| Dmel\_rhi   |     6740 | C            | A            |
| Dmel\_rhi   |     6936 | A            | C            |
| Dmel\_rhi   |     6969 | C            | T            |
| Dmel\_rhi   |     6974 | G            | T            |
| Dmel\_rhi   |     6993 | T            | A            |
| Dmel\_rhi   |     6999 | C            | G            |
| Dmel\_rhi   |     7001 | C            | A            |
| Dmel\_rhi   |     7059 | C            | T            |
| Dmel\_rhi   |     7060 | A            | C            |
| Dmel\_rhi   |     7071 | T            | C            |
| Dmel\_rhi   |     7357 | A            | T            |
| Dmel\_rhi   |     7368 | T            | C            |
| Dmel\_rhi   |     7369 | T            | A            |
| Dmel\_rhi   |     7414 | G            | C            |
| Dmel\_rhi   |     7420 | G            | A            |
| Dmel\_rhi   |     7540 | T            | C            |
| Dmel\_rhi   |     7570 | A            | G            |
| Dmel\_rhi   |     7576 | C            | T            |
| Dmel\_rhi   |     7846 | G            | A            |
| Dmel\_rhi   |     8048 | A            | T            |
| Dmel\_rhi   |     8052 | G            | T            |
| Dmel\_rhi   |     8056 | A            | G            |
| Dmel\_rhi   |     8092 | A            | C            |
| Dmel\_rhi   |     8257 | G            | A            |
| Dmel\_rhi   |     8520 | G            | A            |
| Dmel\_rhi   |     8978 | T            | C            |
| Dmel\_rpl32 |      129 | A            | G            |
| Dmel\_rpl32 |      178 | T            | C            |
| Dmel\_rpl32 |      196 | A            | G            |
| Dmel\_rpl32 |      240 | T            | C            |
| Dmel\_rpl32 |      256 | T            | C            |
| Dmel\_rpl32 |      284 | C            | A            |
| Dmel\_rpl32 |      352 | G            | A            |
| Dmel\_rpl32 |      370 | G            | C            |
| Dmel\_rpl32 |      515 | C            | T            |
| Dmel\_rpl32 |      517 | A            | G            |
| Dmel\_rpl32 |      571 | C            | T            |
| Dmel\_rpl32 |      680 | A            | G            |
| Dmel\_rpl32 |      780 | A            | T            |
| Dmel\_rpl32 |      808 | G            | A            |
| Dmel\_rpl32 |      810 | C            | A            |
| Dmel\_rpl32 |      812 | A            | G            |
| Dmel\_rpl32 |      825 | T            | C            |
| Dmel\_rpl32 |      841 | T            | A            |
| Dmel\_rpl32 |      856 | G            | A            |
| Dmel\_rpl32 |      889 | C            | T            |
| Dmel\_rpl32 |     1444 | C            | T            |
| Dmel\_rpl32 |     1510 | C            | A            |
| Dmel\_rpl32 |     1513 | A            | G            |
| Dmel\_rpl32 |     1551 | A            | G            |
| Dmel\_rpl32 |     1586 | T            | A            |
| Dmel\_rpl32 |     1681 | A            | T            |
| Dmel\_rpl32 |     1692 | C            | T            |
| Dmel\_rpl32 |     1707 | A            | C            |
| Dmel\_rpl32 |     1709 | T            | C            |
| Dmel\_rpl32 |     2413 | T            | G            |
| Dmel\_rpl32 |     2449 | T            | A            |
| Dmel\_rpl32 |     2597 | C            | T            |
| Dmel\_rpl32 |     2640 | T            | G            |
| Dmel\_rpl32 |     2656 | T            | G            |
| Dmel\_rpl32 |     2693 | T            | A            |
| Dmel\_rpl32 |     2707 | T            | A            |
| Dmel\_rpl32 |     2767 | A            | G            |
| Dmel\_rpl32 |     2818 | T            | C            |
| Dmel\_rpl32 |     2834 | G            | A            |
| Dmel\_rpl32 |     2851 | C            | G            |
| Dmel\_rpl32 |     3474 | C            | T            |
| Dmel\_rpl32 |     3482 | A            | G            |
| Dmel\_rpl32 |     3507 | C            | T            |
| Dmel\_rpl32 |     3546 | T            | A            |
| Dmel\_rpl32 |     3548 | A            | C            |
| Dmel\_rpl32 |     3614 | T            | A            |
| Dmel\_rpl32 |     3616 | G            | T            |
| Dmel\_rpl32 |     3627 | G            | A            |
| Dmel\_rpl32 |     3784 | T            | G            |
| Dmel\_rpl32 |     3788 | A            | G            |
| Dmel\_rpl32 |     3882 | A            | G            |
| Dmel\_rpl32 |     3915 | T            | C            |
| Dmel\_rpl32 |     3928 | T            | C            |
| Dmel\_rpl32 |     3932 | G            | A            |
| Dmel\_rpl32 |     3939 | A            | C            |
| Dmel\_rpl32 |     3944 | G            | T            |
| Dmel\_rpl32 |     3949 | C            | T            |
| Dmel\_rpl32 |     4017 | T            | C            |
| Dmel\_rpl32 |     4209 | A            | C            |
| Dmel\_rpl32 |     4257 | C            | A            |
| Dmel\_rpl32 |     4297 | C            | A            |
| Dmel\_rpl32 |     4311 | G            | T            |
| Dmel\_rpl32 |     4383 | C            | T            |
| Dmel\_rpl32 |     4488 | T            | A            |
| Dmel\_rpl32 |     4635 | G            | A            |
| Dmel\_rpl32 |     4804 | T            | A            |
| Dmel\_tj    |      174 | T            | A            |
| Dmel\_tj    |      175 | A            | C            |
| Dmel\_tj    |      196 | A            | C            |
| Dmel\_tj    |      322 | T            | A            |
| Dmel\_tj    |      405 | A            | C            |
| Dmel\_tj    |      673 | C            | T            |
| Dmel\_tj    |      862 | T            | C            |
| Dmel\_tj    |     1575 | T            | C            |
| Dmel\_tj    |     1779 | A            | G            |
| Dmel\_tj    |     1792 | C            | G            |
| Dmel\_tj    |     1796 | C            | G            |
| Dmel\_tj    |     1961 | T            | A            |
| Dmel\_tj    |     1965 | C            | G            |
| Dmel\_tj    |     2049 | A            | C            |
| Dmel\_tj    |     2218 | G            | T            |
| Dmel\_tj    |     2222 | T            | G            |
| Dmel\_tj    |     2247 | G            | C            |
| Dmel\_tj    |     2272 | A            | G            |
| Dmel\_tj    |     2274 | C            | T            |
| Dmel\_tj    |     2390 | A            | C            |
| Dmel\_tj    |     2398 | A            | G            |
| Dmel\_tj    |     2537 | T            | C            |
| Dmel\_tj    |     2569 | G            | C            |
| Dmel\_tj    |     2579 | A            | G            |
| Dmel\_tj    |     2615 | T            | C            |
| Dmel\_tj    |     2648 | A            | G            |
| Dmel\_tj    |     2804 | T            | C            |
| Dmel\_tj    |     3032 | C            | A            |
| Dmel\_tj    |     3092 | A            | G            |
| Dmel\_tj    |     3416 | C            | T            |
| Dmel\_tj    |     3578 | G            | A            |
| Dmel\_tj    |     3581 | G            | C            |
| Dmel\_tj    |     3656 | T            | C            |
| Dmel\_tj    |     3737 | T            | A            |
| Dmel\_tj    |     3773 | G            | A            |
| Dmel\_tj    |     3860 | G            | T            |
| Dmel\_tj    |     3911 | C            | T            |
| Dmel\_tj    |     4096 | C            | A            |
| Dmel\_tj    |     4319 | C            | T            |
| Dmel\_tj    |     4348 | A            | C            |
| Dmel\_tj    |     4380 | G            | A            |
| Dmel\_tj    |     4418 | T            | C            |
| Dmel\_tj    |     4614 | A            | G            |
| Dmel\_tj    |     4911 | A            | C            |
| Dmel\_tj    |     5787 | G            | A            |
| Dmel\_tj    |     5965 | A            | C            |
| Dmel\_tj    |     5975 | A            | G            |
| Dmel\_tj    |     6187 | T            | A            |
| Dmel\_tj    |     6189 | T            | C            |
| Dmel\_tj    |     6280 | T            | G            |
| Dmel\_tj    |     6286 | T            | C            |
| Dmel\_tj    |     6407 | A            | C            |
| Dmel\_tj    |     6542 | T            | G            |
| Dmel\_tj    |     6578 | T            | G            |
| Dmel\_tj    |     6583 | T            | A            |
| Dmel\_tj    |     7003 | T            | C            |
| Dmel\_tj    |     7023 | A            | G            |
| Dmel\_tj    |     7120 | T            | A            |
| Dmel\_tj    |     7121 | T            | C            |
| Dmel\_tj    |     7124 | G            | C            |
| Dmel\_tj    |     7139 | C            | A            |

Dsim genes

``` r
s<-read.table("/Volumes/Temp2/filip/2040A/results/diagnostic_SNPs/define/withPop_Dsim-scg.SNPs")
names(s)<-c("gene","position","Dsim-variant","Dmel-variant")
knitr::kable(s)
```

| gene        | position | Dsim-variant | Dmel-variant |
| :---------- | -------: | :----------- | :----------- |
| Dsim\_rhi   |       51 | T            | G            |
| Dsim\_rhi   |       53 | C            | G            |
| Dsim\_rhi   |       58 | T            | G            |
| Dsim\_rhi   |       59 | T            | A            |
| Dsim\_rhi   |       78 | A            | C            |
| Dsim\_rhi   |       83 | A            | G            |
| Dsim\_rhi   |      116 | G            | T            |
| Dsim\_rhi   |      153 | C            | T            |
| Dsim\_rhi   |      154 | A            | T            |
| Dsim\_rhi   |      161 | G            | A            |
| Dsim\_rhi   |      164 | C            | G            |
| Dsim\_rhi   |      573 | A            | G            |
| Dsim\_rhi   |      845 | A            | G            |
| Dsim\_rhi   |      848 | T            | C            |
| Dsim\_rhi   |      854 | C            | T            |
| Dsim\_rhi   |      892 | A            | T            |
| Dsim\_rhi   |      897 | T            | C            |
| Dsim\_rhi   |      927 | G            | A            |
| Dsim\_rhi   |      973 | A            | G            |
| Dsim\_rhi   |      986 | G            | T            |
| Dsim\_rhi   |     1117 | A            | G            |
| Dsim\_rhi   |     1130 | A            | G            |
| Dsim\_rhi   |     1136 | T            | C            |
| Dsim\_rhi   |     1137 | A            | G            |
| Dsim\_rhi   |     1158 | A            | G            |
| Dsim\_rhi   |     1194 | G            | A            |
| Dsim\_rhi   |     1304 | T            | C            |
| Dsim\_rhi   |     1315 | A            | G            |
| Dsim\_rhi   |     1332 | C            | A            |
| Dsim\_rhi   |     1343 | C            | T            |
| Dsim\_rhi   |     1347 | G            | A            |
| Dsim\_rhi   |     1366 | T            | C            |
| Dsim\_rhi   |     1394 | T            | G            |
| Dsim\_rhi   |     1402 | C            | T            |
| Dsim\_rhi   |     1813 | C            | G            |
| Dsim\_rhi   |     1826 | C            | T            |
| Dsim\_rhi   |     1835 | T            | C            |
| Dsim\_rhi   |     1942 | T            | C            |
| Dsim\_rhi   |     2125 | T            | C            |
| Dsim\_rhi   |     2129 | A            | G            |
| Dsim\_rhi   |     2136 | A            | T            |
| Dsim\_rhi   |     2560 | C            | A            |
| Dsim\_rhi   |     2578 | G            | A            |
| Dsim\_rhi   |     2588 | G            | A            |
| Dsim\_rhi   |     2653 | G            | A            |
| Dsim\_rhi   |     2717 | T            | G            |
| Dsim\_rhi   |     2721 | A            | T            |
| Dsim\_rhi   |     2725 | A            | G            |
| Dsim\_rhi   |     2770 | C            | A            |
| Dsim\_rhi   |     2810 | C            | A            |
| Dsim\_rhi   |     2811 | G            | A            |
| Dsim\_rhi   |     2921 | C            | T            |
| Dsim\_rhi   |     3081 | T            | G            |
| Dsim\_rhi   |     3093 | A            | T            |
| Dsim\_rhi   |     3199 | G            | C            |
| Dsim\_rhi   |     3297 | C            | T            |
| Dsim\_rhi   |     3350 | C            | T            |
| Dsim\_rhi   |     3431 | G            | A            |
| Dsim\_rhi   |     3507 | G            | A            |
| Dsim\_rhi   |     3508 | A            | T            |
| Dsim\_rhi   |     3634 | C            | A            |
| Dsim\_rhi   |     3818 | G            | A            |
| Dsim\_rhi   |     3845 | C            | T            |
| Dsim\_rhi   |     4348 | C            | T            |
| Dsim\_rhi   |     4425 | C            | T            |
| Dsim\_rhi   |     4570 | A            | C            |
| Dsim\_rhi   |     4706 | C            | T            |
| Dsim\_rhi   |     4723 | C            | A            |
| Dsim\_rhi   |     4753 | A            | G            |
| Dsim\_rhi   |     4755 | G            | A            |
| Dsim\_rhi   |     4798 | T            | C            |
| Dsim\_rhi   |     4829 | A            | C            |
| Dsim\_rhi   |     4836 | T            | A            |
| Dsim\_rhi   |     4857 | C            | T            |
| Dsim\_rhi   |     4925 | G            | A            |
| Dsim\_rpl32 |      141 | C            | G            |
| Dsim\_rpl32 |      158 | T            | C            |
| Dsim\_rpl32 |      174 | G            | A            |
| Dsim\_rpl32 |      225 | C            | T            |
| Dsim\_rpl32 |      285 | T            | A            |
| Dsim\_rpl32 |      299 | T            | A            |
| Dsim\_rpl32 |      336 | C            | A            |
| Dsim\_rpl32 |      352 | C            | A            |
| Dsim\_rpl32 |      395 | A            | G            |
| Dsim\_rpl32 |      543 | T            | A            |
| Dsim\_rpl32 |      561 | A            | T            |
| Dsim\_rpl32 |      562 | G            | A            |
| Dsim\_rpl32 |      576 | C            | A            |
| Dsim\_tj    |       86 | G            | T            |
| Dsim\_tj    |      101 | C            | A            |
| Dsim\_tj    |      270 | T            | G            |
| Dsim\_tj    |      274 | G            | T            |
| Dsim\_tj    |      299 | C            | G            |
| Dsim\_tj    |      324 | G            | A            |
| Dsim\_tj    |      326 | T            | C            |
| Dsim\_tj    |      341 | A            | C            |
| Dsim\_tj    |      384 | A            | T            |
| Dsim\_tj    |      406 | T            | C            |
| Dsim\_tj    |      445 | C            | A            |
| Dsim\_tj    |      453 | G            | A            |
| Dsim\_tj    |      592 | C            | T            |
| Dsim\_tj    |      624 | C            | G            |
| Dsim\_tj    |      634 | G            | A            |
| Dsim\_tj    |      670 | C            | T            |
| Dsim\_tj    |      703 | G            | A            |
| Dsim\_tj    |      859 | C            | T            |
| Dsim\_tj    |     1087 | A            | C            |
| Dsim\_tj    |     1471 | T            | C            |
| Dsim\_tj    |     1697 | A            | G            |
| Dsim\_tj    |     1700 | C            | G            |
| Dsim\_tj    |     1775 | C            | T            |
| Dsim\_tj    |     1856 | A            | T            |
| Dsim\_tj    |     1892 | A            | G            |
| Dsim\_tj    |     1979 | T            | G            |
| Dsim\_tj    |     2030 | T            | C            |
| Dsim\_tj    |     2430 | T            | C            |
| Dsim\_tj    |     2459 | C            | A            |
| Dsim\_tj    |     2491 | A            | G            |
| Dsim\_tj    |     2529 | C            | T            |
| Dsim\_tj    |     2709 | A            | G            |
| Dsim\_tj    |     2754 | T            | A            |
| Dsim\_tj    |     2763 | C            | G            |
| Dsim\_tj    |     2769 | G            | A            |
| Dsim\_tj    |     2781 | C            | T            |
| Dsim\_tj    |     2783 | T            | C            |
| Dsim\_tj    |     3058 | C            | A            |
| Dsim\_tj    |     3192 | T            | C            |
| Dsim\_tj    |     3261 | C            | T            |
| Dsim\_tj    |     3389 | A            | C            |
| Dsim\_tj    |     3437 | G            | C            |
