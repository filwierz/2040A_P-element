data\_overview
================
Filip Wierzbicki
3/21/2023

Documentation of the processing steps in the command line

``` bash
#strain data:
awk '{print $1}' SRR8247551.sra.fastq|paste - - - -|awk 'length($2)==140'|tr "\t" "\n"|gzip -c > wXD1.fq.gz #this dataset has mixed length, keep only 140nt reads
for i in *1_mod.fastq.gz;do n=${i%_1_mod.fastq.gz};gzip -cd ${i} ${n}_2_mod.fastq.gz|gzip -c > ../${n}.fq.gz;done #merging PE data into 1 file
for i in *gz;do gzip -cd $i|awk '{print $1}'|cut -c-100|gzip -c > trimmed/${i};done #trims the fastq file to 100bp

#divyas data:
cd /Volumes/Temp2/filip/2040A/data/divya
for i in *_1.fq.gz;do n=${i%_1.fq.gz}; gzip -cd ${i} ${n}_2.fq.gz|gzip -c > merged_PE/${n}.fq.gz ;done #merging PE data into 1 file
for i in *gz;do n=${i#*sim_};n=${n/G/g};n=${n/R/r};n=${n%.fq.gz};n=${n%_*};n=${n/naive/r1g0};r=${n%g*};g=${n#*g};mv $i Dsim_M_t25_g${g}_${r}.fq.gz;done #renaming files
for i in *gz;do gzip -cd $i|awk '{print $1}'|cut -c-100|gzip -c > trimmed/${i};done #trims the fastq file to 100bp

#roberts data:
cd /Volumes/Temp2/filip/2040A/data/hotcold
for i in *gz;do n=${i/hot/t28-18_g};n=${n/cold/t20-10_g};n=${n/base/tX_g0};mv $i Dsim_S_${n};done #renaming files

#pool all files into 1 directory:
cd /Volumes/Temp2/filip/2040A/data/divya/merged_PE/trimmed
for i in *gz;do ln -s /Volumes/Temp2/filip/2040A/data/divya/merged_PE/trimmed/${i} /Volumes/Temp2/filip/2040A/data/pool/${i};done
cd /Volumes/Temp2/filip/2040A/data/hotcold/
for i in *gz;do ln -s /Volumes/Temp2/filip/2040A/data/hotcold/${i} /Volumes/Temp2/filip/2040A/data/pool/${i};done
cd /Volumes/Temp2/filip/2040A/data/filip/
for i in *gz;do ln -s /Volumes/Temp2/filip/2040A/data/filip/${i} /Volumes/Temp2/filip/2040A/data/pool/${i};done


#overview table:
cd /Volumes/Temp2/filip/2040A/data/pool
#basis:
for i in *gz;do echo ${i%.fq.gz}|awk -F "_" -v a="$i" '{print $1,$2,$3,$4,$5,a}';done|sort -k6,6 > /Volumes/Temp2/filip/2040A/results/overview/samples/sampleIDs.txt 
##trimmed read length:
for i in *gz;do gzip -cd $i|head -2|tail -1|awk -v a="$i" '{print length($1),a}';done|sort -k2,2 > /Volumes/Temp2/filip/2040A/results/overview/samples/trimmedRLs.txt
##read length:
cd /Volumes/Temp2/filip/2040A/data/filip
for i in *gz;do gzip -cd $i|head -2|tail -1|awk -v a="$i" '{print length($1),a}';done > /Volumes/Temp2/filip/2040A/results/overview/samples/origRLs-filip.txt
cd /Volumes/Temp2/filip/2040A/data/divya/merged_PE
for i in *gz;do gzip -cd $i|head -2|tail -1|awk -v a="$i" '{print length($1),a}';done > /Volumes/Temp2/filip/2040A/results/overview/samples/origRLs-divya.txt
cd /Volumes/Temp2/filip/2040A/data/hotcold
for i in *gz;do gzip -cd $i|head -2|tail -1|awk -v a="$i" '{print length($1),a}';done > /Volumes/Temp2/filip/2040A/results/overview/samples/origRLs-robert.txt #required manual editing based on table S1 from  https://doi.org/10.1093/molbev/msac141 since data is trimmed already

##merging files:
cd /Volumes/Temp2/filip/2040A/results/overview/samples
cat origRLs-*|sort -k2,2 > originalRLs.txt
join -1 2 -2 2 originalRLs.txt trimmedRLs.txt|awk '{print $2"("$3")",$1}' > RLs.txt
join -1 6 -2 2 sampleIDs.txt RLs.txt|awk '{print $2,$3,$4,$5,$6,$7,$1}' > overview.forR
for i in *.gz;do gzip -cd $i|paste - - - -|wc -l|awk -v a="$i" '{print $1,a}';done > /Volumes/Temp2/filip/2040A/results/overview/samples/readcounts.txt

#quality check
fastqc -o /Volumes/Temp2/filip/2040A/results/overview/fastqc -d /Volumes/Temp2/filip/2040A/results/overview/fastqc/temp -t 10 Dmel_M_t25_g34_r1.fq.gz Dmel_M_t25_g34_r2.fq.gz Dmel_S_t25_g34_r1.fq.gz Dmel_S_t25_g34_r2.fq.gz Dsim_M_t25_g0_r1.fq.gz Dsim_M_t25_g10_r2.fq.gz Dsim_M_t25_g10_r3.fq.gz Dsim_M_t25_g10_r4.fq.gz Dsim_M_t25_g1_r2.fq.gz Dsim_M_t25_g1_r3.fq.gz Dsim_M_t25_g1_r4.fq.gz Dsim_M_t25_g20_r2.fq.gz Dsim_M_t25_g20_r3.fq.gz Dsim_M_t25_g20_r4.fq.gz Dsim_M_t25_g34_r2.fq.gz Dsim_M_t25_g34_r3.fq.gz Dsim_M_t25_g34_r4.fq.gz Dsim_M_t25_g40_r2.fq.gz Dsim_M_t25_g40_r3.fq.gz Dsim_M_t25_g40_r4.fq.gz Dsim_M_t25_g48_r2.fq.gz Dsim_M_t25_g48_r3.fq.gz Dsim_M_t25_g48_r4.fq.gz Dsim_M_t25_g63_r2.fq.gz Dsim_M_t25_g63_r3.fq.gz Dsim_M_t25_g63_r4.fq.gz Dsim_S_t20-10_g100_r1.fq.gz Dsim_S_t20-10_g100_r3.fq.gz Dsim_S_t20-10_g100_r5.fq.gz Dsim_S_t20-10_g10_r1.fq.gz Dsim_S_t20-10_g10_r3.fq.gz Dsim_S_t20-10_g10_r5.fq.gz Dsim_S_t20-10_g20_r1.fq.gz Dsim_S_t20-10_g20_r3.fq.gz Dsim_S_t20-10_g20_r5.fq.gz Dsim_S_t20-10_g30_r1.fq.gz Dsim_S_t20-10_g30_r3.fq.gz Dsim_S_t20-10_g30_r5.fq.gz Dsim_S_t20-10_g40_r1.fq.gz Dsim_S_t20-10_g40_r3.fq.gz Dsim_S_t20-10_g40_r5.fq.gz Dsim_S_t20-10_g50_r1.fq.gz Dsim_S_t20-10_g50_r3.fq.gz Dsim_S_t20-10_g50_r5.fq.gz Dsim_S_t20-10_g60_r1.fq.gz Dsim_S_t20-10_g60_r3.fq.gz Dsim_S_t20-10_g60_r5.fq.gz Dsim_S_t20-10_g70_r1.fq.gz Dsim_S_t20-10_g70_r3.fq.gz Dsim_S_t20-10_g70_r5.fq.gz Dsim_S_t20-10_g80_r1.fq.gz Dsim_S_t20-10_g80_r3.fq.gz Dsim_S_t20-10_g80_r5.fq.gz Dsim_S_t20-10_g90_r1.fq.gz Dsim_S_t20-10_g90_r3.fq.gz Dsim_S_t20-10_g90_r5.fq.gz Dsim_S_t28-18_g10_r1.fq.gz Dsim_S_t28-18_g10_r3.fq.gz Dsim_S_t28-18_g10_r5.fq.gz Dsim_S_t28-18_g20_r1.fq.gz Dsim_S_t28-18_g20_r3.fq.gz Dsim_S_t28-18_g20_r5.fq.gz Dsim_S_t28-18_g30_r1.fq.gz Dsim_S_t28-18_g30_r3.fq.gz Dsim_S_t28-18_g30_r5.fq.gz Dsim_S_t28-18_g40_r1.fq.gz Dsim_S_t28-18_g40_r3.fq.gz Dsim_S_t28-18_g40_r5.fq.gz Dsim_S_t28-18_g50_r1.fq.gz Dsim_S_t28-18_g50_r3.fq.gz Dsim_S_t28-18_g50_r5.fq.gz Dsim_S_t28-18_g60_r1.fq.gz Dsim_S_t28-18_g60_r3.fq.gz Dsim_S_t28-18_g60_r5.fq.gz Dsim_S_tX_g0_r1.fq.gz Dsim_S_tX_g0_r3.fq.gz Dsim_S_tX_g0_r5.fq.gz
```

Generating the final table

``` r
t<-read.table("/Volumes/Temp2/filip/2040A/results/overview/samples/overview.forR")
t$V3<-gsub("t","",as.character(t$V3))
t$V4<-gsub("g","",as.character(t$V4))
t$V5<-gsub("r","",as.character(t$V5))
names(t)<-c("species","variant","temp.","gen.","rep.","rl.(trimmed)","file")

r<-read.table("/Volumes/Temp2/filip/2040A/results/overview/samples/readcounts.txt")
names(r)<-c("rc","file")
r$rc<-paste (round(r$rc/1000000,1),"M",sep = "")
names(r)<-c("rc.","file")

tr<-merge(r,t,by="file")
tr<-subset(tr,select=c("species","variant","temp.","gen.","rep.","rl.(trimmed)","rc.","file"))
knitr::kable(tr)
```

| species | variant | temp. | gen. | rep. | rl.(trimmed) | rc.    | file                            |
| :------ | :------ | :---: | :--- | :--- | :----------- | :----- | :------------------------------ |
| Dmel    | M       |  25   | 34   | 1    | 100(100)     | 60.8M  | Dmel\_M\_t25\_g34\_r1.fq.gz     |
| Dmel    | M       |  25   | 34   | 2    | 100(100)     | 72.9M  | Dmel\_M\_t25\_g34\_r2.fq.gz     |
| Dmel    | S       |  25   | 34   | 1    | 100(100)     | 78.7M  | Dmel\_S\_t25\_g34\_r1.fq.gz     |
| Dmel    | S       |  25   | 34   | 2    | 100(100)     | 74.7M  | Dmel\_S\_t25\_g34\_r2.fq.gz     |
| Dsim    | M       |  25   | 0    | 1    | 125(100)     | 34.6M  | Dsim\_M\_t25\_g0\_r1.fq.gz      |
| Dsim    | M       |  25   | 1    | 2    | 125(100)     | 54.5M  | Dsim\_M\_t25\_g1\_r2.fq.gz      |
| Dsim    | M       |  25   | 1    | 3    | 125(100)     | 53.7M  | Dsim\_M\_t25\_g1\_r3.fq.gz      |
| Dsim    | M       |  25   | 1    | 4    | 125(100)     | 63.3M  | Dsim\_M\_t25\_g1\_r4.fq.gz      |
| Dsim    | M       |  25   | 10   | 2    | 125(100)     | 62.4M  | Dsim\_M\_t25\_g10\_r2.fq.gz     |
| Dsim    | M       |  25   | 10   | 3    | 125(100)     | 46.9M  | Dsim\_M\_t25\_g10\_r3.fq.gz     |
| Dsim    | M       |  25   | 10   | 4    | 125(100)     | 44.7M  | Dsim\_M\_t25\_g10\_r4.fq.gz     |
| Dsim    | M       |  25   | 20   | 2    | 125(100)     | 46.1M  | Dsim\_M\_t25\_g20\_r2.fq.gz     |
| Dsim    | M       |  25   | 20   | 3    | 125(100)     | 57.9M  | Dsim\_M\_t25\_g20\_r3.fq.gz     |
| Dsim    | M       |  25   | 20   | 4    | 125(100)     | 54.9M  | Dsim\_M\_t25\_g20\_r4.fq.gz     |
| Dsim    | M       |  25   | 34   | 2    | 125(100)     | 68.4M  | Dsim\_M\_t25\_g34\_r2.fq.gz     |
| Dsim    | M       |  25   | 34   | 3    | 125(100)     | 47.1M  | Dsim\_M\_t25\_g34\_r3.fq.gz     |
| Dsim    | M       |  25   | 34   | 4    | 125(100)     | 50M    | Dsim\_M\_t25\_g34\_r4.fq.gz     |
| Dsim    | M       |  25   | 40   | 2    | 125(100)     | 45.1M  | Dsim\_M\_t25\_g40\_r2.fq.gz     |
| Dsim    | M       |  25   | 40   | 3    | 125(100)     | 53.7M  | Dsim\_M\_t25\_g40\_r3.fq.gz     |
| Dsim    | M       |  25   | 40   | 4    | 125(100)     | 51.4M  | Dsim\_M\_t25\_g40\_r4.fq.gz     |
| Dsim    | M       |  25   | 48   | 2    | 125(100)     | 47.3M  | Dsim\_M\_t25\_g48\_r2.fq.gz     |
| Dsim    | M       |  25   | 48   | 3    | 125(100)     | 54.7M  | Dsim\_M\_t25\_g48\_r3.fq.gz     |
| Dsim    | M       |  25   | 48   | 4    | 125(100)     | 50M    | Dsim\_M\_t25\_g48\_r4.fq.gz     |
| Dsim    | M       |  25   | 63   | 2    | 100(100)     | 68.5M  | Dsim\_M\_t25\_g63\_r2.fq.gz     |
| Dsim    | M       |  25   | 63   | 3    | 100(100)     | 73.8M  | Dsim\_M\_t25\_g63\_r3.fq.gz     |
| Dsim    | M       |  25   | 63   | 4    | 100(100)     | 90.3M  | Dsim\_M\_t25\_g63\_r4.fq.gz     |
| Dsim    | S       | 20-10 | 10   | 1    | 100(100)     | 69.4M  | Dsim\_S\_t20-10\_g10\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 10   | 3    | 100(100)     | 81.2M  | Dsim\_S\_t20-10\_g10\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 10   | 5    | 100(100)     | 68.3M  | Dsim\_S\_t20-10\_g10\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 100  | 1    | 150(100)     | 82M    | Dsim\_S\_t20-10\_g100\_r1.fq.gz |
| Dsim    | S       | 20-10 | 100  | 3    | 150(100)     | 79.5M  | Dsim\_S\_t20-10\_g100\_r3.fq.gz |
| Dsim    | S       | 20-10 | 100  | 5    | 150(100)     | 75.8M  | Dsim\_S\_t20-10\_g100\_r5.fq.gz |
| Dsim    | S       | 20-10 | 20   | 1    | 100(100)     | 85.3M  | Dsim\_S\_t20-10\_g20\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 20   | 3    | 100(100)     | 50M    | Dsim\_S\_t20-10\_g20\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 20   | 5    | 100(100)     | 48.8M  | Dsim\_S\_t20-10\_g20\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 30   | 1    | 100(100)     | 58.4M  | Dsim\_S\_t20-10\_g30\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 30   | 3    | 100(100)     | 65.4M  | Dsim\_S\_t20-10\_g30\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 30   | 5    | 100(100)     | 66.2M  | Dsim\_S\_t20-10\_g30\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 40   | 1    | 120(100)     | 38.7M  | Dsim\_S\_t20-10\_g40\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 40   | 3    | 120(100)     | 93.8M  | Dsim\_S\_t20-10\_g40\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 40   | 5    | 120(100)     | 43.5M  | Dsim\_S\_t20-10\_g40\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 50   | 1    | 150(100)     | 77.4M  | Dsim\_S\_t20-10\_g50\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 50   | 3    | 150(100)     | 82.1M  | Dsim\_S\_t20-10\_g50\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 50   | 5    | 150(100)     | 77.7M  | Dsim\_S\_t20-10\_g50\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 60   | 1    | 150(100)     | 26.6M  | Dsim\_S\_t20-10\_g60\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 60   | 3    | 150(100)     | 24.5M  | Dsim\_S\_t20-10\_g60\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 60   | 5    | 150(100)     | 66M    | Dsim\_S\_t20-10\_g60\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 70   | 1    | 150(100)     | 79.7M  | Dsim\_S\_t20-10\_g70\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 70   | 3    | 150(100)     | 90.7M  | Dsim\_S\_t20-10\_g70\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 70   | 5    | 150(100)     | 101.2M | Dsim\_S\_t20-10\_g70\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 80   | 1    | 150(100)     | 69.3M  | Dsim\_S\_t20-10\_g80\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 80   | 3    | 150(100)     | 74.8M  | Dsim\_S\_t20-10\_g80\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 80   | 5    | 150(100)     | 84.9M  | Dsim\_S\_t20-10\_g80\_r5.fq.gz  |
| Dsim    | S       | 20-10 | 90   | 1    | 150(100)     | 79.9M  | Dsim\_S\_t20-10\_g90\_r1.fq.gz  |
| Dsim    | S       | 20-10 | 90   | 3    | 150(100)     | 84M    | Dsim\_S\_t20-10\_g90\_r3.fq.gz  |
| Dsim    | S       | 20-10 | 90   | 5    | 150(100)     | 80.1M  | Dsim\_S\_t20-10\_g90\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 10   | 1    | 100(100)     | 68.2M  | Dsim\_S\_t28-18\_g10\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 10   | 3    | 100(100)     | 54M    | Dsim\_S\_t28-18\_g10\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 10   | 5    | 100(100)     | 89.4M  | Dsim\_S\_t28-18\_g10\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 20   | 1    | 100(100)     | 175.9M | Dsim\_S\_t28-18\_g20\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 20   | 3    | 100(100)     | 170M   | Dsim\_S\_t28-18\_g20\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 20   | 5    | 100(100)     | 158.1M | Dsim\_S\_t28-18\_g20\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 30   | 1    | 100(100)     | 81.5M  | Dsim\_S\_t28-18\_g30\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 30   | 3    | 100(100)     | 73.1M  | Dsim\_S\_t28-18\_g30\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 30   | 5    | 100(100)     | 44.4M  | Dsim\_S\_t28-18\_g30\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 40   | 1    | 100(100)     | 117.9M | Dsim\_S\_t28-18\_g40\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 40   | 3    | 100(100)     | 151.7M | Dsim\_S\_t28-18\_g40\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 40   | 5    | 100(100)     | 100.8M | Dsim\_S\_t28-18\_g40\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 50   | 1    | 100(100)     | 58.7M  | Dsim\_S\_t28-18\_g50\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 50   | 3    | 100(100)     | 60M    | Dsim\_S\_t28-18\_g50\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 50   | 5    | 100(100)     | 72.9M  | Dsim\_S\_t28-18\_g50\_r5.fq.gz  |
| Dsim    | S       | 28-18 | 60   | 1    | 100(100)     | 97.3M  | Dsim\_S\_t28-18\_g60\_r1.fq.gz  |
| Dsim    | S       | 28-18 | 60   | 3    | 100(100)     | 55.9M  | Dsim\_S\_t28-18\_g60\_r3.fq.gz  |
| Dsim    | S       | 28-18 | 60   | 5    | 100(100)     | 63.1M  | Dsim\_S\_t28-18\_g60\_r5.fq.gz  |
| Dsim    | S       |   X   | 0    | 1    | 100(100)     | 147.1M | Dsim\_S\_tX\_g0\_r1.fq.gz       |
| Dsim    | S       |   X   | 0    | 3    | 100(100)     | 211.7M | Dsim\_S\_tX\_g0\_r3.fq.gz       |
| Dsim    | S       |   X   | 0    | 5    | 100(100)     | 244.2M | Dsim\_S\_tX\_g0\_r5.fq.gz       |

``` r
write.table(tr,"/Volumes/Temp2/filip/2040A/results/overview/samples/final-overview.txt",quote=FALSE,sep="\t",row.names = FALSE)
```
