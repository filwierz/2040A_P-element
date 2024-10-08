TE\_abundance
================
Filip Wierzbicki
4/11/2023

``` bash
#mapping
##dmel
nohup sh -c 'for i in Dmel_*;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dmel.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/TEs_3scg/${n}.sort.bam;done' &
##dsim
nohup sh -c 'for i in Dsim_*;do n=${i%.fq.gz};bwa bwasw -t 10 /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dsim.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp2/filip/2040A/map/TEs_3scg/${n}.sort.bam;done' &

#analysis:
cd /Volumes/Temp2/filip/2040A/map/TEs_3scg
conda activate deviaTE_env
for i in *bam;do ln -s /Volumes/Temp2/filip/2040A/map/TEs_3scg/${i} /Volumes/Temp2/filip/2040A/results/deviate/${i};done
cd /Volumes/Temp2/filip/2040A/results/deviate
for i in *bam;do samtools index $i;done
##dmel
nohup sh -c 'for i in Dmel*bam;do cat /Volumes/Temp2/filip/2040A/ref/TEnames.txt|while read TE;do deviaTE_analyse --input $i --single_copy_genes Dmel_tj,Dmel_rpl32,Dmel_rhi --library /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dmel.fasta --family $TE;done;done' &
##dsim
nohup sh -c 'for i in Dsim*bam;do cat /Volumes/Temp2/filip/2040A/ref/TEnames.txt|while read TE;do deviaTE_analyse --input $i --single_copy_genes Dsim_tj,Dsim_rpl32,Dsim_rhi --library /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dsim.fasta --family $TE;done;done' &
###note that those loops are very time consuming. If this becomes a bottleneck, I can write a storm script that submits multiple commands to parallelize the process
##summary for the P-element:
for i in *PPI251;do awk '$2=="insertions/haploid:"' $i|awk -v a="$i" '{print $3,a}';done|sed 's/.sort.bam.PPI251//g'|awk -F "_" '{print $1,$2,$3,$4,$5}' > forR/PPI241-copynr.forR
#batch2:
cd /Volumes/Temp2/filip/2040A/data/batch2/fastq
nohup sh -c 'for i in *.fq.gz;do n=${i%.fq.gz};bwa bwasw -t 20 /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dmel.fasta $i|samtools sort -@ 4 -m 3G - > /Volumes/Temp3/filip/2040A/map/TEs_3scg/batch2/${n}.sort.bam;done' &
cat /Volumes/Temp2/filip/2040A/data/batch2/metadata_batch2.txt|while read b n;do ln -s /Volumes/Temp3/filip/2040A/map/TEs_3scg/batch2/${b}.fastq.sort.bam /Volumes/Temp2/filip/2040A/results/deviate/batch2/${n}.sort.bam;done
/Volumes/Temp2/filip/2040A/results/deviate/batch2
nohup sh -c 'for i in *bam;do deviaTE_analyse --input $i --single_copy_genes Dmel_tj,Dmel_rpl32,Dmel_rhi --library /Volumes/Temp2/filip/2040A/ref/TEs-scg/teseqs-3scg-dmel.fasta --family PPI251;done'&
```

``` r
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
t<-read.table("/Volumes/Temp2/filip/2040A/results/deviate/forR/PPI241-copynr.forR")
t2<-read.table("/Volumes/Temp2/filip/2040A/results/deviate/batch2/forR/PPI251-copynr.forR")
t<-rbind(t,t2)
names(t)<-c("copies","species","variant","temperature","generation","replicate")
t$type<-paste(t$temperature,t$variant,t$species,sep="_")

base<-subset(t,type=="tX_S_Dsim")
baco<-base
baco$temperature<-c("t20-10")
baco$type<-paste(baco$temperature,baco$variant,baco$species,sep="_")

baho<-base
baho$temperature<-c("t28-18")
baho$type<-paste(baho$temperature,baho$variant,baho$species,sep="_")

#t<-subset(t,type!="t25_M_Dmel")
#t<-subset(t,type!="t25_S_Dmel")
t<-subset(t,type!="tX_S_Dsim")

t<-rbind(baco,baho,t)

t$temperature<-gsub("t","",t$temperature)
t$generation<-as.numeric(gsub("g","",t$generation))
t$type1<-paste(t$replicate,t$type,sep="_")
t$type2<-paste(t$species,t$variant,sep="_")

labs <- c("2040G", "2040A")
names(labs) <- c("M","S")

##consistent replicate labels
repl<-read.table("/Volumes/Temp2/filip/2040A/ref/replicate_labels.txt")
names(repl)<-c("id","rep")
t$id<-paste(t$species,t$variant,t$replicate,sep="_")
t<-left_join(t,repl,by="id")
cols <- c("20-10"="#3393FF", "28-18"="#FF3333", "25"="#10CB40")
g<-ggplot(t, aes(x=generation, y=copies,color=temperature,by=type1)) + geom_line()+facet_grid(species ~ variant,scales="free",labeller = labeller(variant=labs))+ylab("insertions per haploid") +scale_color_manual(values = cols)
plot(g)
```

![](TE_abundance_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
ggsave("/Volumes/Temp2/filip/2040A/results/figures/invasion_dynamics.pdf",width=8,height=6)
ggsave("/Volumes/Temp2/filip/2040A/results/figures/invasion_dynamics.png",width=8,height=6)
```
