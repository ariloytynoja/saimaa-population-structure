#bcftools view -R first_ten.txt saimaa_posmask_98.vcf.gz | vcftools --vcf - --012 --out saimaa_test
#gzip saimaa_test.012
#gzip saimaa_test.012.pos

library(IRanges)
library(tidyr)
library(ggplot2)
library(data.table)
library(dplyr)
library(scales)  

f <- read.table("data/norppa_test.fa.fai")

m <- as.data.frame(fread('gunzip -cq data/norppa_test.posmask.bed.gz',sep="\t"))
colnames(m) <- c("ctg","start","end")
m$start <- m$start+1
m$end <- m$end+1

d <- as.data.frame(t(fread('gunzip -cq data/saimaa_test.012.gz'))) # 012-genotypes generated with vcftools, transposed
d <- d[-1,]
colnames(d) <-read.table("data/saimaa_test.012.indv",header=F)$V1

pos <- as.data.frame(fread('gunzip -cq data/saimaa_test.012.pos.gz',sep="\t"))
colnames(pos) <- c("ctg","pos")

d <- cbind.data.frame(pos,d)

ns  <- read.table("data/ids_ns.txt")$V1      # lists of sample ids per region
cs  <- read.table("data/ids_cs.txt")$V1
ss  <- read.table("data/ids_ss.txt")$V1
all <- c(ns,cs,ss)

###

ww <- 250000
dc <- c()

for(ct in f[,1]){

    d1 <- d[d$ctg==ct,]            # select one contig at time

    d2 <- cbind(d1,floor(d1$pos/ww)+1)            # add a new column indicating the window (each 250kbp in size)
    colnames(d2)[dim(d2)[2]] <- "window"          # label it "window"

    if(dim(d2)[1]>2){
        
        d3 <- aggregate(d2[,3:dim(d)[2]],         # count heterozygous sites per sample and per window
                    by=list(d2$window),FUN=function(x){sum(x==1)})
        colnames(d3)[1] <- "window"

        m1 <- m[m[,1]==ct,]                       # select pos. masked regions for this contig
        m2 <- IRanges(start=m1$start,end=m1$end)  # convert into IRanges

        w1 <- seq(1,f$V2[f$V1==ct],ww)            # make 250kbp windows for this contig
        w2 <- IRanges(start=w1,width=ww)          # convert into IRanges

        s <- rep(0,dim(d3)[1])
        t <- rep(0,dim(d3)[1])
        
        for(i in end(w2)/ww){                     # iterate over each window ("i" being the index)
            s[i] <- sum(width( restrict(m2,start(w2)[i],end(w2)[i]) ))   # sum of pos. masked regions
            t[i] <- i                                                    # window index
        }
        
        s <- cbind(t,s)
        colnames(s) <- c("window","nsites")

        d4 <- merge(d3,s,by.x = 1, by.y = 1)      # merge the het. counts and pos. mask counts

        dc <- rbind.data.frame(dc,cbind(ct,d4))   # add them to the final table, written to the disk in the end
    }
}

write.table(dc,"data/saimaa_test_H.tsv",row.names=F,quote=F)

###

pi_for_set <- function(set,d){
    n.alleles <- length(set)*2
    count.alt <- rowSums(d[,set]);
    count.ref <- n.alleles-count.alt
    count.ref * count.alt * 2 / ( n.alleles * (n.alleles-1) )
}

pi_for_two_sets <- function(set1,set2,d){
    n.alleles1 <- length(set1)*2
    n.alleles2 <- length(set2)*2
    count.alt1 <- rowSums(d[,set1]);
    count.alt2 <- rowSums(d[,set2]);
    count.ref1 <- n.alleles1-count.alt1
    count.ref2 <- n.alleles2-count.alt2
    ( count.ref1 * count.alt2 + count.alt1 * count.ref2 ) / ( n.alleles1 * n.alleles2 )
}

###

ww <- 250000
dc <- c()

ns.c  <- match(ns,gsub("X","",colnames(d)))  # columns for each subpoops
cs.c  <- match(cs,gsub("X","",colnames(d)))
ss.c  <- match(ss,gsub("X","",colnames(d)))
all.c <- match(all,gsub("X","",colnames(d)))

for(ct in f[,1]){

    d1 <- d[d$ctg==ct,]                           # select one contig at time

    d2 <- cbind(d1,floor(d1$pos/ww)+1)            # add a new column indicating the window (each 250kbp in size)
    colnames(d2)[dim(d2)[2]] <- "window"          # label it "window"
   
    if(dim(d2)[1]>2){

        pis <- cbind.data.frame(ns=pi_for_set(ns.c,d2),  # compute pi for NS
                                cs=pi_for_set(cs.c,d2),
                                ss=pi_for_set(ss.c,d2),
                                all=pi_for_set(all.c,d2),
                                ns.cs=pi_for_two_sets(ns.c,cs.c,d2), # compute dxy between NS and CS 
                                ns.ss=pi_for_two_sets(ns.c,ss.c,d2), 
                                cs.ss=pi_for_two_sets(cs.c,ss.c,d2)
                               ) 
        
        d3 <- aggregate(pis,by=list(d2$window),FUN=sum)  # sum per window
        colnames(d3)[1] <- "window"
        
        m1 <- m[m[,1]==ct,]                       # select pos. masked regions for this contig
        m2 <- IRanges(start=m1$start,end=m1$end)  # convert into IRanges

        w1 <- seq(1,f$V2[f$V1==ct],ww)            # make 250kbp windows for this contig
        w2 <- IRanges(start=w1,width=ww)          # convert into IRanges

        s <- rep(0,dim(d3)[1])
        t <- rep(0,dim(d3)[1])
        
        for(i in end(w2)/ww){                     # iterate over each window ("i" being the index)
            s[i] <- sum(width( restrict(m2,start(w2)[i],end(w2)[i]) ))   # sum of pos. masked regions
            t[i] <- i                                                    # window index
        }

        s <- cbind(t,s)
        colnames(s) <- c("window","nsites")

        d4 <- merge(d3,s,by.x = 1, by.y = 1)      # merge the pi/dxy counts and pos. mask counts

        dc <- rbind.data.frame(dc,cbind(ct,d4))  # add them to the final table, written to the disk in the end
    }
}

write.table(dc,"data/saimaa_test_pi.tsv",row.names=F,quote=F)

###

het <- read.table("data/saimaa_test_H.tsv",head=T)

het[,3:100] <- het[,3:100]/het[,101]

het.bins <- apply(het[het$nsites>100000,3:100], 2, 
                    function(x) hist(x,breaks=c(-0.001,seq(0.0002,0.0052,0.0002),0.027),plot=F)$counts )

het.bins <- cbind.data.frame(bins=c(seq(0,0.00519,0.0002),0.025),het.bins)

het.bins.m <- gather(het.bins, variable, value, -bins)
colnames(het.bins.m) <- c("bin","smp","count")
                   
het.bins.m$smp2="other"
het.bins.m$smp2[het.bins.m$smp=="X2575"]="2575"

ggplot()+
 geom_line(data=subset(het.bins.m,smp2=="other"),aes(x=bin,y=count,group=smp,col=smp2),size=0.25)+
 geom_point(data=subset(het.bins.m,smp2=="other"),aes(x=bin,y=count,group=smp,col=smp2),shape=20,size=0.25)+
 geom_boxplot(data=subset(het.bins.m,smp2!="2575"),aes(x=bin,y=count,group=bin),color="darkgray",outlier.size = 0.25)+
 geom_line(data=subset(het.bins.m,smp2=="2575"),aes(x=bin,y=count,group=smp,col=smp2),size=0.25,linetype="solid")+
 geom_point(data=subset(het.bins.m,smp2=="2575"),aes(x=bin,y=count,group=smp,col=smp2),shape=20,size=4)+
 theme_classic()+theme(legend.position="top")+
 xlim(-0.0001,0.0051)+
 ylab("#Windows")+xlab("Nucleotide diversity")+labs(color="")+xlab("Nucleotide diversity")

###

pi <- read.table("data/saimaa_test_pi.tsv",head=T)
pi[,3:9] <- pi[,3:9]/pi[,10]

pi.bins <- apply(pi[pi$nsites>100000,3:9], 2,
                   function(x) hist(x,breaks=c(-0.001,seq(0.0001,0.0052,0.0001),0.027),plot=F)$counts )
               
pi.bins <- cbind.data.frame(bins=c(seq(0,0.00519,0.0001),0.025),pi.bins)

pi.bins.m <- gather(pi.bins, variable, value, -bins)
colnames(pi.bins.m) <- c("bin","pop","count")

pi.bins.means = pi.bins.m %>% filter(pop %in% c("ns","cs","ss")) %>% group_by(bin) %>% summarise(mean(count))
colnames(pi.bins.means)=c("bin","mean.count")
pi.bins.means$pop="mean"

ggplot()+
 geom_line(data=subset(pi.bins.m,pop %in% c("ns","cs","ss")),aes(x=bin,y=count,color=pop),size=0.25,linetype="solid")+
 geom_point(data=subset(pi.bins.m,pop %in% c("ns","cs","ss")),aes(x=bin,y=count,color=pop),shape=19,size=1)+
 geom_line(data=pi.bins.means,aes(x=bin,y=mean.count,color=pop),size=1,linetype="solid")+
 geom_point(data=pi.bins.means,aes(x=bin,y=mean.count,color=pop),shape=19,size=2)+
 geom_point(data=subset(pi.bins.m,pop=="all"),aes(x=bin,y=count,color=pop),shape=19,size=2)+
 geom_line(data=subset(pi.bins.m,pop=="all"),aes(x=bin,y=count,color=pop),size=1)+
 theme_classic()+theme(legend.position="top")+
 xlim(-0.0001,0.0051)+
 ylab("#Windows")+xlab("")+labs(color="")+
 scale_color_manual(breaks=c("ns", "cs", "ss", "mean", "all"),
                       labels=c("NS", "CS", "SS", "Average w. regions", "Total Saimaa"),
                       values=hue_pal()(5))
