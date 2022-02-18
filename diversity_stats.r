library(IRanges)
library(reshape2)
library(ggplot2)
library(arrangements)

library(scales)  
col6=hue_pal()(6)
col3=hue_pal()(3)

system(paste(
    "bcftools query -f '%CHROM\t%POS[\t%GT]\n' -H data/saimaa_all_posmask_98.vcf.gz |",
    " sed 's/0|0/0/g;s/1|1/2/g;s/0|1/1/g;s/1|0/1/g;s/# //;s/:GT//g;s/\\[[0-9]*\\]//g' |",
    " bgzip -c > data/saimaa_all_posmask_98.tsv.gz")
)

# compute the number of base differences for four chromosomes (given as number of REF alleles: 0, 1 or 2)
#
pi_for_four <- function(set,d){
    rowMeans(apply(set,1,function(x){
        (d[,x[1]]+d[,x[2]]==1)*1/2 +
        (d[,x[1]]+d[,x[2]]==2)*2/3 +
        (d[,x[1]]+d[,x[2]]==3)*1/2
    }))
}

ns=read.table("data/ns_keep.txt")$V1
cs=read.table("data/cs_keep.txt")$V1
ss=read.table("data/ss_keep.txt")$V1
as=c(ns,cs,ss)

# Saimaa IDs are numbers and "X" is added when reading in
#
nsx=paste0("X",ns)
csx=paste0("X",cs)
ssx=paste0("X",ss)
asx=paste0("X",as)

f=read.table("data/norppa_12122017.fa.fai")

zz=gzfile("data/norppa_12122017.posmask.bed.gz",'rt') 
m=read.table(zz)
colnames(m)=c("ctg","start","end")
m$start=m$start+1
m$end=m$end+1

zz=gzfile("data/saimaa_all_posmask_98.tsv.gz",'rt') 
d=read.table(zz,header=T)

# compute the sum of heterozygous sites (GT==1) per sample per window
# here, done in 250k windows
#
ww=250000
dc = c()
for(ctg in f[,1]){
    d1=d[as.character(d[,1])==ctg,]

    d2=cbind(d1,floor(d1$POS/ww)+1)
    colnames(d2)[dim(d2)[2]]="window"

    if(dim(d2)[1]>2){
        
        # H for a single sample
        d3=aggregate(d2[,3:dim(d)[2]],by=list(d2$window),FUN=function(x){sum(x==1)})
        colnames(d3)[1]="window"

        m1=m[m[,1]==ctg,]
        m2=IRanges(start=m1$start,end=m1$end)

        w1=seq(1,f$V2[f$V1==ctg],ww)
        w2=IRanges(start=w1,width=ww)

        s=rep(0,dim(d3)[1])
        t=rep(0,dim(d3)[1])
        for(i in end(w2)/ww){
            s[i]=sum(width( restrict(m2,start(w2)[i],end(w2)[i]) ))
            t[i]=i
        }
        s = cbind(t,s)
        colnames(s)=c("window","nsites")

        d4=merge(d3,s,by.x = 1, by.y = 1)

        dc=rbind.data.frame(dc,cbind(ctg,d4))
    }
}

# the file contains the sums of heterozygous sites per sample per window
# 1st column: contig
# 2nd column: window number
# columns 3-100: 98 Saimaa samples
# 103rd column: number of positive-masked sites in the window
#
write.table(dc,paste0("data/saimaa_windows_posmask_",ww,".tsv"),row.names=F,quote=F)


#zz=gzfile("data/saimaa_all_posmask_98.tsv.gz",'rt') 
#d=read.table(zz,header=T)

ww=250000
dc = c()
for(ctg in f[,1]){
    d1=d[as.character(d[,1])==ctg,]

    d2=cbind(d1,floor(d1$POS/ww)+1)
    colnames(d2)[dim(d2)[2]]="window"

    if(dim(d2)[1]>2){

        ns.c=match(ns,gsub("X","",colnames(d)))
        cs.c=match(cs,gsub("X","",colnames(d)))
        ss.c=match(ss,gsub("X","",colnames(d)))
        as.c=match(as,gsub("X","",colnames(d)))
        
        # pi for a single population
        set1=combinations(ns.c,2)
        set2=combinations(cs.c,2)
        set3=combinations(ss.c,2)
        pis=cbind.data.frame(ns=pi_for_four(set1,d2),
                             cs=pi_for_four(set2,d2),
                             ss=pi_for_four(set3,d2))
        
        # dxy for two populations
        set1=expand.grid(ns.c,cs.c)
        set2=expand.grid(ns.c,ss.c)
        set3=expand.grid(cs.c,ss.c)
        pis=cbind.data.frame(pis,
                             ns.cs=pi_for_four(set1,d2),
                             ns.ss=pi_for_four(set2,d2),
                             cs.ss=pi_for_four(set3,d2))

        d3=aggregate(pis,by=list(d2$window),FUN=sum)
        colnames(d3)[1]="window"

        m1=m[m[,1]==ctg,]
        m2=IRanges(start=m1$start,end=m1$end)

        w1=seq(1,f$V2[f$V1==ctg],ww)
        w2=IRanges(start=w1,width=ww)

        s=rep(0,dim(d3)[1])
        t=rep(0,dim(d3)[1])
        for(i in end(w2)/ww){
            s[i]=sum(width( restrict(m2,start(w2)[i],end(w2)[i]) ))
            t[i]=i
        }
        s = cbind(t,s)
        colnames(s)=c("window","nsites")

        d4=merge(d3,s,by.x = 1, by.y = 1)

        dc=rbind.data.frame(dc,cbind(ctg,d4))
    }
}
write.table(dc,paste0("data/3subpops_saimaa_windows_posmask_",ww,".tsv"),row.names=F,quote=F)

ggplot_pop1 <- function(dat) {
    ggplot() +
    geom_line(data = dat,aes(x,y,group=smp,col=pop),size=0.25,alpha=0.5) +
    xlim(c(-0.001,0.01)) + 
    xlab("Pi") + ylab("") + 
    theme(axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size = 0.25),
          axis.ticks.x = element_line(size = 0.25),
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=5),
          text = element_text(size=8),
         )
}

ggplot_pop2 <- function(dat) {
    ggplot() +
    geom_line(data = dat,aes(x,y,group=pop,col=pop),size=0.5) +
    xlim(c(-0.001,0.01)) + 
    xlab("Pi") + ylab("") + 
    theme(axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size = 0.25),
          axis.ticks.x = element_line(size = 0.25),
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=5),
          text = element_text(size=8),
         )
}

ggplot_pop3 <- function(dat) {
    dat$npop=rep(1,dim(dat)[1])
    dat$npop[dat$pop %in% c("ns.cs","ns.ss","cs.ss")]=2
    ggplot() +
    geom_line(data = subset(dat,npop==2),aes(x,y,group=pop,col=pop),linetype="solid",size=0.5) +
    geom_line(data = subset(dat,npop==1),aes(x,y,group=pop,col=pop),linetype="dashed",size=0.5) +
    xlim(c(-0.001,0.006)) + 
    xlab("Pi") + ylab("") + 
    theme(axis.line.y = element_blank(), 
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size = 0.25),
          axis.ticks.x = element_line(size = 0.25),
          axis.text.y = element_text(size=5),
          axis.text.x = element_text(size=5),
          text = element_text(size=8),
         )
}


theta_density <- function(dat,smp,pop) {
    d=density(dat$pi,na.rm=T,from=0,adj=0.5)
    d$x=c(min(d$x),d$x,max(d$x))
    d$y=c(0,d$y,0)
    data.frame(smp=rep(smp,length(d$x)),
               pop=rep(pop,length(d$x)),
               mean=rep(mean(dat$pi,na.rm=T),length(d$x)),
               x=d$x,y=d$y)
}

theta_sort <- function(dat) {
    tm = aggregate(dat$mean,list(dat$smp),FUN=mean)
    tn = aggregate(rep(1,dim(dat)[1]),list(dat$smp),FUN=sum)
    dat=cbind(dat,ord=rep(rank(tm$x),tn$x))
}


theta_pop <- function(dat,pop) {
    res=c()
    for(smp in unique(dat$smp)){
        res=rbind(res,theta_density(dat[dat$smp==smp,],smp,pop))
    }
    dat=theta_sort(res)
    dat
}

ww=250000
data250k.gt=read.table(paste0("data/saimaa_windows_posmask_",ww,".tsv"),head=T)

data250k.gt2=data250k.gt
data250k.gt2[,3:100]=data250k.gt2[,3:100]/data250k.gt2[,101]
data250k.gt2=melt(data250k.gt2,id.vars=c(1:2,101))
colnames(data250k.gt2)=c("ctg","window","nsites","smp","pi")

densn250k.gt=theta_pop(data250k.gt2[data250k.gt2$nsites>150000,],"saimaa")

set1=densn250k.gt[densn250k.gt$smp %in% nsx,]
set1$pop="ns"
set2=densn250k.gt[densn250k.gt$smp %in% csx,]
set2$pop="cs"
set3=densn250k.gt[densn250k.gt$smp %in% ssx,]
set3$pop="ss"
setA250=rbind.data.frame(set1,set2,set3)

# H for all but one sample
# (sample 2572 excluded as it is nearly completely homozygous for contig 0)
#
plot1 <- ggplot_pop1(setA250[setA250$smp!="X2572",])
plot1 <- plot1+
    coord_cartesian(xlim=c(-0.0001,0.0075))+
    theme_classic()+
    theme(
        text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position="top",
        legend.key.width = unit(3,"mm"),
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.text=element_text(size=10,margin = margin(r = 8, l=-2, unit = "pt")))+
    guides(color=guide_legend(title="",nrow = 1,override.aes=list(size=3,alpha=1)))+
    ylab("")+xlab("Nucleotide diversity")+
    scale_color_manual(breaks=c("ns", "cs", "ss", "ns.cs", "ns.ss", "cs.ss"),
                       labels=c("North", "Central", "South Saimaa", "NS-CS", "NS-SS", "CS-SS"),
                       values=col3)

plot1

ww=250000
pi250k.3p=read.table(paste0("data/3subpops_saimaa_windows_posmask_",ww,".tsv"),head=T)
pi250k.3p[,3:8]=pi250k.3p[,3:8]/pi250k.3p[,9]

pi250k.3p2 = melt(pi250k.3p,id.vars = c(1,2,9))
colnames(pi250k.3p2)=c("ctg","window","nsites","smp","pi")

pop="saimaa"
densn250k.3p=theta_pop(pi250k.3p2[pi250k.3p2$nsites>150000,],pop)
densn250k.3p$pop=densn250k.3p$smp

# plot within-pop pi only 
#
plot2 <- ggplot_pop2(densn250k.3p[densn250k.3p$smp %in% c("ns","cs","ss"),])

plot2 <- plot2 +
    coord_cartesian(xlim=c(-0.0001,0.0075))+
    theme_classic()+
    theme(
        text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position="top",
        legend.key.width = unit(3,"mm"),
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.text=element_text(size=10,margin = margin(r = 8, l=-4, unit = "pt")))+
    guides(color=guide_legend(title="",nrow = 1,override.aes=list(size=3)))+
    ylim(0,700)+
    ylab("")+xlab("Nucleotide divergence")+
    scale_color_manual(breaks=c("ns", "cs", "ss", "ns.cs", "ns.ss", "cs.ss"),
                       labels=c("North", "Central", "South Saimaa", "NS-CS", "NS-SS", "CS-SS"),
                       values=c(col3,col6[2:4]))

plot2

# plot within pop pi and across-pop dxy
#
plot3 <- ggplot_pop3(densn250k.3p)

plot3 <- plot3 +
    coord_cartesian(xlim=c(-0.0001,0.0075))+
    guides(color=guide_legend(title="",nrow = 1,override.aes=list(size=2)))+
    theme_classic()+
    theme(
        text = element_text(size=8),
        axis.text = element_text(size=7),
        legend.position="top",
        legend.key.width = unit(1.5,"mm"),
        axis.line.y = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.margin=margin(t = 0, unit='cm'),
        legend.text=element_text(size=6,margin = margin(r = 2, l=-4, unit = "pt")))+
    ylim(0,700)+
    ylab("")+xlab("Nucleotide diversity")+
    scale_color_manual(breaks=c("ns", "cs", "ss", "ns.cs", "ns.ss", "cs.ss"),
                       labels=c("NS", "CS", "SS", "NS-CS", "NS-SS", "CS-SS"),
                       values=c(col3,col6[2:4]))

plot3
