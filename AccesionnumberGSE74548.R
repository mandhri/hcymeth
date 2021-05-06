## ----packages-----------------------------------------------------------------
suppressPackageStartupMessages({
    library("plyr")
    library("R.utils")
    library("missMethyl")
    library("limma")
    library("topconfects")
    library("minfi")
    library("IlluminaHumanMethylation450kmanifest")
    library("RColorBrewer")
    library("IlluminaHumanMethylation450kanno.ilmn12.hg19")
    library("eulerr")
    library("plyr")
    library("gplots")
    library("reshape2")
    library("beeswarm")
  })
# Annotation
ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
myann <- data.frame(ann450k[,c("UCSC_RefGene_Name","Regulatory_Feature_Group")])
promoters <- grep("Prom",myann$Regulatory_Feature_Group)


## ----load_functions-----------------------------------------------------------
source("https://raw.githubusercontent.com/markziemann/ART_methylation/master/meth_functions.R")
myranks<-function(x) {
  x$score <- sign(x$logFC)/log10(x$adj.P.Val)
  y <- x[,"score",drop=FALSE]
  y$rn <- x$Row.names
  return(y)
}

# heatmap for continuous data
heatmap_c<-function(dm,name,mx,n, groups) {
  my_palette <- colorRampPalette(c("blue", "white", "red"))(25)
  topgenes <-  rownames(head(dm[order(dm$P.Value),],n))
  ss <- mx[which(rownames(mx) %in% topgenes),]
  mycols <- colorRampPalette(c("white","yellow", "orange", "red","darkred"))(n = length(groups))
  colCols <- mycols[order(groups)]
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.6,trace="none",cexCol=0.4,
  ColSideColors=colCols ,  col=my_palette, main=name)
}

# heatmap for continuous data - topconfects
make_heatmap2_c <- function(confects,name,mx,n, groups) {
  topgenes <-  head(confects$table$name,n)
  my_palette <- colorRampPalette(c("blue", "white", "red"))(25)
  ss <- mx[which(rownames(mx) %in% topgenes),]
  mycols <- colorRampPalette(c("white","yellow", "orange", "red","darkred"))(n = length(groups))
  colCols <- mycols[order(groups)]
  heatmap.2(ss,scale="row",margin=c(10, 10),cexRow=0.6,trace="none",cexCol=0.4,
  ColSideColors=colCols ,  col=my_palette, main=name)  
}



## ----data---------------------------------------------------------------------

dir.create("GSE74548")

ARRAY_SAMPLESHEET="GSE74548/GSE74548_Sample_description_BProof.txt"
# only download it if it is not present on the system
if ( !file.exists(ARRAY_SAMPLESHEET ) ) {
    DLFILE=paste(ARRAY_SAMPLESHEET,".gz",sep="")
    download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74548/suppl/GSE74548_Sample_description_BProof.txt.gz",
        destfile = DLFILE)
    gunzip(DLFILE)
}

ARRAY_DATA="GSE74548/GSE74548_RAW.tar"
# only download it if it is not present on the system
if ( !dir.exists("GSE74548/IDAT") ) {
  dir.create("GSE74548/IDAT")
  download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE74548&format=file",
    destfile = ARRAY_DATA)
    untar(exdir = "GSE74548/IDAT", tarfile = ARRAY_DATA)
}

baseDir <- "GSE74548"
R.utils::gunzip("GSE74548/IDAT/GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz",overwrite=TRUE, remove=FALSE)
targets <- read.metharray.sheet(baseDir,pattern="csv")
download.file("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE74nnn/GSE74548/matrix/GSE74548_series_matrix.txt.gz",destfile = "GSE74548_series_matrix.txt.gz")

gse<- getGEO(filename = "GSE74548_series_matrix.txt.gz")

targets <- pData(phenoData(gse))
targets <- targets[order(rownames(targets)),]
mybase <- unique(gsub("_Red.idat.gz" ,"", gsub("_Grn.idat.gz", "" ,list.files("./GSE74548",pattern = "GSM",recursive = TRUE))))
mybase <- paste("GSE74548/", mybase, sep = "")
# sample number discrepancy. 311 IDAT FILES but GEO states 174 subjects.
gsm <-sapply(strsplit(mybase,"_"),"[[",1)
gsm <-gsub("GSE74548/IDAT/","",gsm)
#only using sampels described in meta data
targets$Basename<- mybase[which(gsm %in% rownames(targets))]
rgSet <- read.metharray.exp(targets = targets)
rgSet
head(targets)



## ----norm,fig.cap="Figure 1. Normalisation of bead-array data with SWAN."-----
mSet <- preprocessRaw(rgSet)
mSetSw <- SWAN(mSet,verbose=FALSE)
par(mfrow=c(1,2), cex=0.8)
densityByProbeType(mSet[,1], main = "Raw")
densityByProbeType(mSetSw[,1], main = "SWAN")


## ----filterprobes-------------------------------------------------------------
# include sex chromosomes
detP <- detectionP(rgSet)
keep <- rowSums(detP < 0.01) == ncol(rgSet)
mSetSw <- mSetSw[keep,]

# exclude SNP probes
mSetSw <- mapToGenome(mSetSw)
mSetSw_nosnp <- dropLociWithSnps(mSetSw)
dim(mSetSw)
dim(mSetSw_nosnp)
mSetSw <- mSetSw_nosnp

# exclude sex chromosomes
keep <- !(featureNames(mSetSw) %in% ann450k$Name[ann450k$chr %in% c("chrX","chrY")])
mSetFlt <- mSetSw[keep,]
head(mSetFlt)
dim(mSetFlt)


## ----beta_m_vals--------------------------------------------------------------
# include sex chromosomes
meth <- getMeth(mSetSw)
unmeth <- getUnmeth(mSetSw)
Mval <- log2((meth + 100)/(unmeth + 100))
beta <- getBeta(mSetSw)

# exclude sex chromosomes
meth <- getMeth(mSetFlt)
unmeth <- getUnmeth(mSetFlt)
Mval_flt <- log2((meth + 100)/(unmeth + 100))
beta_flt <- getBeta(mSetFlt)


## ----sex chromosome diagnostic------------------------------------------------
cgx<-rownames(Locations[which(Locations$chr %in% "chrX"),])
cgy<-rownames(Locations[which(Locations$chr %in% "chrY"),])
mvx<-Mval[which(rownames(Mval) %in% cgx),]
mvy<-Mval[which(rownames(Mval) %in% cgy),]


targets_m<-rownames(subset(targets,`gender:ch1`=="male"))
targets_f<-rownames(subset(targets,`gender:ch1`=="female"))


Mvalm<-Mval[,colnames(Mval)%in%targets_m]
Mvalf<-Mval[,colnames(Mval)%in%targets_f]

Mvalf
dim(Mvalf)

head(Mvalm)
mvxm<-Mvalm[which(rownames(Mvalm) %in% cgx),]
mvym<-Mvalm[which(rownames(Mvalm) %in% cgy),]
mvxm

mvxf<-Mvalf[which(rownames(Mvalf) %in% cgx),]
mvyf<-Mvalf[which(rownames(Mvalf) %in% cgy),]

plot(colMeans(mvx),colMeans(mvy),col="gray")
points(colMeans(mvxf),colMeans(mvyf),col="red")
points(colMeans(mvxm),colMeans(mvym),col="blue")



## ----scree,fig.cap="Figure 2. Screeplot shows contribution of top principal components to the overall variation in the experiment."----
colnames(Mval)<-sapply(strsplit(colnames(Mval),"_"),"[[",1)
colnames(Mval_flt)<-sapply(strsplit(colnames(Mval_flt),"_"),"[[",1)
head(Mval)
head(Mval_flt)
par(mfrow=c(2,1))  
myscree(Mval,main="incl sex chr")
myscree(Mval_flt,main="excl sex chr")


## ----mds1,fig.width = 8 ,fig.height = 8,fig.cap="Figure 3. MDS plot coloured by targets groups"----
head(targets)
dim(targets)
targets$group <- factor(targets$source_name_ch1)
sample_group<-factor(targets$group)
targets$sex <- factor(targets$`gender:ch1`)
dim(sample_group)

colour_palette=brewer.pal(n = length(levels(targets$group)), name = "Paired")
colours <- colour_palette[as.integer(factor(targets$source_name_ch1))]
plot(1,axes = FALSE,xlab="",ylab="",main="treatment groups")
legend("center",legend=levels(targets$group),pch=16,cex=1.2,col=colour_palette)
mydist <- plotMDS(Mval, labels=targets$geo_accession,col=colours,main="sex chromosomes included")
mydist_flt <- plotMDS(Mval_flt, labels=targets$geo_accession,col=colours,main="sex chromosomes excluded")


## ----hist1--------------------------------------------------------------------
hcy<- targets$characteristics_ch1.10
hcy<-strsplit(as.character(targets$characteristics_ch1.10) , " ")
hcy<-sapply(hcy , "[" ,5)
hcy<-as.numeric(hcy)
hist(hcy,breaks = 20,xlab = "hcy levels")
hcy_groups<-cut(hcy,breaks = c(0,10,16,30),labels = c("low","medium","high"))
table(hcy_groups)


## ----mds2,fig.width = 8 ,fig.height = 8,fig.cap="Figure 3. MDS plot coloured by homocysteine levels"----
head(targets)
targets$sex <- factor(targets$`gender:ch1`)
targets$hcy_groups <- hcy_groups
sample_group<-hcy_groups

colour_palette=brewer.pal(n = length(levels(sample_group)), name = "Paired")
colours <- colour_palette[as.integer(factor(hcy_groups))]
plot(1,axes = FALSE,xlab="",ylab="",main="hcy levels")
legend("center",legend=levels(hcy_groups),pch=16,cex=1.2,col=colour_palette)
plotMDS(mydist, labels=targets$geo_accession,col=colours,main="sex chromosomes included")
plotMDS(mydist_flt, labels=targets$geo_accession,col=colours,main="sex chromosomes excluded")


## ----hist2--------------------------------------------------------------------
age<-targets$characteristics_ch1.5
age<-strsplit(as.character(targets$characteristics_ch1.5) , " ")
age<-sapply(age, "[" ,4) 
age<-as.numeric(age)
hist(age,breaks = 20,xlab = "age levels")
age_groups<-cut(age,breaks = c(0,65,70,75),labels = c("old","older","oldest"))
table(age_groups)


## ----mds3,fig.width = 8 ,fig.height = 8,fig.cap="Figure 3. MDS plot coloured by age"----
head(targets)
targets$age <- factor(age_groups)
sample_group <-factor(age_groups)
colour_palette=brewer.pal(n = length(levels(sample_group)), name = "Paired")
colours <- colour_palette[as.integer(factor(age_groups))]
plot(1,axes = FALSE,xlab="",ylab="",main="age")
legend("center",legend=levels(age_groups),pch=16,cex=1.2,col=colour_palette)
plotMDS(mydist, labels=targets$geo_accession,col=colours,main="sex chromosomes included")
plotMDS(mydist_flt,labels=targets$geo_accession,col=colours,main="sex chromosomes excluded")


## ----hist3--------------------------------------------------------------------
folate_levels<-targets$characteristics_ch1.8
folate_levels<-strsplit(as.character(targets$characteristics_ch1.8) , " ")
folate_levels<-sapply(folate_levels, "[" ,5)
folate_levels<-as.numeric(folate_levels)
hist(folate_levels,breaks = 20,xlab = "folate levels")
folate_levels<-cut(folate_levels,breaks = c(0,20,35,50,70,90),labels = c("lowest","low","medium","high","highest"))
table(folate_levels)
folate_levels


## ----mds4,fig.width = 8 ,fig.height = 8,fig.cap="Figure 3. MDS plot coloured by folate levels"----
targets$folate_levels <- folate_levels
sample_group<-folate_levels

colour_palette=brewer.pal(n = length(levels(sample_group)), name = "Paired")
colours <- colour_palette[as.integer(factor(targets$folate_levels))]
plot(1,axes = FALSE,xlab="",ylab="",main="Folate levels")
legend("center",legend=levels(targets$folate_levels),pch=16,cex=1.2,col=colour_palette)
plotMDS(mydist, labels=targets$geo_accession,col=colours,main="sex chromosomes included")
plotMDS(mydist_flt, labels=targets$geo_accession,col=colours,main="sex chromosomes excluded")


## ----hist4--------------------------------------------------------------------
vitb12_levels<-targets$characteristics_ch1.9
vitb12_levels<-strsplit(as.character(targets$characteristics_ch1.9) , " ")
vitb12_levels<-sapply(vitb12_levels ,"[",6)
vitb12_levels<-as.numeric(vitb12_levels)
hist(vitb12_levels,breaks = 20,xlab = "vitb12 levels")
vitb12_levels<-cut(vitb12_levels,breaks = c(100,300,500,700,900,1115),labels = c("lowest","low","medium","high","highest"))
table(vitb12_levels)
vitb12_levels


## ----mds5,fig.width = 8 ,fig.height = 8,fig.cap="Figure 3. MDS plot coloured by vitaminb12 levels"----
targets$vitb12_levels <- vitb12_levels
sample_group<-vitb12_levels

colour_palette=brewer.pal(n = length(levels(sample_group)), name = "Paired")
colours <- colour_palette[as.integer(factor(targets$vitb12_levels))]
plot(1,axes = FALSE,xlab="",ylab="",main="vitb12 levels")
legend("center",legend=levels(targets$vitb12_levels),pch=16,cex=1.2,col=colour_palette)
plotMDS(mydist, labels=targets$geo_accession,col=colours,main="sex chromosomes included")
plotMDS(mydist_flt, labels=targets$geo_accession,col=colours,main="sex chromosomes excluded")


## ----dm1A,fig.height=8,fig.width=8--------------------------------------------
samplesheet<-targets[grep("follow-up",targets$source_name_ch1),]
sex <- factor(samplesheet$`gender:ch1`)

age<-strsplit(as.character(samplesheet$characteristics_ch1.5) , " ")
age<-sapply(age,"[",4) 
age<-as.numeric(age)
groups<-factor(samplesheet$source_name_ch1,levels = c("Buffy coat, placebo, follow-up","Buffy coat, FA/vB12, follow-up"))
mx <-Mval_flt
name="placebo_vs_supplement_follow-up"
head(myann)
head(samplesheet)

design <- model.matrix(~age+ sex +groups)
    mxs <- mx[,which(colnames(mx) %in% rownames(samplesheet) )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=4, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma1a <- dma[order(dma$P.Value),]
    head(dm)
    head(dma1a)
mxs
dim(design)
dma1a_d<-nrow(subset(dm,adj.P.Val<0.05,logFC<0))
dma1a_u<-nrow(subset(dm,adj.P.Val<0.05,logFC>0))

confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
head(confects)

colCols <- as.numeric(as.factor(groups))

colCols

make_volcano(dma1a,name = "placebo_vs_supplement_follow-up(sex,age)",mx=Mval_flt)
rownames(dma1a)<-dma1a[,1]
colnames(dma1a)
head(dma1a)
make_beeswarms(dm=dma1a ,name="placebo_vs_supplement_follow-up(sex,age)" , mx=beta_flt , groups=groups , n= 15)
make_beeswarms_confects(confects = confects ,name="placebo_vs_supplement_follow-up(sex,age)" , mx=beta_flt , groups=groups , n= 15)

make_heatmap(dm=dma1a , name="placebo_vs_supplement_baseline(sex,age)" , mx=mxs ,n = 50, groups=groups)
make_heatmap2(confects = confects , name="placebo_vs_supplement_baseline(sex,age)",mx=mxs ,n = 50, groups=groups)


## ----dm2A---------------------------------------------------------------------
samplesheet<-targets[grep("baseline",targets$source_name_ch1),]
sex <- factor(samplesheet$`gender:ch1`)
age<-strsplit(samplesheet$characteristics_ch1.5)," ")
age<s-sapply(age,"[",4)
age<-as.numeric(age)
groups<-factor(samplesheet$source_name_ch1,levels = c("Buffy coat, placebo, baseline","Buffy coat, FA/vB12, baseline"))
mx <-Mval_flt
name="placebo_vs_supplement_baseline"
head(myann)
head(samplesheet)
design <- model.matrix(~age+ sex +groups)
    mxs <- mx[,which(colnames(mx) %in% rownames(samplesheet) )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=4, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma2a <- dma[order(dma$P.Value),]
    head(dm)
    head(dma2a)
    dma2a
str(mxs)
dim(design)

dma2a_d<-nrow(subset(dm,adj.P.Val<0.05,logFC<0))
dma2a_u<-nrow(subset(dm,adj.P.Val<0.05,logFC>0))


confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
head(confects)

make_volcano(dma2a,name = "placebo_vs_supplement_baseline(sex,age)",mx=Mval_flt)
rownames(dma2a)<-dma2a[,1]
str(dma2a)
make_beeswarms(dm=dma2a ,name="placebo_vs_supplement_baseline(sex,age)" , mx=beta_flt , groups=groups , n= 15)
make_beeswarms_confects(confects = confects ,name="placebo_vs_supplement_baseline(sex,age)" , mx=beta_flt , groups=groups , n= 15)

make_heatmap(dm=dma2a , name="plb_supl_baseline" , mx=mxs ,n = 50, groups=groups)
make_heatmap2(confects = confects , name="placebo_vs_supplement_baseline(sex,age)",mx=mxs ,n = 50, groups=groups)


## ----dm3A,fig.height=8,fig.width=8--------------------------------------------
samplesheet<-targets
hcy<- samplesheet$characteristics_ch1.10
hcy<-strsplit(as.character(samplesheet@characteristics_ch1.10)," ")
hcy<-sapply(hcy,"[",5)
hcy<- as.numeric(hcy)
head(samplesheet)
samplesheet<-samplesheet[which(!is.na(hcy)),]
hcy<-hcy[which(!is.na(hcy))]
groups<-cut(hcy,breaks = c(0,10,16,30),labels = c("low","medium","high"))
groups<-as.integer(groups)

sex <- factor(samplesheet$`gender:ch1`)
head(samplesheet)

age<-strsplit(as.character(samplesheet@characteristics_ch1.5), " ")
age<-sapply(age,"[",4)
age<-as.numeric(age)
mx <-Mval_flt
name="hcy_levels "
head(myann)
head(mx)

design <- model.matrix(~ age+ sex + groups)
mxs <- mx[,which( colnames(mx) %in% rownames(samplesheet) )]
fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=4, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma3a <- dma[order(dma$P.Value),]
    head(dm)
    head(dma3a)
 dma3a_d<-nrow(subset(dm,adj.P.Val<0.05,logFC<0))
 dma3a_u<-nrow(subset(dm,adj.P.Val<0.05,logFC>0)) 
 
confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
head(confects)

make_volcano(dma3a,name = "hcy_levels(sex,age)",mx=Mval_flt)
rownames(dma3a)<-dma3a[,1]
make_beeswarms(dm=dma3a ,name="hcy_levels(sex,age)" , mx=beta_flt , groups=groups , n= 15)
make_beeswarms_confects(confects = confects ,name="hcy_levels(sex,age)" , mx=beta_flt , groups=groups , n= 15)

heatmap_c(dm=dma3a , name="hcy_levels(sex,age)" , mx=mxs ,n = 50, groups=groups)
make_heatmap2_c(confects = confects , name="hcy_levels(sex,age)",mx=mxs ,n = 50, groups=groups)


## ----dm4A---------------------------------------------------------------------
samplesheet<-targets
folate_levels<- samplesheet$characteristics_ch1.8
folate_levels<-strsplit(as.character(samplesheet$characteristics_ch1.8)," ")
folate_levels<-sapply(folate_levels,"[",5)
folate_levels<- as.numeric(folate_levels)
head(samplesheet)

sex <- factor(samplesheet$`gender:ch1`)
head(samplesheet)
dim(sex)
groups<-factor(samplesheet$source_name_ch1,levels = c("Buffy coat, placebo, baseline","Buffy coat, FA/vB12, baseline"))
mx <-Mval_flt
name="folate_levels"
head(myann)

design <- model.matrix(~ age+ sex + folate_levels)
    mxs <- mx[,which( colnames(mx) %in% rownames(samplesheet) )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=4, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma4a <- dma[order(dma$P.Value),]
    head(dm)
    head(dma4a)
    dma4a_d<-nrow(subset(dm,adj.P.Val<0.05,logFC<0))
    dma4a_u<-nrow(subset(dm,adj.P.Val<0.05,logFC>0))

confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
head(confects)

make_volcano(dma4a,name = "folate_levels(sex,age)",mx=Mval_flt)
rownames(dma4a)<-dma4a[,1]
make_beeswarms(dm=dma4a ,name="folate_levels(sex,age)" , mx=beta_flt , groups=groups , n= 15)
make_beeswarms_confects(confects = confects ,name="folate_levels(sex,age)" , mx=beta_flt , groups=groups , n= 15)

make_heatmap(dm=dma4a , name="folate_levels(sex,age)" , mx=mxs ,n = 50, groups=groups)
make_heatmap2(confects = confects , name="folate_levels(sex,age)",mx=mxs ,n = 50, groups=groups)


## ----dm5A---------------------------------------------------------------------
samplesheet<-targets
vitminb12<- samplesheet$characteristics_ch1.9
vitminb12<-strsplit(as.character(samplesheet$characteristics_ch1.9),""
vitminb12<-sapply(vitaminb12b12,"[",6)
vitminb12<- as.numeric(vitminb12)
head(samplesheet)

sex <- factor(samplesheet$`gender:ch1`)
head(samplesheet)
mx <-Mval_flt
name="vitmin_b12_levels"
head(myann)

groups<-factor(samplesheet$source_name_ch1,levels = c("Buffy coat, placebo, baseline","Buffy coat, FA/vB12, baseline"))
design <- model.matrix(~ age+ sex + vitminb12)
    mxs <- mx[,which( colnames(mx) %in% rownames(samplesheet) )]
    fit.reduced <- lmFit(mxs,design)
    fit.reduced <- eBayes(fit.reduced)
    summary(decideTests(fit.reduced))
    dm <- topTable(fit.reduced,coef=4, number = Inf)
    dma <- merge(myann,dm,by=0)
    dma5a <- dma[order(dma$P.Value),]
    head(dm)
    head(dma5a)
    dma5a_d<- nrow(subset(dm,adj.P.Val<0.05,logFC<0))
    dma5a_u<-nrow(subset(dm,adj.P.Val<0.05,logFC>0))
 
confects <- limma_confects(fit.reduced, coef=3, fdr=0.05)
head(confects)

make_volcano(dma5a,name = "vitmin_b12_levels(sex,age)",mx=Mval_flt)
rownames(dma5a)<-dma5a[,1]
make_beeswarms(dm=dma5a ,name="vitmin_b12_levels" , mx=beta_flt , groups=groups , n= 15)
make_beeswarms_confects(confects = confects ,name="vitmin_b12_levels" , mx=beta_flt , groups=groups , n= 15)

make_heatmap(dm=dma5a , name="vitmin_b12_levels(sex,age)" , mx=mxs ,n = 50, groups=groups)
make_heatmap2(confects = confects , name="vitmin_b12_levels(sex,age)",mx=mxs ,n = 50, groups=groups)


## ----summarised dms-----------------------------------------------------------

downs<-c(dma1a_d,dma2a_d,dma3a_d, dma4a_d,dma5a_d)
ups<-c(dma1a_u,dma2a_u, dma3a_u,dma4a_u,dma5a_u)
my_summary<-data.frame(downs,ups)
row.names(my_summary)<-c("placebo vs supplement follow-up(sex,age)",
                         "placebo vs supplement Baseline(sex,age)",
                         "hcy levels(sex,age)",
                         "folate levels(sex,age)",
                         "vitb12 levels(sex,age)")
library(knitr)
kable(my_summary)


## ----venn1,fig.cap="Comparison of probes in related to placebo vs supplement in follow-up(sex,age)", fig.width = 8 ,fig.height = 8----
v1 <- list("plb v sup fl(up)" =dma1a_u , 
           "hcy(up)" =dma3a_d ,
           "plb v sup fl(dwn)" =dma1a_d ,
           "hcy(dwn)" =dma3a_d)
head(v1)
plot(euler(v1, shape = "ellipse"), quantities = TRUE)


## ----venn2,fig.cap="Comparison of probes in related to hcy(sex) vs folate level(sex,age)", fig.width = 8 ,fig.height = 8----
v2 <- list("fol(up)" =dma4a_u, 
           "hcy(up)" =dma3a_u,
           "fol(dwn)" =dma4a_d,
           "hcy(dwn)" =dma3a_d)
head(v2)
str(v2)
plot(euler(v2, shape = "ellipse"), quantities = TRUE)



## ----venn3,fig.cap="Comparison of probes in related to hcy(age) vs folate level(sex,age)", fig.width = 8 ,fig.height = 8----
v3 <- list("plb v sup fl(up)" =dma1a_u, 
           "fol(up)" = dma4a_u ,
           "plb v sup fl(dwn)" =dma1a_d ,
           "fol(dwn)" =dma4a_d)
head(v3)
plot(euler(v3, shape = "ellipse"), quantities = TRUE)


## ----venn4,fig.cap="Comparison of probes in related to hcy vs folate level(sex,age)", fig.width = 8 ,fig.height = 8----
v4 <- list("hcy(up)" =dma3a_u, 
           "fol(up)" = dma4a_u ,
           "hcy(dwn)" =dma3a_d ,
           "fol(dwn)" =dma4a_d)
head(v4)
plot(euler(v4, shape = "ellipse"), quantities = TRUE)


## ----venn5,fig.cap="Comparison of probes in related to hcy vs vitb12 level(sex,age)", fig.width = 8 ,fig.height = 8----
v5<- list("hcy(up)" =dma3a_u, 
           "vitb12(up)" = dma5a_u ,
           "hcy(dwn)" =dma3a_d ,
           "vitb12(dwn)" =dma5a_d)
head(v5)
plot(euler(v5, shape = "ellipse"), quantities = TRUE)


## ----spearman, fig.width = 8 ,fig.height = 8----------------------------------
mycontrasts <- list("placebo_vs_supplement-follow-up(sex,age)"=dma1a,"placebo vs supplement Baseline(sex,age)"=dma2a,"hcy levels(sex,age)"=dma3a,"folate levels(sex,age)"=dma4a,"vitb12 levels(sex,age)"=dma5a)
lapply(mycontrasts, head)
myrnks<-lapply(X = mycontrasts, FUN = myranks)
str(myrnks)
df <- join_all(myrnks,by="rn")

rownames(df) <- df$rn
df$rn=NULL
colnames(df) <- names(mycontrasts)
head(df)
mycors <- cor(df,method = "spearman")
my_palette <- colorRampPalette(c("darkred","red", "orange", "yellow","white"))(n = 25)
heatmap.2(mycors,scale="none",margin=c(10, 10),cexRow=0.8,trace="none",cexCol=0.8,
    col=my_palette,main="Spearman correlations")
mycors


## ----genesets, fig.width = 8 ,fig.height = 8----------------------------------
library("mitch")
# gene sets
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip", 
    destfile="ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip",overwrite = TRUE)
genesets <- gmt_import("ReactomePathways.gmt")


## ----1danalysis, fig.width = 8 ,fig.height = 8--------------------------------
placebo_vs_supplement_followup_sex_age_mitch <- run_mitch_1d(dma= dma1a, name="placebo_vs_supplement_followup(sex,age)_mitch")
head(placebo_vs_supplement_followup_sex_age_mitch,50)
placebo_vs_supplement_Baseline_sex_age_mitch <- run_mitch_1d(dma= dma2a, name="placebo_vs_supplement_Baseline_(sex,age)_mitch")
head(placebo_vs_supplement_Baseline_sex_age_mitch,50)
hcy_levels_sex_age_mitch <- run_mitch_1d(dma= dma3a, name="hcy_levels(sex,age)_mitch")
head(hcy_levels_sex_age_mitch,50)
folate_levels_sex_age_mitch <- run_mitch_1d(dma= dma4a, name="folate_levels(sex,age)_mitch")
head(folate_levels_sex_age_mitch,50)
vitb12_levels_sex_age_mitch <- run_mitch_1d(dma= dma5a, name="vitb12 levels(sex,age)_mitch")
head(vitb12_levels_sex_age_mitch,50)


## ----mitch, fig.width = 8 ,fig.height = 8-------------------------------------
xl <- list("placebo vs supplement follow-up(sex,age)"=dma1a,"placebo vs supplement Baseline(sex,age)"=dma2a, 
"hcy_levels(sex_age)"=dma3a,"folate_levels(sex_age)"=dma4a,"vitb12_levels(sex_age)"= dma5a)
xxl <- lapply(X = xl,run_mitch_rank)  
xxxl <- lapply(xxl,function(xxl) { xxl$genenames <- rownames(xxl) ; xxl} )
xxll <- join_all(xxxl,by="genenames")
rownames(xxll) <- xxll$genenames
xxll$genenames=NULL
colnames(xxll) <- names(xl)
head(xxll)
capture.output(
        res <- mitch_calc(xxll,genesets = genesets,priority = "significance")
        , file = "/dev/null", append = FALSE,
        type = c("output", "message"), split = FALSE)
head(res$enrichment_result,20)
unlink("estill_multi_mitch.pdf")
     capture.output(
        mitch_plots(res,outfile="estill_multi_mitch.pdf")
        , file = "/dev/null", append = FALSE,
        type = c("output", "message"), split = FALSE)
     
length(which(res$enrichment_result$p.adjustMANOVA<0.05))
hm<-as.matrix(res$enrichment_result[1:40,4:8])
rownames(hm)<-res$enrichment_result[1:40,1]
heatmap.2(hm,trace ="none",scale = "none",dendrogram = "none",margins = c(7,25),cexCol = 0.7)


## ----sessioninfo--------------------------------------------------------------
sessionInfo()

