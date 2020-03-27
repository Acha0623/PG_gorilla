## LD Dacey
library("snpMatrix")

## Gbb
gbb<- read.plink("Gbb")
ld<- ld.snp(gbb, dep=500)
do.rsq2 <- ("rsq2" %in% names(ld))
snp.names <- attr(ld, "snp.names")
name<- c("#M1", "#M2", "rsq2", "Dprime", "lod")
r.maybe <- ld$rsq2
max.depth <- dim(ld$dprime)[2]
res<-matrix(NA,ncol=3,nrow=length(snp.names)*max.depth)
count<-1
for (i.snp in c(1:(length(snp.names) - 1))) {
  for (j.snp in c((i.snp + 1):length(snp.names))) {
    step <- j.snp - i.snp
    if (step > max.depth) {
      break
    }
    res[count,]<-c(snp.names[i.snp], snp.names[j.snp],r.maybe[i.snp, step])
    count<-count+1;
  }
}
resNum<-as.numeric(res)#converts to numeric values
dim(resNum)<-dim(res)
dis<- resNum[,2]-resNum[,1] #calculates distances between sites
disbin<- cut(dis, br=seq(0,2000000,by=10000))
#cuts distances into bins of size 1K from distance "0" to "2000000"
gbb<- tapply(resNum[,3], disbin, mean, na.rm=T)
ressum<- tapply(resNum[,3], disbin, sum, na.rm=T)
func<- function(x)#makes a function "x"
  length(x[!is.na(x)])
#gives the lengths of what is not (denoted by "!") "NA" in "x"
reslength<- tapply(resNum[,3], disbin, func)
#gives the length in each bin of 1K for each value of r2
write.table(reslength, file="Gbblenght.txt")
write.table(ressum, file="Gbbsum.txt")

## Gbg
gbg<- read.plink("Gbg")
ld<- ld.snp(gbg, dep=500)
do.rsq2 <- ("rsq2" %in% names(ld))
snp.names <- attr(ld, "snp.names")
name<- c("#M1", "#M2", "rsq2", "Dprime", "lod")
r.maybe <- ld$rsq2
max.depth <- dim(ld$dprime)[2]
res<-matrix(NA,ncol=3,nrow=length(snp.names)*max.depth)
count<-1
for (i.snp in c(1:(length(snp.names) - 1))) {
  for (j.snp in c((i.snp + 1):length(snp.names))) {
    step <- j.snp - i.snp
    if (step > max.depth) {
      break
    }
    res[count,]<-c(snp.names[i.snp], snp.names[j.snp],r.maybe[i.snp, step])
    count<-count+1;
  }
}
resNum<-as.numeric(res)#converts to numeric values
dim(resNum)<-dim(res)
dis<- resNum[,2]-resNum[,1] #calculates distances between sites
disbin<- cut(dis, br=seq(0,2000000,by=10000))
#cuts distances into bins of size 1K from distance "0" to "2000000"
gbg<- tapply(resNum[,3], disbin, mean, na.rm=T)
ressum<- tapply(resNum[,3], disbin, sum, na.rm=T)
func<- function(x)#makes a function "x"
  length(x[!is.na(x)])
#gives the lengths of what is not (denoted by "!") "NA" in "x"
reslength<- tapply(resNum[,3], disbin, func)
#gives the length in each bin of 1K for each value of r2
write.table(reslength, file="Gbglength.txt")
write.table(ressum, file="Gbgsum.txt")


## Ggg
ggg<- read.plink("Ggg")
ld<- ld.snp(ggg, dep=500)
do.rsq2 <- ("rsq2" %in% names(ld))
snp.names <- attr(ld, "snp.names")
name<- c("#M1", "#M2", "rsq2", "Dprime", "lod")
r.maybe <- ld$rsq2
max.depth <- dim(ld$dprime)[2]
res<-matrix(NA,ncol=3,nrow=length(snp.names)*max.depth)
count<-1
for (i.snp in c(1:(length(snp.names) - 1))) {
  for (j.snp in c((i.snp + 1):length(snp.names))) {
    step <- j.snp - i.snp
    if (step > max.depth) {
      break
    }
    res[count,]<-c(snp.names[i.snp], snp.names[j.snp],r.maybe[i.snp, step])
    count<-count+1;
  }
}
resNum<-as.numeric(res)#converts to numeric values
dim(resNum)<-dim(res)
dis<- resNum[,2]-resNum[,1] #calculates distances between sites
disbin<- cut(dis, br=seq(0,2000000,by=10000))
#cuts distances into bins of size 1K from distance "0" to "2000000"
ggg<- tapply(resNum[,3], disbin, mean, na.rm=T)
ressum<- tapply(resNum[,3], disbin, sum, na.rm=T)
func<- function(x)#makes a function "x"
  length(x[!is.na(x)])
#gives the lengths of what is not (denoted by "!") "NA" in "x"
reslength<- tapply(resNum[,3], disbin, func)
#gives the length in each bin of 1K for each value of r2
write.table(reslength, file="Ggglength.txt")
write.table(ressum, file="Gggsum.txt")

#saves output as a .pdf
pdf("LDdecay_chr21.pdf") 
plot(gbb, type="l", col='darkred', ylim = c(0,1), ylab = 'R2', xlab = 'distance (10kbp)', lwd = 2, cex.lab=1.5,cex.axis=1.5)
lines(gbg, col='chartreuse3', lwd = 1.5)
lines(ggg, col='dodgerblue2', lwd = 1.5)
legend("topright",legend=c("Gbb","Gbg","Ggg"),
       col=c("darkred","chartreuse3","dodgerblue2"), lty = 1, lwd = 1.5,cex=1.5)
print('Done LD decay')
