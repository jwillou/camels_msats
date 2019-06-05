setwd("/Volumes/jwillou/camels/data/")
library(scales)
library(MsatAllele)

#read in data
d1 = read.table("c_saudi.txt", sep="\t", header=T)
d2 = read.table("c_pakistan.txt", sep="\t", header=T)
d1$id = paste(d1$location, d1$id, sep="_")
d2$id = paste(d2$location, d2$id, sep="_")

#list of priemrs and indexing column
d1.pi   = data.frame(primer = colnames(d1)[seq(3, ncol(d1), 2)], index = seq(3, ncol(d1), 2))
d2.pi   = data.frame(primer = colnames(d2)[seq(3, ncol(d1), 2)], index = seq(3, ncol(d1), 2))
primers = unique(c(as.character(d1.pi$primer), as.character(d2.pi$primer)))

#####combine into single data frame####
out = data.frame(id = c(as.character(d1$id), as.character(d2$id)), location = c(d1$location, d2$location))
for(p in primers){
  #d1
  if(is.na(match(p, d1.pi$primer))){
    t1 = cbind(rep(NA, nrow(d1)), rep(NA, nrow(d1)))
  }else(
    t1 = cbind(d1[,d1.pi$index[match(p, d1.pi$primer)]],d1[,(d1.pi$index[match(p, d1.pi$primer)]+1)])
  )
  #d2
  if(is.na(match(p, d2.pi$primer))){
    t2 = cbind(rep(NA, nrow(d2)), rep(NA, nrow(d2)))
  }else(
    t2 = cbind(d2[,d2.pi$index[match(p, d2.pi$primer)]],d2[,(d2.pi$index[match(p, d2.pi$primer)]+1)])
  )
  t = rbind(t1, t2)
  colnames(t) = c(colnames(d1)[d1.pi$index[match(p, d1.pi$primer)]], colnames(d1)[(d1.pi$index[match(p, d1.pi$primer)]+1)])
  if(is.na(colnames(t))){
    colnames(t) = c(colnames(d2)[d2.pi$index[match(p, d2.pi$primer)]], colnames(d2)[(d2.pi$index[match(p, d2.pi$primer)]+1)])
  }
  out = cbind(out, t)
}
data = as.data.frame(out)
#increase YWLL44_F so that alleles are 3 digits later
data$YWLL44_F = data$YWLL44_F + 100
data$YWLL44_R = data$YWLL44_R + 100

####generate graphs to look at####
for(p in 3:length(seq(3, ncol(data), 2))){
  t = c(data[,p], data[,(p+1)])
  t = sort(t)  
  yliml = min(t, na.rm=T)-5
  ylimu = max(t, na.rm=T)+5
  plot(-100, -100, pch=19, col=alpha("red", 0.5), cex=1, main=strsplit(colnames(data)[p], "_")[[1]][1], ylim=c(yliml, ylimu), xlim=c(1, length(t)))
  for(y in seq(yliml, ylimu, 2)){
    abline(h=y, col="grey50", lty=2)
  }
  points(x=seq(1, length(t), 1), y=t, pch=19, col=alpha("red", 0.5), cex=1)
}

####bin alleles####
#set of data object (sdata)
sdata = NULL
for(p in seq(3,35,2)){
  t = cbind(as.character(data$id), rep(p, nrow(data)), rep(strsplit(colnames(data)[p], "_")[[1]][1], nrow(data)), data[,p], data[,(p+1)])
  sdata = rbind(sdata, t)
}
colnames(sdata) = c("sample", "panel", "marker", "a1", "a2")
sdata[is.na(sdata[,5]),4] = NA
sdata[is.na(sdata[,4]),5] = NA
sdata = sdata[complete.cases(sdata),]
write.table(sdata, "sdata.txt", col.names=T, row.names=F, sep="\t", quote=F, append=F)
db = read.frag.sizes(in.file="/Volumes/jwillou/camels/data/sdata.txt","06-04-19",1)

for(p in as.character(unique(db$Marker))){
  AlleleCum(db, p)
  AlleleHist(db,p)
}
write.PG.file.loc(db, outfile="/Volumes/jwillou/camels/data/camel_binned")  







