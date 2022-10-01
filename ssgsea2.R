setwd("Documents")
getwd()
library("ComplexHeatmap")
library("readxl")
library("circlize")
dAta<-read.csv("C:\\Users\\RUCHIRA\\Documents\\CG4.csv",row.names = 1)
print(dAta)
summary(dAta)

## counts per matrix(CPM)
c5pm<- dAta
for (i in 1:ncol(dAta)) {
  c5pm[,i]=(dAta[,i]/sum(dAta[,i]))*1000000  
}
print(c5pm)

##log of cpm
cLog<-log2(c5pm+1)
print(cLog)

## z value of log
library("matrixStats")

zVal<-(cLog - rowMeans(cLog)) / rowSds(as.matrix(cLog))[row(cLog)]
zVal[is.na(zVal)]=0
print(zVal)


## variance
varS<-apply(zVal, 1,varDiff)
print(varS)
varS = sort(varS,decreasing = T)
top50 = varS[1:50]
print(top50)
zVal2<-zVal[names(top50),]

caData<-matrix(rnorm(100,0,5),nrow = 10,ncol = 10)
print(caData)

##heatmap
library(ComplexHeatmap)
HM<-heatmap(caData, xlab = "Species",ylab = "Range")
print(HM)

colnames(caData)<-paste0("gene",1:10)
rownames(caData)<-paste0("species",1:10)


##heatmap(caData,xlab = "Sample",ylab="gene")
Zs <- (caData - mean(caData, na.rm = TRUE)) / sd(caData, na.rm = TRUE) 
Zs

#counts per million matrix
#CPMM = apply(dAta, 2, function(x) (x/sum(x))*1000000)
#print(CPMM)

# creating empty matrix
bmat = matrix(NA, ncol= 4, nrow= nrow(caData))
rownames(bmat)= rownames(caData)
colnames(bmat)= c('test','control','pval','log2FC')
print(bmat)

#l2cpm = log2(mat +1)
#print(l2cpm)

for(i in 1:nrow(caData)){
  vector1 = as.numeric(caData[i,1:4])
  vector2 = as.numeric(caData[i,5:10])
  
  res = t.test(vector1, vector2, paired = F, alternative = 'two.sided')
  bmat[i,1]= res$estimate[[1]]
  bmat[i,2]= res$estimate[[2]]
  bmat[i,3]= res$p.value
  bmat[i,4]=mat[i,1]- mat[i,2]
  
}
#storing in datafrmae
bmat= as.data.frame(bmat)
num = which(is.nan(bmat$pval))
mat[num,'pval']=1

# for volcano plot
# install EnhancedVolcano
library(EnhancedVolcano)
EnhancedVolcano(bmat,lab = rownames(bmat),x='log2FC',y= 'pval')  # error in data

