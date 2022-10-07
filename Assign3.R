
library(readxl)
#library(circlize)
# Reading csv file
dFile<-read.csv("GSE186280.csv",row.names=1,header = T)
print(dFileFile)

# calculating count per matrix of data file
cmat = dFile
for(i in 1:ncol(dFile)){
  cmat[,i]=(dFile[,i]/sum(dFile[,i]))*1000000
}
print(cmat)
cmat[is.na(cmat)]=0

# Calculating log of count per matrix(CPM)
logC = log2(cmat+1)
summary(logC)

#Calculating z-score of logC
z_score = (logC-rowMeans(logC))/rowSds(as.matrix(logC))
row(logC)
z_score[is.na(z_score)]=0
print(z_score)

#Variance
vargenes = apply(z_scores,1,var)
print(vargenes)

vargenes = sort(vargenes,decreasing = T)
top50 = vargenes[1:50]
pmat = z_score[names(top50),]

Hmat = as.matrix(pmat)
#heatmap
#library(ComplexHeatmap)
#heatmap(Hmat)
#creating empty matrix
mat = matrix(NA, ncol= 4, nrow= nrow(dFile))
rownames(mat)= rownames(dFile)
colnames(mat)= c('T','Control','Pvalue','log2FC')
print(mat)

for(i in 1:nrow(dFile)){
  vect1 = as.numeric(dFile[i,1:5])
  vect2 = as.numeric(dFile[i,6:10])
  res = t.test(vect1, vect2, alternative = 'two.sided', paired = F)
  mat[i,1]= res$estimate[[1]]
  mat[i,2]= res$estimate[[2]]
  mat[i,3]= res$p.value
  mat[1,4]=mat[i,1]- mat[i,2]  
  
}

mat= as.data.frame(mat)
num = which(is.nan(mat$Pvalue))
mat[num,'Pvalue']=1
# Differential Expression Gene Analysis(DEGE)
# Plotting Volcano graph
library("EnhancedVolcano")
EnhancedVolcano(mat,lab = rownames(mat),x='log2FC',y= 'Pvalue')



