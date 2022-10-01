

rData = 
vad<-function(rData){
   #library(matrixStats)
  #rData = readline("Enter your data:")
  cpmMat = rData 
  for (i in 1:ncol(rData)) {

    cpmMat[,i] = (rData[,i]/sum(rData[,i]))*1000000
    }
    cpmMat[,i] = log2(cpmMat[,i] +1)
    rlog = log2(cpmMat+1)
    
    fData = cpmMat
    for (i in 1:ncol(data)) {
      zScore = (cpmMat - rowMeans(fData))/rowSds(as.matrix(cpmMat))[row(cpmMat)]
      z_scores[is.na(z_scores)]=0
      zSScore = as.matrix(zScore)
      print(zSScore)
      }
}
vad(rData)
# library(ComplexHeatmap)
pdf("heatMap.pdf",width = 10,height = 10)
heatmap(zSScore)

dev.off()









