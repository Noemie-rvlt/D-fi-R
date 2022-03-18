setwd("~/R/défi R")

# Importation des données
expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1, sep = "\t", header = T)

#Choisir ses paramètres : choisir N, method (kmeans ou HCL) et distance (euclidean ou correlation)
N = 10
method = "HCL"
distance="euclidean"

# Ma fonction :)
plotGenes <- function(expData, title = "", yMin = 0, yMax = NULL, meanProfile = TRUE){
  
  # Check function parameters
  if(is.null(yMax)){
    
    print("You must specify a maximal value for Y axis")
    
  }else{
    
    # Representation of the first expression profile
    plot(1:ncol(expData), expData[1,], col = "grey", type = "l",
         ylim = c(floor(yMin), ceiling(yMax)),
         xlab = "Time point", ylab = "Gene expression level",
         main = title)
    
    # Add expression profile for other genes
    for(i in 2:nrow(expData)){
      
      lines(1:ncol(expData), expData[i,], col = "grey")
      
      # end of for()  
    }
    
    # Average expression profile
    if(meanProfile == TRUE){
      expMean = apply(expData, 2, mean)
      lines(1:ncol(expData), expMean, col = "red", 
            lwd = 1.5, lty = "dashed")
    }
    
    # end of else()   
  }
  
  # end of function plotGenes()  
}

# Clustering 
if((method=="kmeans")&(distance=="euclidean")){
plotGenes(expData, yMax = 100000)
for(i in 1:N){
  cluster <- expData[which(resKmeans$cluster == i),]
  plotGenes(cluster, yMax = 100000)
  if(nrow(cluster)>1){
    heatmap(as.matrix(cluster))}
}}else if((method=="HCL")&(distance=="euclidean")){
  d<- dist(expData)
  resHCL <- hclust(d)
  plot(resHCL)
  for(i in 1:N){
    cluster <- expData[which(cutree(resHCL, k = N) == i),]
    plotGenes(cluster, yMax = 100000)
    if(nrow(cluster)>1){
      heatmap(as.matrix(cluster))}
  }
}else if((method=="kmeans")&(distance=="correlation")){
  cor(expData)
  cor(t(expData))
  matDist <- as.dist(1 - cor(t(expData)))
  resKmeans2 <- kmeans(matDist, centers = N)
  plotGenes(expData, yMax = 100000)
  for(i in 1:N){
    cluster <- expData[which(resKmeans2$cluster == i),]
    plotGenes(cluster, yMax = 100000)
    if(nrow(cluster)>1){
      heatmap(as.matrix(cluster))}
  }
}else if((method=="HCL")&(distance=="correlation")){
  cor(expData)
  cor(t(expData))
  matDist <- as.dist(1 - cor(t(expData)))
  resHCL2 <- hclust(matDist)
  plot(resHCL2)
for(i in 1:N){
  cluster <- expData[which(cutree(resHCL2, k = N) == i),]
  if(nrow(cluster)>1){
    heatmap(as.matrix(cluster))}
}}
