setwd("~/R/défi R")

# Importation des données
expData <- read.table("cell-cycle_SCERE_DUO.txt", row.names = 1, sep = "\t", header = T)

# Comment sont les données ?
sum(is.na(expData))
# --> Il n'y pas de NA
boxplot(expData)

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

# Clustering avec la méthode k-means
N = 4
resKmeans <- kmeans(expData, centers = N)
plotGenes(expData, yMax = 100000)
cluster1 <- expData[which(resKmeans$cluster == 1),]
plotGenes(cluster1, yMax = 100000)
if(nrow(cluster1)>1){heatmap(as.matrix(cluster1))}
cluster2 <- expData[which(resKmeans$cluster == 2),]
plotGenes(cluster2, yMax = 100000)
cluster3 <- expData[which(resKmeans$cluster == 3),]
plotGenes(cluster3, yMax = 100000)
cluster4 <- expData[which(resKmeans$cluster == 4),]


# Classification HCL
d<- dist(expData)
resHCL <- hclust(d)
plot(resHCL)
cluster1 <- expData[which(cutree(resHCL, k = N) == 1),]
plotGenes(cluster1, yMax = 100000)
cluster2 <- expData[which(cutree(resHCL, k = N) == 2),]
plotGenes(cluster2, yMax = 100000)
cluster3 <- expData[which(cutree(resHCL, k = N) == 3),]
plotGenes(cluster3, yMax = 100000)
cluster4 <- expData[which(cutree(resHCL, k = N) == 4),]
plotGenes(cluster4, yMax = 1000000)

# Création d'une matrice de distance
cor(expData)
cor(t(expData))
matDist <- as.dist(1 - cor(t(expData)))

# Kmeans et correlation
resKmeans2 <- kmeans(matDist, centers = 4)
plotGenes(expData, yMax = 100000)
cluster1 <- expData[which(resKmeans2$cluster == 1),]
heatmap(as.matrix(cluster1))
plotGenes(cluster1, yMax = 100000)
cluster2 <- expData[which(resKmeans2$cluster == 2),]
plotGenes(cluster2, yMax = 100000)
heatmap(as.matrix(cluster2))
cluster3 <- expData[which(resKmeans2$cluster == 3),]
plotGenes(cluster3, yMax = 100000)
heatmap(as.matrix(cluster3))
cluster4 <- expData[which(resKmeans2$cluster == 4),]
plotGenes(cluster4, yMax = 100000)
heatmap(as.matrix(cluster4))

# HCL et correlation
resHCL2 <- hclust(matDist)
plot(resHCL2)
cluster1 <- expData[which(cutree(resHCL2, k = N) == 1),]
plotGenes(cluster1, yMax = 100000)
heatmap(as.matrix(cluster1))
cluster2 <- expData[which(cutree(resHCL2, k = N) == 2),]
plotGenes(cluster2, yMax = 100000)
heatmap(as.matrix(cluster2))
cluster3 <- expData[which(cutree(resHCL2, k = N) == 3),]
plotGenes(cluster3, yMax = 100000)
heatmap(as.matrix(cluster3))
cluster4 <- expData[which(cutree(resHCL2, k = N) == 4),]
plotGenes(cluster4, yMax = 1000000)
heatmap(as.matrix(cluster4))










