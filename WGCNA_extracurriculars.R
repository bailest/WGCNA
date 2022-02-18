#Module and trait specific analyses
#Much closer look for further curation and insight
#Create module heatmaps (up to 3 at a time)
sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="blue";
plotMat(t(scale(trans.exprData[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="turquoise";
plotMat(t(scale(trans.exprData[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="brown";
plotMat(t(scale(trans.exprData[,dynamicColors==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )

#Shows relationship of expression between entire module heatmap and the module eigengene (ME)
#(ME can be considered the most representative gene expression profile of the module)
sizeGrWindow(8,7);
which.module="blue"
ME=MEs[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(trans.exprData[,dynamicColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

#Pairwise scatterplot of eigengenes and a trait
l=datTraits$Larvae
png(filename = "Infected pairwise correlation.png", width = 45,height = 30, units = "in", res = 600)
plotMEpairs(MEs,y=l)
dev.off()

#Determine module significance to a particular chosen trait
GS1=as.numeric(cor(Infected,trans.exprData, use="p"))
GeneSignificance=abs(GS1)
#Determine genes of interest in a trait based on high significance and intramodular connectivity
datKME=signedKME(trans.exprData, MEs, outputColumnName="MM.")
FilterGenes= abs(GS1)> .2 & abs(datKME$MM.brown)>.8
sigGenes=data.frame(dimnames(data.frame(trans.exprData))[[2]][FilterGenes])
write.csv(sigGenes, file = "Infected significant genes.csv")
# Next module significance is defined as average gene significance.
png(filename = "Infected significant modules.png", width = 45,height = 30, units = "in", res = 600)
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,dynamicColors)
dev.off()

#Dendrogram depicting hierarchical clustering of samples to traits
sampleTree2 = hclust(dist(trans.exprData), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

#Plotting gene significance vs intramodular connectivity
#Genes highly significant within a trait and highly connected to the module are important
ADJ1=abs(cor(trans.exprData,use="p"))^6
Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors)
colorlevels=unique(dynamicColors)
png(filename = "Larvae sig v connectivity.png", width = 50,height = 10, units = "in", res = 600)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (dynamicColors==whichmodule);
  verboseScatterplot(Alldegrees1$kWithin[restrict1],
                     GeneSignificance[restrict1], col=dynamicColors[restrict1],
                     main=whichmodule,
                     xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
}
dev.off()

#Module membership measure^6 vs intramodular connectivity
which.color="blue";
restrictGenes=dynamicColors==which.color
verboseScatterplot(Alldegrees1$kWithin[ restrictGenes],
                   (datKME[restrictGenes, paste("MM.", which.color, sep="")])^6,
                   col=which.color,
                   xlab="Intramodular Connectivity",
                   ylab="(Module Membership)^6")

