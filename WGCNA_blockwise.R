setwd("~/R")
#Packages to install
#install.packages(c("dynamicTreeCut", "cluster", "flashClust", "Hmisc", "reshape", "foreach", "doParallel") ) 
#source("http://bioconductor.org/biocLite.R") 
#biocLite("impute")
#install.packages("WGCNA")
####
#BiocManager::install(c("GO.db","preprocessCore","impute"))
library(WGCNA)
library(flashClust)

options(stringsAsFactors = FALSE)
#enableWGCNAThreads()
#allowWGCNAThreads()
#^Allows WGCNA to multithread
exprData <- read.delim(".txt", row.names=1, header = TRUE)

#gene.names=rownames(exprData)
#^Input data should be nothing but expression values labelled by sample and gene id
trans.exprData=t(exprData)
#^WGCNA expects gene ids to be in columns, this line transposes the data

#Removing genes with 0 variance
gsg = goodSamplesGenes(trans.exprData, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(trans.exprData)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(trans.exprData)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  trans.exprData = trans.exprData[gsg$goodSamples, gsg$goodGenes]
}
gene.names=colnames(trans.exprData)

#Loading in trait data
datTraits = read.csv(".csv")
#form a data frame analogous to expression data that holds the clinical traits
rownames(datTraits) = datTraits$Sample
datTraits$Sample = NULL
table(rownames(datTraits)==rownames(trans.exprData)) 
#^should return TRUE if datasets align correctly, otherwise your names are out of order
#The row labels in trait file have to perfectly match with the column labels of expression data

#Change the "networkType" argument from signed to unsigned or hybrid if you want to; signed is preferred
#Dont forget to change signed everywhere else it appears from here on out
powers = c(c(1:10), seq(from = 12, to=20, by=2));
sft=pickSoftThreshold(trans.exprData,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "signed")
#Plot results
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
pdf("Connectivity.pdf",width=5,height=5)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
#Red line corresponds to using an R^2 cut-off
abline(h=0.80,col="red")
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#Generating adjacency and TOM similarity matrices based on the selected softpower
#Choose softpower based on the lowest number that reaches the red line from the Scale independence plot just generated
#if it does not reach that read line, default to 18
softPower = 18;

#define maxBlockSize based on how many transcripts your data set has
bwnet = blockwiseModules(trans.exprData, maxBlockSize = 40000,
                         power = softPower, TOMType = "signed", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = .75,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "TOM-blockwise",
                         verbose = 3)

moduleLabels = bwnet$colors 
moduleColors = labels2colors(bwnet$colors) 
MEs = bwnet$MEs 
geneTree = bwnet$dendrograms[[1]] 

#for creating this plot, plotDendroAndColors() needs to be called individually for each
#block that was created. if there are 40,000 genes per block and we have 160,000 genes
#then we have 4 blocks and this must be called 4 times as shown below. Adjust accordingly
pdf("Genetree.pdf",width=15,height=5)
plotDendroAndColors(bwnet$dendrograms[[1]], moduleColors[bwnet$blockGenes[[1]]],
                    "Module colors", 
                    dendroLabels = F, hang = 0.03,
                    
                    addGuide = TRUE, guideHang = 0.05)
plotDendroAndColors(bwnet$dendrograms[[2]], moduleColors[bwnet$blockGenes[[2]]],
                    "Module colors", #
                    dendroLabels = F, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(bwnet$dendrograms[[3]], moduleColors[bwnet$blockGenes[[3]]],
                    "Module colors", #
                    dendroLabels = F, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

plotDendroAndColors(bwnet$dendrograms[[4]], moduleColors[bwnet$blockGenes[[4]]],
                    "Module colors", #
                    dendroLabels = F, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Correlate traits --------------------------------------------------------
#Define number of genes and samples
nGenes = ncol(trans.exprData)
nSamples = nrow(trans.exprData)
#Recalculate MEs with color labels
MEs0 = moduleEigengenes(trans.exprData, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#Print correlation heatmap between modules and traits
textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
par(mar= c(6, 8.5, 3, 3))

#display the corelation values with a heatmap plot

pdf("Module relationships.pdf",width=25,height=15)
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs),
               ySymbols= names(MEs),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= TRUE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
dev.off()
#For a closer look, if you wanted to pull out genes belonging to a certain module, you can use the following command:
#names(trans.exprData)[dynamicColors=="brown”]

#Creates text files for each color module and a folder to place them in
dir.create("Modules")
module_colors= setdiff(unique(dynamicColors), "")
for (color in module_colors){
  module=gene.names[which(dynamicColors==color)]
  write.table(module, paste("~/Modules/module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}
#Create spreadsheets that contain the significance and correlation values displayed on the "Module Relationships" heatmap
write.table(moduleTraitPvalue, "moduleTraitPvalue.txt", sep="\t")
write.table(moduleTraitCor, "moduleTraitCor.txt", sep="\t")

#Heatmap showing expression of all genes groups by module designation
pdf("heatmap.pdf",width=25,height=15)
module.order <- unlist(tapply(1:ncol(trans.exprData),as.factor(dynamicColors),I))
m<-t(t(trans.exprData[,module.order])/apply(trans.exprData[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",RowSideColors=dynamicColors[module.order])
dev.off()

#Visualize the Tom plot. Raise the dissimilarity matrix to the power of 4 to bring out the module structure
png(filename = "TOM.png", width = 6,height = 6, units = "in", res = 600)
dissTOM = 1-TOM;
plotTOM = dissTOM^7;
diag(plotTOM) = NA;
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
dev.off()

#Multi-dimensional scaling plot
pdf("mds plot.pdf",width=12,height=12)
cmd1=cmdscale(as.dist(dissTOM),2)
par(mfrow=c(1,1))
plot(cmd1, col=as.character(dynamicColors), main="MDS plot",
     xlab="Scaling Dimension 1", ylab="Scaling Dimension 2")
dev.off()

#Plot trait-module relationships
# Recalculate module eigengenes
MEs = moduleEigengenes(trans.exprData, dynamicColors)$eigengenes
# Isolate trait data for plotting
#change the following to match the traits defined in your trait file
Control = as.data.frame(datTraits$Control)
names(Control) = "Control"
Chloropyrifos = as.data.frame(datTraits$Chloropyrifos)
names(Chloropyrifos) = "Chloropyrifos"
Propoxur = as.data.frame(datTraits$Propoxur)
names(Propoxur) = "Propoxur"
Fipronil = as.data.frame(datTraits$Fipronil)
names(Fipronil) = "Fipronil"
Amitraz = as.data.frame(datTraits$Amitraz)
names(Amitraz) = "Amitraz"
Heatshock = as.data.frame(datTraits$Heatshock)
names(Heatshock) = "Heatshock"
Rapid.Cold.Hardening = as.data.frame(datTraits$Rapid.Cold.Hardening)
names(Rapid.Cold.Hardening) = "Rapid Cold Hardening"
Dehydration = as.data.frame(datTraits$Dehydration)
names(Dehydration) = "Dehydration"
Starvation = as.data.frame(datTraits$Starvation)
names(Starvation) = "Starvation"


# Add the traits to existing module eigengenes
MET = orderMEs(cbind(MEs, Control,	Chloropyrifos,	Propoxur,	Fipronil,	Amitraz,	Heatshock,	Rapid.Cold.Hardening,	Dehydration,	Starvation))
# Plot the relationships among the eigengenes and the trait

pdf("Eigengenes heat and dend.pdf",width=30,height=30)
par(cex = 1)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)
dev.off()
#Plot the dendgrogram
pdf("Eigengenes dend.pdf",width=20,height=20)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0), plotHeatmaps = FALSE)
dev.off()
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf("Eigengenes heat.pdf",width=20,height=20)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
dev.off()

#####################################################################################################
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

