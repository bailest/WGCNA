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

#calclute the adjacency matrix
adj= adjacency(trans.exprData,type = "signed", power = softPower);
#turn adjacency matrix into topological overlap to minimize the effects of noise and spurious associations
TOM=TOMsimilarityFromExpr(trans.exprData,networkType = "signed", TOMType = "signed", power = softPower);
colnames(TOM) =rownames(TOM)
dissTOM=1-TOM
geneTree = flashClust(as.dist(dissTOM),method="average");
#plot the resulting clustering tree (dendrogram)
pdf("Genetree.pdf",width=15,height=5)
plot(geneTree, xlab="", sub="",cex=0.3);
minModuleSize = 20;
#Module identification using dynamic tree cut
#if this first dynamicmods line fails, comment it out and use the other one
dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
#dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, method="hybrid", deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize);
#the following command gives the module labels and the size of each module. Label 0 is reserved for unassigned genes
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(trans.exprData[,restGenes], power = softPower)
colnames(diss1) =rownames(diss1)
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
diag(diss1) = NA;
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
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
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
