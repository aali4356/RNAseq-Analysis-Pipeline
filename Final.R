library(dplyr)
library(ggplot2)
library(tidyr)
library(plotly)
library(infer) 
library(ggfortify)
library(broom)
library(tibble)
library(gridExtra)
library(rsample)
library(Biobase)
library(reshape2)
library(plyr)
library(e1071)
library(ShortRead)
library(rlang)
library(tidyverse)
library(rnaseqGene)

# For the second set of data
library(edgeR)
library(DESeq2)
library(limma)
# library(Glimma)
# library(gplots)
library(org.Mm.eg.db)
# library(RColorBrewer)
library(GEOquery)
library(tibble)
library(tidyverse)
library(gplots)
library(Homo.sapiens)
library(topGO)


# Load Data ---------------------------------------------------------------


gset <- getGEO("GSE145652", GSEMatrix =TRUE, getGPL=FALSE)

sampleinfo = pData(gset[[1]])[, 1:3]
sampleinfo$status = factor(c(rep('CON',3), rep('IGA', 3)))
sampleinfo = as.data.frame(sampleinfo)


# load counts table from GEO
urld <- "https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts"
path <- paste(urld, "acc=GSE145652", "file=GSE145652_raw_counts_GRCh38.p13_NCBI.tsv.gz", sep="&");
tbl <- as.matrix(data.table::fread(path, header=T, colClasses="integer"), rownames=1)

# pre-filter low count genes
# keep genes with at least 2 counts > 10
keep <- rowSums( tbl >= 10 ) >= 2
tbl <- tbl[keep, ]
countsData <- tbl
# colnames(countsData) <- c('CON1', 'CON2', 'CON3', 'IGA1', 'IGA2', 'IGA3')
# log transform raw counts
# instead of raw counts can display vst(as.matrix(tbl)) i.e. variance stabilized counts
log.normalized.counts <- log10(tbl + 1)
colnames(log.normalized.counts) <- c('CON1', 'CON2', 'CON3', 'IGA1', 'IGA2', 'IGA3') 



# Differential Expression Analysis ----------------------------------------


# Cleaning data for DESeq
meta.data <- data.frame(
  id = colnames(log.normalized.counts),
  condition = c(rep("control", times = 3), rep("IgA", times = 3))
)

# Establishing variables specifically for this dataset
conditions <- factor(c(rep("control", times = 3), rep("IgA", times = 3)))
design <- model.matrix(~0 + conditions)
colnames(design) <- levels(conditions)

# columns are not exact
col.data <- sampleinfo
col.data <- col.data[,c(1,3)]
# making sure the row names in col.data matches to column names in countsData
colnames(countsData) <- rownames(col.data)

# are they in the same order?
all(colnames(countsData) == rownames(col.data))


dds <- DESeqDataSetFromMatrix(countData = countsData,
                              colData = col.data,
                              design = ~ status)

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Make sure that we have high reads
dds

# comparing other levels to untreated
dds$status <- relevel(dds$status, ref = "CON")

dds <- DESeq(dds)
res <- results(dds)

# DESeq Results
head(res)

# For Cytoscape
DESeq.analysis.table <- data.frame(res)

library(org.Hs.eg.db)

geneID <- rownames(res)
geneEnsembl <- mapIds(org.Hs.eg.db, keys = geneID, column = "ENSEMBL", keytype = "ENTREZID")
geneID <- rownames(res)
geneSymbol <- mapIds(org.Hs.eg.db, keys = geneID, column = "SYMBOL", keytype = "ENTREZID")

DESeq.analysis.table <- data.frame(ENTREZ = rownames(res),
                                   symbol_name = unname(geneSymbol), 
                                   ensembl_id = unname(geneEnsembl),
                                   DESeq.analysis.table)
# Sort by 'adj.P.Val' column in ascending order (lowest values first)
top.DESeq <- DESeq.analysis.table[order(DESeq.analysis.table$padj), ]
head(DESeq.analysis.table[order(DESeq.analysis.table$padj), ] , 15)
write.csv(DESeq.analysis.table, "DESeq.analysis.table.csv", row.names=FALSE)


# Limma differential expression analysis
# linear model fit off of Bayes moderation of standard errors
fit <- lmFit(countsData, design)
contrast.matrix <- makeContrasts(IgA - control, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Get the results and select differentially expressed miRNAs
results <- topTable(fit2, adjust="BH", number=Inf)


# Displaying results from limma analysis
geneID <- rownames(results)
geneEnsembl <- mapIds(org.Hs.eg.db, keys = geneID, column = "ENSEMBL", keytype = "ENTREZID")
geneID <- rownames(results)
geneSymbol <- mapIds(org.Hs.eg.db, keys = geneID, column = "SYMBOL", keytype = "ENTREZID")

limma.analysis.table2 <- data.frame(ENTREZ = rownames(results),
                                   symbol_name = unname(geneSymbol), 
                                   ensembl_id = unname(geneEnsembl),
                                   results)
head(limma.analysis.table2)
dim(limma.analysis.table2)
write.csv(limma.analysis.table2, "analysis.table.csv", row.names=FALSE)


# Volcano Plot of DEG -----------------------------------------------------

get_upregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange>=1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

get_downregulated <- function(df){
  
  key <- intersect(rownames(df)[which(df$log2FoldChange<=-1)], rownames(df)[which(df$pvalue<=0.05)])
  
  results <- as.data.frame((df)[which(rownames(df) %in% key),])
  return(results)
}

library(EnhancedVolcano)
res1 <- na.omit(res)

min_width <- min(res1$log2FoldChange)
max_width <- max(res1$log2FoldChange)
max_height <- -log10(min(res1[res1$pvalue>0, 5]))

up <- subset(res1, res1$log2FoldChange > 1 & res1$pvalue <= 0.05)
up <- up[order(-up$log2FoldChange),]
up_list <- head(rownames(up), n=10L)

down <- subset(res1, res1$log2FoldChange < -1 & res1$pvalue <= 0.05)
down <- down[order(down$log2FoldChange),]
down_list <- head(rownames(down), n=10L)

plot_top_20 <- c(up_list, down_list)
EnhancedVolcano(res1,
                lab=rownames(res1),
                x="log2FoldChange",
                y="pvalue",
                selectLab=plot_top_20,
                drawConnectors=TRUE,
                FCcutoff=1.0,
                pCutoff=0.05,
                title="Volcano Plot",
                subtitle="IgAN vs. Control",
                # legendVisible=F,
                caption = paste0('Total Genes = ', nrow(res1)),
                xlim=c(min_width, max_width),
                ylim=c(0, max_height))


# Potential Biomarkers/Transcription Factors ------------------------------

# Fetch the data from the Human TFs database
tf_data <- read_csv("http://humantfs.ccbr.utoronto.ca/download/v_1.01/DatabaseExtract_v_1.01.csv")

# Extract the list of TFs
known_TFs <- tf_data$`Ensembl ID`

# Print the first few TFs
head(known_TFs)

# Sorting the differentially expressed genes by adj.P.Val and logFC
sorted_analysis_table <- DESeq.analysis.table[order(DESeq.analysis.table$padj, -abs(DESeq.analysis.table$log2FoldChange)),]

# Select the top N genes as potential biomarkers
topN <- 300
potential_biomarkers <- sorted_analysis_table[1:topN,]

# Intersect the differentially expressed genes with known transcription factors
TF_intersect <- intersect(na.omit(sorted_analysis_table$ensembl_id), known_TFs)

# Select these transcription factors from the differentially expressed genes
potential_TF_biomarkers <- sorted_analysis_table[sorted_analysis_table$ensembl_id %in% TF_intersect,]

head(potential_TF_biomarkers,15)
write.csv(potential_TF_biomarkers, "potential_TF_biomarkers.csv", row.names=FALSE)


# T-test and SAM ----------------------------------------------------------

# Pulling Likelihood Ratio Test P-values
res.pvalue.df <- data.frame(
  ID = top.DESeq$symbol_name,
  p.value = top.DESeq$pvalue,
  p.adj = top.DESeq$padj
)
head(res.pvalue.df)


get.t.test <- function(g.data, g.row){
  
  t.result <- try(t.test(g.data[g.row , 1:3], g.data[g.row , 4:6])$p.value)
  
  if(grepl(pattern = "Error", x = t.result)){ 
    return(NA)
  } else {
    t.result
  }
  
}

t.test.results <- matrix(0, ncol = 2, nrow = nrow(countsData))
for (gene in 1:nrow(countsData)) {
  # print(gene)
  t.test.results[gene,1] <- rownames(countsData)[gene]
  t.test.results[gene,2] <- get.t.test(countsData, gene)
}

t.test.results <- na.omit(t.test.results)
t.test.ordered <- t.test.results[order(t.test.results[,2]),]

t.test.df <- data.frame(
  ID = t.test.ordered[, 1],
  p.value = as.numeric(t.test.ordered[, 2]),
  p.adj = p.adjust(as.numeric(t.test.ordered[,2]), method = "BH")
)
head(t.test.df)

# SAM Analysis
library(samr)
gene.data <- as.matrix(countsData)

l <- c(1,1,1,2,2,2)
g <- factor(l)    

sg = sub("0","2",l)                 # grouping
sm= list(x = gene.data, y= sg, logged2=TRUE)           # SAM input matrix
st = samr(sm,resp.type="Two class unpaired",nperm=100)  # sam test

samr.pvalues <- samr.pvalues.from.perms(st$tt, st$ttstar)

samr.pvalues.df <- as.data.frame(samr.pvalues)
tmp.pval <- samr.pvalues.df[order(samr.pvalues.df$samr.pvalues),]
samr.pvalues.df$samr.padj <- p.adjust(tmp.pval, method = "BH")

samr.pvalues.df[,1] <- tmp.pval

# Still no significant p-values here
head(samr.pvalues.df)

# Histograms
# P-value Plot ------------------------------------------------------------

# T-Test
ttest.pvalue.plot <- 
  ggplot(data = t.test.df, aes(x = p.value)) + 
  geom_histogram() + 
  ggtitle("Non-adjusted T-test P-values")

ttest.padj.plot <- 
  ggplot(data = t.test.df, aes(x = p.adj)) + 
  geom_histogram() + 
  ggtitle("Adjusted T-test P-values")

# DeSeq2 LRT
res.pvalue.plot <- 
  ggplot(data = res.pvalue.df, aes(x = p.value)) + 
  geom_histogram() + 
  ggtitle("Non-adjusted LRT P-values")

res.padj.plot <- 
  ggplot(data = res.pvalue.df, aes(x = p.adj)) + 
  geom_histogram() + 
  ggtitle("Adjusted LRT P-values")

# SAM Analysis Test
samr.pvalues.plot <- 
  ggplot(data = samr.pvalues.df, aes(x = samr.pvalues)) + 
  geom_histogram() + 
  ggtitle("Non-adjusted SAM P-values")

samr.padj.plot <- 
  ggplot(data = samr.pvalues.df, aes(x = samr.padj)) + 
  geom_histogram() + 
  ggtitle("Adjusted SAM P-values")


# Plotting Data
grid.arrange(ttest.pvalue.plot,
             ttest.padj.plot,
             res.pvalue.plot,
             res.padj.plot,
             samr.pvalues.plot,
             samr.padj.plot,
             ncol = 2)

data.frame(
  T.Test.pvalue = sum(t.test.df$p.value < 0.05, na.rm=TRUE),
  T.Test.pvalue = sum(t.test.df$p.adj < 0.05, na.rm=TRUE),
  DeSeq2.pvalue = sum(res.pvalue.df$p.value < 0.05, na.rm=TRUE),
  DeSeq2.padj = sum(res.pvalue.df$p.adj < 0.05, na.rm=TRUE),
  SAM.pvalue = sum(samr.pvalues.df$samr.pvalues < 0.05, na.rm=TRUE),
  SAM.padj = sum(samr.pvalues.df$samr.padj < 0.05, na.rm=TRUE)
)


# Hierarchical clusters ---------------------------------------------------

# data(airway)
# dds=DESeq2::DESeqDataSet(airway, design = ~cell+dex)
dds <- dds[ rowSums(counts(dds)) > 1, ] 
dds <- estimateSizeFactors(dds)
countdata <- countsData
countdata <- countdata[rowSums(countdata) > 0, ]
countdata.deseq <- sweep(countdata, 2 , sizeFactors(dds), FUN = "/")
rld=rlog(dds)
vars <- apply(assay(rld), 1, var)
# labels <- paste(colData(airway)[,"cell"],
#                 ifelse(colData(airway)[,"dex"]=="trt","T","U"),
#                 sep=":")
corr.dist <- function(x) {
  as.dist(1-cor(t(x))) 
}

select <- order(vars, decreasing = T)[1:300]
temp <- assay(rld[select,])
# colnames(temp) <- labels
temp <- data.frame(temp)
# temp$gene <- row.names(temp)
data <- temp
data <- na.omit(data) # listwise deletion of missing
data <- scale(data) # standardize variables 

# Ward Hierarchical Clustering
d <- dist(data, method = "euclidean") # distance matrix
dev.off()
hier.clust <- hclust(d, method="ward.D2")
# par(mar = rep(2, 4))
plot(hier.clust, lab = F) # display dendogram
groups <- cutree(hier.clust, k=4) # cut tree into 4 clusters
# draw dendogram with red borders around the 4 clusters
rect.hclust(hier.clust, k=4, border="red") 

d.sample <- dist(t(data), method = "euclidean")
fit <- hclust(d.sample, method = "ward.D2")
plot(fit, lab = F)
groups <- cutree(fit, k = 3)
rect.hclust(fit, k = 3, border = "red")


# K-Means -----------------------------------------------------------------


library (NbClust)
library (cluster)
# library (clustertend) #depreciated
library(hopkins)
library (factoextra)

library (fpc)

# Overall view for k-means clustering approach
plot(data)

hopkins.test <- 
  hopkins(data, m = nrow(data) - 1)

fviz_dist(dist(data), show_labels = FALSE)

# estimates the average distance between clusters
fviz_nbclust(data, pam, method = "silhouette")

# Looking to see when the rate of the increase of the gap statistic slows down
# Ideally we would try to maximize the gap statistic
fviz_nbclust(data, pam, method = "gap_stat")

# Seems like K = 6 is a good point to cut off

# To double check, I like to see the number of clusters in PCA
# It hints at the variablity of data

fviz_pca_ind(prcomp(data), 
             title = "PCA",
             palette = "jco",
             geom = "point",
             legend = "bottom")
# No conclusive data, but there seems to be a general/consistent spread 


clusternum <- NbClust((data), distance = "euclidean", method = "kmeans")

# Clustering
pam.res4 <- pam(data, 4, metric = "euclidean", stand = FALSE)
pam.res3 <- pam(data, 3, metric = "euclidean", stand = FALSE)

fviz_cluster(pam.res3, data = data, geom = "point", 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE)
fviz_cluster(pam.res4, data = data, geom = "point", 
             ellipse.type = "euclid", star.plot = TRUE, repel = TRUE)
# 4 clusters seem to be better

# fviz_silhouette(pam.res5, palette = "jco")
fviz_silhouette(pam.res4, palette = "jco")

# 4 clusters is best

hier.clust
clusters.cut <- cutree(hier.clust, k = 5)

# group 1: from hierarchical cluster
hc.group <- subset(data, clusters.cut == 1)
hc.group.genes <- noquote(row.names(hc.group))

# group 2: from k-means
kmeans.group <- subset(data, pam.res4$clustering == 1)
kmeans.group.genes <- noquote(row.names(kmeans.group))

# Outputting the genes for analysis
# Setting Directory to source file (deleting this when moving code to Markdown)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Making text files to copy list to DAVID
write(hc.group.genes, "hc.group.txt")
write(kmeans.group.genes, "kmeans.group.txt")

# Pulling Data from DAVID analysis (pulled from downloaded txt files)
hc.group.output.df <- read.csv("hc.group.output.txt", sep = "\t")
kmeans.group.output.df <- read.csv("kmeans.group.output.txt", sep = "\t")

# Hierarchical cluster analysis (checking proteins / middle column)
hc.group.output.df$Name

# K-means analysis (checking proteins / middle column)
kmeans.group.output.df$Name



# GSEA --------------------------------------------------------------------

library(clusterProfiler)

res <- res[res$baseMean > 50,]

res <- res[order(-res$stat),]
gene.list <- res$stat
names(gene.list) <- row.names(res)

# The GSE analysis function
gse <- gseGO(gene.list,
             ont = "BP", 
             keyType = "ENTREZID",
             OrgDb = "org.Hs.eg.db",
             eps = 1e-300)
gse.df <- as.data.frame(gse)

# Top 5 gene sets
head(gse)

gseaplot(gse, geneSetID = 1) #better position for enrichment score, more left leaning

gseaplot(gse, geneSetID = 2)


# Cytoscape / Network -----------------------------------------------------

genematrix <- countsData

#Check For Outliers
IAC=cor(genematrix,use="p")
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
cluster1=hclust(as.dist(1-IAC))

topN <- 1000
potential_biomarkers <- sorted_analysis_table[1:topN,]

keepGenesExpr = rank(-rowMeans(genematrix))<=1000
# keepGenesExpr <- rownames(potential_biomarkers)
mydata<-genematrix[rownames(genematrix) %in% rownames(potential_biomarkers),]
#we need to transpose our matrix for the format WGCNA expects
dataExpr <- as.matrix(t(mydata))
dataExpr <- t(mydata)

dataExpr[1:6, 1:1000] <- as.numeric(as.character(dataExpr[1:6, 1:1000]))
library(WGCNA)
net = blockwiseModules(dataExpr, power = 5,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold =10, mergeCutHeight = 0.2,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase="TOM", verbose=3, ds=3)

#now we get to see our modules!
mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
#save modules
genes=colnames(dataExpr)
moduleColors=labels2colors(net$colors)
mymodules=cbind(genes,moduleColors)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# save(mymodules, file = "modules.txt")
#correlate with phenotype
phenotype=list(0,0,0,1,1,1)
#get Eigengenes
MEs0=moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, phenotype, use = "p")
moduleTraitPvalue=corPvalueStudent(moduleTraitCor, 11)


#Select the gene modules
genes = colnames(dataExpr)

#if you want export specific colors, substitute the second modulecolors by above modules
inModule = is.finite(match(moduleColors, moduleColors))
modGenes = genes[inModule]

#Select the corresponding topologic overlap 
adjacency = adjacency(dataExpr, power = 9, type = "unsigned")
TOM = TOMsimilarity(adjacency)

modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

exportNetworkToCytoscape(
  modTOM,
  edgeFile = "edge.txt",
  nodeFile = "node.txt",
  weighted = TRUE,
  threshold = 0.75,
  nodeNames = modGenes,
  nodeAttr = moduleColors[inModule])


whole.network.results <- read.csv("Whole_Network_node.csv", sep = ",")
whole.network.results <- as.data.frame(whole.network.results)
filtered.network.results <- whole.network.results[,c(30,31,23,12,13,14,15,17)]

# table of central genes per cluster
filtered.network.results[filtered.network.results$MCODE..Node.Status..1. == "Seed",]

# Sort by neighborhood connectivity and 
neighborhood.conn <-
  filtered.network.results[order(filtered.network.results$MCODE..Score..1.,
                                 decreasing = TRUE), ]
head(neighborhood.conn, 15)

# MCODE score (interconnectivity)
mcode.score <-
  filtered.network.results[order(filtered.network.results$NeighborhoodConnectivity), ]
head(mcode.score, 15)
