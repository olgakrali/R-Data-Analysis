
library(SingleCellExperiment)
library(scater)
library(scran)  
library(pcaMethods)
library(SC3)
library(pheatmap)
library(mclust)
library(mvoutlier)
library(limma)
library(scRNA.seq.funcs)
library(sva) 
library(matrixStats)
library(M3Drop)
library(RColorBrewer)
library(TSCAN)
library(scImpute) 


data <- read.table("dataset.csv", header = TRUE, sep="\t") 
labelb<-read.table("LabelsBatch1.csv", header = TRUE, sep=",")

n <- data[,1] #first column contains the names of the cells
#transpose data 
data1 <- as.data.frame(t(data[,-1]))
dim(data1)

colnames(data1) <-n
colnames(data1) #cell names
rownames(data1) #genes


counts <- as.matrix(data1) 
class(counts) <- "numeric" #otherwise the matrix will contain characters

sce <- SingleCellExperiment(
  assays = list(counts = counts),
  colData = labelb
)

expressed_genes <- rowSums(counts(sce) > 0) > 0  #drop genes without expression in any cell
sce <- sce[expressed_genes, ] 


#Expression QC
#Plot of genes with higher expression 
#Add spike in as control genes

erccgenes <- grep ("^ERCC", rownames(data1),data1, value = TRUE)

#find MT genes in our dataset
mitochondria <- grep ("^MT", rownames(data1), data1, value = TRUE)



isSpike(sce, "ERCC") <- rownames(sce) %in% erccgenes
isSpike(sce, "MT") <- rownames(sce) %in% mitochondria
  
#ntroduce spike in in the single cell experiment
reads <- calculateQCMetrics(
  sce,
  feature_controls = list(
    ERCC = isSpike(sce, "ERCC"),
    MT = isSpike(sce, "MT")
  )
)


#Distribution of total counts, and number of genes
hist(
  reads$total_counts,
  breaks = 100
)

abline(v = 280000, col = "red") #right of this line the total counts frequency is relatively higher
filter_by_total_counts <- (reads$total_counts > 280000)
table(filter_by_total_counts)

hist(
  reads$total_features,              #most of the cells have 3000-5000 detected genes
  breaks = 100
)
abline(v = 1800, col = "red")

filter_by_expr_features <- (reads$total_features > 1800)
table(filter_by_expr_features)


plotQC(reads, type = "highest-expression") #identifies the top 50 genes from which the most
#reads come from


#relationship between spike ins(ERCC, MT) and endogenous RNAs 

plotPhenoData(
  reads,
  aes_string(
    x = "total_features",          
    y = "pct_counts_MT",
    colour = "batch"
  )
)

filter_by_MT <- reads$pct_counts_MT < 0.8  # remove cells with unusual number of reads in MT genes
table(filter_by_MT)


plotPhenoData(
  reads,
  aes_string(
    x = "total_features",
    y = "pct_counts_ERCC",
    colour = "batch"
  )
)

filter_by_ERCC <- reads$pct_counts_ERCC < 0.25

#Use the filters created above and store it as a new variable (we did manual manipulation)
#to the filters, maybe it is not 100% accurate)

#Filter cells, I leaves 416/446 cells
reads$use <- (
  filter_by_expr_features &
    filter_by_total_counts &
    filter_by_ERCC &
    filter_by_MT
)

table(reads$use)

#How cells are distributed after cell filtering

reads <- plotPCA(
  reads,
  size_by = "total_features",
  shape_by = "use",
  pca_data_input = "pdata",
  detect_outliers = TRUE,
  return_SCE = TRUE
)

table(reads$outlier) #identifies 1 outlier. Cannot remove it. Probably the authors have already removed the outliers

#find differences between automatic and manual way to identify outliers and filter cells
auto <- colnames(reads)[reads$outlier]
man <- colnames(reads)[!reads$use]
venn.diag <- vennCounts(
  cbind(colnames(reads) %in% auto,
        colnames(reads) %in% man)
)
vennDiagram(                         
  venn.diag,
  names = c("Automatic", "Manual"),  #which one do we trust?
  circle.col = c("red", "green")
)




#Filters detectable genes (more than 1 transcript in at least 2 cells) (18517)

filter_genes <- apply(
  counts(reads),
  1,
  function(x) length(x[x > 1]) >= 2
)
rowData(reads)$use <- filter_genes

table(filter_genes) #some genes are undetectable (3110) 

dim(reads[rowData(reads)$use, colData(reads)$use])


#create new assay
assay(reads, "logcounts_raw") <- log2(counts(reads) + 1)

logcounts (reads) <- log2(calculateCPM(reads, use.size.factors = FALSE) + 1)

reducedDim(reads) <- NULL



##QC filtering (continue analysis with filtering; 18517 genes and 416 cells)##

reads.qc <- reads[rowData(reads)$use, colData(reads)$use] #filtered genes and cells respectively


genes <- rowSums(counts(reads.qc)>0)>0   #apply this, in case there is again zero gene 
reads.qc2 <- reads.qc[genes, ]          # expression after applying QC (it excludes 4 more genes)
endog_genes <- !rowData(reads.qc2)$is_feature_control

plotPCA(                    #before the quality control
  reads[endog_genes, ],
  exprs_values = "counts",  
  colour_by = "cell_type",
  size_by = "total_features",
  shape_by = "tissue"
)

plotPCA(
  reads[endog_genes, ],
  exprs_values = "logcounts_raw",  #better distinction of the different cell types than counts
  colour_by = "cell_type",        #but better try to use logcounts instead
  size_by = "total_features",
  shape_by = "tissue"
)

plotPCA(                      #after quality control
  reads.qc2[endog_genes, ],
  exprs_values = "counts",        #repeat with logcounts, and logcounts_raw
  colour_by = "cell_type",
  size_by = "total_features",
  shape_by = "tissue"
)


#tSNE plots (what could be the better value for perplexity?)


plotTSNE(            
  reads[endog_genes, ],
  exprs_values = "logcounts_raw",   #check counts and logcounts(cpm)
  perplexity = 130,
  colour_by = "cell_type",
  size_by = "total_features",
  shape_by = "tissue",
  rand_seed = 123456,
  check_duplicates = FALSE
)



plotTSNE(             
  reads.qc2[endog_genes, ],
  exprs_values = "logcounts_raw",   #check counts and logcounts(cpm)
  perplexity = 130,
  colour_by = "cell_type",
  size_by = "total_features",
  shape_by = "tissue",
  rand_seed = 123456,
  check_duplicates = FALSE
)

# identify correlations between principal components and total features
#total features
plotQC(
  reads.qc2[endog_genes, ],
  type = "find-pcs",
  exprs_values = "logcounts_raw",
  variable = "total_features"
)


plotQC(                         #identifies the best explanatory variables
  reads.qc2[endog_genes, ],
  type = "expl",
  exprs_values = "logcounts_raw",
  variables = c(
    "total_features",
    "total_counts",
    "batch",
    "tissue",
    "cell_type",
    "pct_counts_MT"
  )
)

#Combat batch effect removal 

combat_data <- logcounts(reads.qc2)
mod_data <- as.data.frame(t(combat_data))

# Basic batch removal
mod = model.matrix(~ 1, data = mod_data)

assay(reads.qc2, "combat") <- ComBat(
  dat = t(mod_data),
  batch = factor(reads.qc2$batch),
  mod = mod,
  par.prior = TRUE,
  prior.plots = FALSE
  
)



#check all the assays (counts, log, logcounts_raw) for the best explanatory variables

for(n in assayNames(reads.qc2)) {
  print(
    plotQC(
      reads.qc2[endog_genes, ],
      type = "expl",
      exprs_values = n,
      variables = c(
        "total_features",
        "total_counts",
        "batch",
        "tissue",
        "cell_type",
        "pct_counts_MT"
      )
    ) +
      ggtitle(n)
  )
}



plotPCA(
  reads.qc2[endog_genes, ],
  exprs_values = "combat",  #check difference with counts
  colour_by = "cell_type",
  size_by = "age",
  shape_by = "tissue"
)

plotTSNE(               
  reads.qc2,
  exprs_values = "combat",   #check counts and logcounts(cpm)
  perplexity = 130,
  colour_by = "cell_type",
  size_by = "age",
  shape_by = "tissue",
  rand_seed = 123456,
  check_duplicates = FALSE
)

combat <- assays(reads.qc2)$combat

#write.csv(combat,"Combatbatch.csv")



#6) Clustering

#After combat correction (use the combat data from the read.qc assay)

set.seed(1234567)

#Create a new sce1, where counts = combat, and logcounts = combat_data (logcounts(reads.qc2))
 
labels <- colData(reads.qc2)

sce1 <- SingleCellExperiment(
  assays = list(counts = combat),
  colData = labels
)

assay(sce1, "logcounts") <- combat_data

plotPCA(sce1, colour_by = "cell_type")

#estimate the number of clusters (we already know that the cell types are 9)
sce1 <- sc3_estimate_k(sce1)
metadata(sce1)$sc3$k_estimation #I got 12 clusters


rowData(sce1)$feature_symbol <- rownames(sce1) #as feature symbols the gene names
sce1<- sc3(sce1, ks = 12, biology = TRUE, gene_filter = FALSE) #remove gene filtering, it does not work otherwise (no other zeros contained)


#result interpretation
sc3_plot_consensus(sce1, k = 12, show_pdata = "cell_type")
sc3_plot_consensus(sce1, k = 12, show_pdata = "age") #let's try with age

sc3_plot_silhouette(sce1, k =12)

sc3_plot_expression(sce1, k = 12, show_pdata = "cell_type")
sc3_plot_expression(sce1, k = 12, show_pdata = "age")


#find marker genes

sc3_plot_markers(sce1, k = 12, show_pdata = "cell_type")

#New PCA plot based on SC3 clustering

plotPCA(sce1, colour_by = "sc3_12_clusters")



#Comparison of the clustering with the original dataset (cell types)
# 0.67, it is relatively good

adjustedRandIndex(colData(sce1)$cell_type, colData(sce1)$sc3_12_clusters)

#t-SNE and k means clustering


sce1 <- plotTSNE(sce1, rand_seed = 1, return_SCE = TRUE)

#use the suggested 12 clusters found by the first clustering (SC3)
colData(sce1)$tSNE_kmeans <- as.character(kmeans(sce1@reducedDims$TSNE, centers = 12)$clust)
plotTSNE(sce1, rand_seed = 1, colour_by = "tSNE_kmeans")

#comparison between tSNE_kmeans and clusters given by the authors
adjustedRandIndex(colData(sce1)$cell_type, colData(sce1)$tSNE_kmeans) #0.4893

#Repeat with the 9 cell types (clusters) from the dataset


colData(sce1)$tSNE_kmeans <- as.character(kmeans(sce1@reducedDims$TSNE, centers = 9)$clust)
plotTSNE(sce1, rand_seed = 1, colour_by = "cell_type")


#SINCERA (hierarchical clustering) based on combat data (22 clusters!)

#as an imput
dat <- apply(combat, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")

num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)

pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = kk,
  kmeans_k = 100,
  show_rownames = FALSE
)



####### Feature selection (new sce without spike ins)

#Quality control can be achieved without the previous methods with M3Drop library
set.seed(1)


#start with reads

cells <- colData(reads)$cell_type

scenew <- M3DropCleanData(
  counts(reads),
  labels = cells,
  min_detected_genes = 100,
  is.counts = TRUE
)

exp_matrix <- scenew$data   #(18768 detectable genes, sligtly larger than QC method)
celllabels <- factor(scenew$labels)
cell_col <- brewer.pal(max(3,length(unique(celllabels))), "Set3")





# HVG (High Variable Genes with Brennecke method)


#pink dots are all these genes with significant biological variability

Brennecke_HVG <- BrenneckeGetVariableGenes(     #use the exp matrix from the previous step
       exp_matrix,
       fdr = 0.01,
       minBiolDisp = 0.5
   )

length(Brennecke_HVG) #1812 HVGs

### Correlated expression (try with another computer)

correlmatrix <- cor(t(exp_matrix), method = "spearman") #for cor between genes; memory limit exceeds!!!

diag(correlmatrix) <- rep(0, times=nrow(exp_matrix))
score <- apply(correlmatrix, 1, function(x) {max(abs(x))}) 
names(score) <- rownames(exp_matrix);
score <- score[order(-score)]
Correlgenes = names(score[1:1500])

# Feature selection with PCA loadings ###

pcaload <- prcomp(log(exp_matrix + 1)/log(2))
plot(pcaload$rotation[,1], pcaload$rotation[,2], pch =16, col = cell_col[as.factor(celllabels)])

score <- rowSums(abs(pcaload$x[,c(1,2)]))
names(score) <- rownames(exp_matrix)
score <- score[order(-score)]
PCA_genes = names(score[1:1500]) 



#9) Pseudo-time construction  (I do not think it helps us somehow!)

set.seed(1)
#TSCAN   (a combination of clustering and pseudotime analysis)

reads
cells <- colData(reads)$cell_type
lcoun <- logcounts(reads)
colnames(lcoun) <- cells

proclcoun <- TSCAN::preprocess(lcoun)
colnames(proclcoun) <- 1:ncol(lcoun)
lcounclust <- TSCAN::exprmclust(proclcoun, clusternum = 9)
TSCAN::plotmclust(lcounclust) #graph 

lcounorderTSCAN <- TSCAN::TSCANorder(lcounclust, orderonly = F)
pseudotime_order_tscan <- as.character(lcounorderTSCAN$sample_name)

cells[lcounclust$clusterid == 9]


colours <- rainbow(n = 9) 
tmp <-
  factor(
    cells[as.numeric(pseudotime_order_tscan)],
    levels = c("neurons", "oligodendrocytes", "hybrid", "OPC", "fetal_replicating", "fetal_quiescent", "astrocytes", "endothelial", "microglia")
  )

plot(
  as.numeric(tmp),
  xlab = "Pseudotime Order",
  ylab = "Timepoint",
  col = colours[tmp],
  pch = 16
)

#we have different cell types not diffent different cell states (I do not know how we could utilize that process in this data)

colnames(lcoun) <- 1:ncol(lcoun)
TSCAN::singlegeneplot(
  lcoun[rownames(lcoun) == "ACP5", ],   #trace a gene from one cell type to the other
  lcounorderTSCAN
)

#### Imputation ####
## Gene drop outs###

# work on reads single cell experiment (unmodified data)

write.csv(counts(reads), "reads.csv")
scimpute(
  count_path = "reads.csv",
  infile = "csv",
  outfile = "txt",
  out_dir = "./",
  Kcluster = 9,
  ncores = 1          #ncores=2 is not supported by windows
)

#Let's inspect the data

imp <- read.table("scimpute_count.txt")
colnames(imp) <- NULL
imp <- SingleCellExperiment(
  assays = list(logcounts = log2(as.matrix(imp) + 1)),
  colData = colData(reads)
)
plotPCA(
  imp,
  colour_by = "cell_type"
)

#How imputed data will respond to clustering
imp <- sc3_estimate_k(imp)
metadata(imp)$sc3$k_estimation   #10 clusters (closer to the cell types number)

rowData(imp)$feature_symbol <- rownames(imp)

imp <- sc3(imp, ks = 10, gene_filter = FALSE)

#clusters (sc3 derived) comparison with cell types (lower value the previous part (sc3 on the transformed data: QC, batch effect removal))
adjustedRandIndex(colData(reads)$cell_type, colData(imp)$sc3_10_clusters) #0.58


plotPCA(
  imp,
  colour_by = "sc3_10_clusters"
)



