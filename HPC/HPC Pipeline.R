###################################################################################################
# Playground

setwd("~/data/IBP_data")

require(org.Ce.eg.db, quietly = T) # for c elegans
require(GO.db)
require(scde, quietly = T)
require(parallel, quietly = T)
require(Matrix)
library(KEGGREST)

getPossibleArguments()

#######################################################################################################################
### Pipeline arguments 
data.path = "./GenomicData/Experiment1.RData"
UMI.type = "tissue"
gene.type = "symbol"
UMI.subset = c("Body wall muscle", "Neurons", "Hypodermis")
n.cores = parallel::detectCores()
gene.set.annot = "GO"
Plot = F
min.reads = 1
min.detected = 3
min.lib.size = 100
number.cell = 600


#######################################################################################################################


#######################################################################################################################
# pagoda.analysis
# INPUT: 
#   - UMI.type:   string indicating whether column must be coded as a cell type ("cell"), as a tissue type ("tissue"), 
#                 or just keep the UMI encoding ("", default).
#   - gene.type:    a string indication which gene annotation to use. Default ("") uses the WormBase annotation, 
#                   "symbol" uses the classical symbol annotation.
#   - UMI.subset:   a vector of character indicating which cell type/tissue type/UMI code to keep for the analysis. 
#                   NA(default) performs the analysis on the whode dataset.
#   - n.cores:    the number of cores to use to do the computing
#   - gene.set.annot:   either a string ("KEGG" for KEGG pathways, or "de novo" for de novo gene set analysis) or a list 
#                       (names = gene set names, values = gene names). Any other value will take the GO annotation as default
#   - Plot:     whether to plot intermediate results or not 
#   - min.lib.size:   Minimum number of genes detected in a cell. Cells with fewer genes will be removed (default: 100)
#   - min.reads:      Minimum number of reads per gene. Genes with fewer reads will be removed (default: 1)
#   - min.detected:   Minimum number of cells a gene must be seen in. Genes not seen in a sufficient number of cells 
#                     will be removed (default: 3)

#-------------------------------------------------------------
### STEP1: Prepare the data 

# Get Cao et al. 2017 dataset
dataset <- readRDS(data.path)
cat("Data loaded...\n")

# Get the column names
if(UMI.type == "cell"){
  cell.names <- dataset$cell.names[colnames(dataset$exp.data),"cell.type"]
  cell.names <- paste0(cell.names,"-",1:length(cell.names)) # Avoid duplicate names
} else if (UMI.type == "tissue"){
  cell.names <- dataset$cell.names[colnames(dataset$exp.data),"tissue"]
  cell.names <- paste0(cell.names,"-",1:length(cell.names)) # Avoid duplicate names
} else if (UMI.type == "neuron"){
  cell.names <- dataset$cell.names[colnames(dataset$exp.data),"neuron.type"]
  cell.names <- paste0(cell.names,"-",1:length(cell.names)) # Avoid duplicate names
} else{
  cell.names <- colnames(dataset$exp.data)
}

# Get the gene annotation
if(gene.type == "symbol"){
  gene.names <- as.character(dataset$gene.names[rownames(dataset$exp.data),])
} else {
  gene.names <- as.character(dataset$gene.names)
}

# Subset the data (if needed)
if (!is.na(UMI.subset[1])){
  if(is.character(UMI.subset)){
    if (length(UMI.subset) < 2) { stop("You cannot subset for less than 2 cell/tissue types !") }
    ind <- which( gsub("-\\d+", "", cell.names) %in% UMI.subset)
    dataset$exp.data <- dataset$exp.data[,ind]
    cell.names <- cell.names[ind]
  } else if (is.numeric(UMI.subset)){
    ind <- sample(1:ncol(dataset$exp.data), size = UMI.subset, replace = F)
    dataset$exp.data <- dataset$exp.data[,ind]
    cell.names <- cell.names[ind]
  }
}

# If an annotation list is given, keep only the genes from that list
if (is.list(gene.set.annot)){
  genes.interest <- unique(unlist(lapply(gene.set.annot, function(x){return(x)})))
  ind <- which(gene.names %in% genes.interest)
  dataset$exp.data <- dataset$exp.data[ind,]
  gene.names <- gene.names[ind]
}
dimension <- dim(dataset$exp.data)

clean.data <- dataset$exp.data
dimnames(clean.data) <- list(gene.names, cell.names)

# Data cleaning: remove poor cells and genes
# Clean for genes
clean.data <- clean.data[rowSums(clean.data) > min.reads, ]
clean.data <- clean.data[rowSums(clean.data > 0) > min.detected, ]
clean.data <- clean.data[, colSums(clean.data > 0) > min.lib.size]

# Clean for cells 
# Problem the library size frequency distribution is not the same bor both tissues
# So, sample the same number of cells for each tissue
ind.1 <- sample(which(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[1]), size = min(number.cell, sum(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[1])), replace = F)
ind.2 <- sample(which(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[2]), size = min(number.cell, sum(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[2])), replace = F)
ind.3 <- sample(which(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[3]), size = min(number.cell, sum(gsub("-\\d*","",colnames(clean.data)) %in% UMI.subset[3])), replace = F)
clean.data <- clean.data[,c(ind.1, ind.2, ind.3)]

dimension <- dim(clean.data)
dimension.names <- dimnames(clean.data)

# Coerce the dgCMatrix to an integer matrix
clean.data <- as.integer(as.matrix(clean.data))
dim(clean.data) <- dimension
dimnames(clean.data) <- dimension.names
rm(dataset) # free memory space 



#-------------------------------------------------------------
### STEP2: Error model fitting
cat("Fitting error models...\n")
cell.group <- gsub("-\\d+", "", colnames(clean.data))
if(chunk != 0){
  i <- 1
  while(i < ncol(clean.data)){
    ind <- i:(i+chunk-1)
    knn.temp <- knn.error.models(clean.data[,ind], groups = cell.group[ind], k = ncol(clean.data)/4, n.cores = n.cores, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 20, min.size.entries = 800, verbose = 0, save.model.plots = Plot)
    if(i == 1){
      knn <- knn.temp
    } else {
      knn <- rbind(knn, knn.temp)
      attr(knn, "groups") <- c(attr(knn,"groups"), attr(knn.temp,"groups"))
      rm(knn.temp)
    }
    i <- i + chunk
    print(i)
  }
} else {
  knn <- knn.error.models(clean.data, groups = cell.group, k = ncol(clean.data)/4, n.cores = n.cores, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 20, min.size.entries = 800, verbose = 0, save.model.plots = Plot)
}
cat("Models where fitted",ifelse(Plot,"(and plotted; see cell.models.pdf in wd)",""),"...\n")


#-------------------------------------------------------------
### STEP3: Normalizing variance and controlling for sequence depth
cat("Normalizing variance...\n")
varinfo <- pagoda.varnorm(knn, counts = clean.data, trim = 3/ncol(clean.data), max.adj.var = 5, n.cores = n.cores, plot = Plot, verbose = F) # Normalizing variance
varinfo <- pagoda.subtract.aspect(varinfo, colSums(clean.data[, rownames(knn)]>0)) # controlling for sequencing depth
cat("Variance normalized...\n")


#-------------------------------------------------------------
### STEP4: weighted PCA
# Get the gene sets
if (gene.set.annot == "GO"){
  # default is using the GO annotation
  require(org.Ce.eg.db, quietly = T) # for c elegans
  ids <- unlist(lapply(mget(rownames(clean.data), org.Ce.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
  rids <- names(ids); names(rids) <- ids 
  gos.interest <- ls(org.Ce.egGO2ALLEGS)
  gl.env <- lapply(mget(gos.interest, org.Ce.egGO2ALLEGS), function(x){ unique(as.character(na.omit(rids[x]))) }) 
} else if(gene.set.annot == "KEGG"){
  load("./Genomic Data/kegg.Ce.glist.RData")
  gl.env <- kegg.Ce.glist
} else if (is.list(gene.set.annot)){
  gl.env = gene.set.annot
} else { 
  stop("You did not provide a correct gene set annotation ")
} 
# Clean the gene set based on the number of genes in each set
gl.env <- clean.gos(gl.env, min.size = 5, max.size = 1000, annot = F) # remove GOs with too few or too many genes
gl.env <- list2env(gl.env) # convert to an environment
cat("Gene list loaded using the ", gene.set.annot, " gene set annotation ...\n")

# Weighted first principal component magnitudes for each gene set in the provided environment
cat("wPCA computing...\n")
system.time(pwpca <- pagoda.pathway.wPCA(varinfo, gl.env, n.components = 1, n.cores = n.cores))
cat("wPCA completed...\n")

# De novo gene clustering (if needed)
if (gene.set.annot == "de novo"){
  clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 50, n.cores = n.cores, plot = Plot)
} else {
  clpca <- NULL
}

#-------------------------------------------------------------
### STEP5: output visualization

cat("Started visualization...\n")
tam <- try(pagoda.top.aspects(pwpca, clpca, return.table = F, plot = Plot, z.score = 1.96))
if(class(tam) == "try-error"){
  stop("No overdispersed gene set!")
}

# Reduce the redundacy in the results
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca=clpca)

if(nrow(tamr$xv) < 2){
  tamr2 <- tam 
} else {
  tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = F, cell.clustering = hc)
}

# Cell clustering
hc <- pagoda.cluster.cells(tamr2, varinfo)


if(gene.set.annot == "GO"){
  require("GO.db")
  sel.go.id <- rownames(tamr2$xv)
  go.id <- gsub("#PC1# ", "", rownames(tamr2$xv))
  go.term <- Term(GOTERM)[go.id]
  rownames(tamr2$xv) <- go.term 
}


if(Plot){
  par(xpd=T)
  l2cols <- c("coral", "lightblue")[as.integer(factor(gsub("-\\d*","", hc$labels), levels = UMI.subset))]
  pagoda.view.aspects(tamr2, cell.clustering = hc, box = T, labCol = NA, margins = c(0.5, 20), col.cols = l2cols, cexRow=1.5,xpd=T)
}

cat("Done !\n")

# getPossibleArguments
# This function displays on the console the possible cell types or tissues types for the pipeline.
getPossibleArguments <- function(){
  subset.arg <- readRDS("./GenomicData/Experiment1.RData")
  cell.names <- subset.arg$cell.names
  cell.names[is.na(cell.names)] <- "NA"
  
  cell.types <- sapply(unique(cell.names[,1]), function(type){sum(cell.names[,1] == type)})
  cell.types <- cell.types[order(names(cell.types))]
  tissue.types <- sapply(unique(cell.names[,2]), function(type){sum(cell.names[,2] == type)})
  tissue.types <- tissue.types[order(names(tissue.types))]
  neuron.types <- sapply(unique(cell.names[,3]), function(type){sum(cell.names[,3] == type)})
  neuron.types <- neuron.types[order(names(neuron.types))]
  
  
  cat("The following cell types can be given:\n\n")
  cat("\tAmount\tCell Name\n\n")
  for (i in 1:length(cell.types)){
    cat(paste0("\t", cell.types[i], "\t", names(cell.types[i]),"\n"))
  }
  cat("==============================\n\n")
  
  cat("The following tissue types can be given:\n\n")
  cat("\tAmount\tTissue Name\n\n")
  for (i in 1:length(tissue.types)){
    cat(paste0("\t", tissue.types[i], "\t", names(tissue.types[i]),"\n"))
  }
  cat("==============================\n\n")
  
  cat("The following neuron types can be given:\n\n")
  cat("\tAmount\tNeuron Name\n\n")
  for (i in 1:length(neuron.types)){
    cat(paste0("\t", neuron.types[i], "\t", names(neuron.types[i]),"\n"))
  }
  cat("==============================\n\n")
}