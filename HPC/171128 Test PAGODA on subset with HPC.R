setwd("~/data/IBP_data/")
# setwd("/media/christophe/Windows/2nd stage - shortcut/1st semester/Integrated Bioinformatics Project/Project/Genomic Data/")

library(scde)
library(Matrix)
# install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")
# install.packages("/media/christophe/Windows/2nd stage - shortcut/1st semester/Integrated Bioinformatics Project/Packages/scde-1.99.2.tar.gz", repos = NULL, type="source")

##################################################################
# Step 1: Prepare the data

### Preparing the data

subset.file <- readRDS("./GenomicData/ASER-ASEL-subset.RData")
names.temp <- list(dimnames(subset.file$data)[[1]], paste0(dimnames(subset.file$data)[[2]],1:ncol(subset.file$data)))
dim.temp <- dim(subset.file$data)
subset.data <- as.integer(as.matrix(subset.file$data))
dim(subset.data) <- dim.temp
dimnames(subset.data) <- names.temp
subset.gnames <- subset.file$gene.names
subset.cnames <- subset.file$phenotype
# change the gene name according to the gene names used in the upcoming gene set analysis 
rownames(subset.data) <- subset.gnames[rownames(subset.data),2]

# Default: cd <- clean.counts(as.matrix(subset.data), min.lib.size = 1800, min.reads = 10, min.detected = 5) deletes everything
clean.data <- clean.counts(counts = subset.data, min.lib.size = 100, min.reads = 1, min.detected = 3) 
# So here keep: 
#   - Cells with at least 100 gene counts
#   - Genes with at least 1 read
#   - Genes must bbe seen in at least 3 cells 


##################################################################
# Step 2: Fitting error models

cell.group <- gsub("\\d+", "", colnames(clean.data))
# First change the gene name according to the gene names used in the upcoming gene set analysis 
system.time(knn <- knn.error.models(clean.data, groups = cell.group, k = ncol(clean.data)/4, n.cores = 4, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10, min.size.entries = 800, verbose = 1))
# min.size.entries had to be changed from 2000 to 800, because not enough genes are present in the cleaned data
# elapsed 10.6 s with 3 cores
# ended with warning: In stats::cor(sqrt(matrix(as.numeric(as.matrix(ca[, ids])), nrow = nrow(ca),:the standard deviation is zero


##################################################################
# Step 3: Normalizing variance and controlling for sequence depth

varinfo <- pagoda.varnorm(knn, counts = clean.data, trim = 3/ncol(clean.data), max.adj.var = 5, n.cores = 1, plot = TRUE) # this makes a plot
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]

# controlling for sequence depth
varinfo <- pagoda.subtract.aspect(varinfo, colSums(clean.data[, rownames(knn)]>0))


##################################################################
# Step 4: Evaluate overdispersion of pre-defined gene sets

# library(org.Hs.eg.db) for human
library(org.Ce.eg.db) # for c elegans
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(clean.data), org.Ce.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
rids <- names(ids); names(rids) <- ids 
# convert GO lists from ids to gene names
gos.interest <- ls(org.Ce.egGO2ALLEGS)
go.env <- lapply(mget(gos.interest, org.Ce.egGO2ALLEGS), function(x) as.character(na.omit(rids[x]))) 
go.env <- clean.gos(go.env, min.size = 5, max.size = 1000, annot = F) # remove GOs with too few or too many genes
go.env <- list2env(go.env) # convert to an environment
length(go.env)

# Now, we can calculate weighted first principal component magnitudes for each GO gene set in the 
# provided environment
system.time(pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 4)) # instead of using clean.gos() we can add arguments min.pathway.size and max.pathway.size
# Elapsed 61;72 s with 4 cores and 1445 gene sets

df <- pagoda.top.aspects(pwpca, return.table = T, plot = TRUE, z.score = 1.96) # this makes a plot
tam <- pagoda.top.aspects(pwpca, clpca = NULL, return.table = F, plot = F, z.score = 1.96, n.cells = NULL, adjust.scores = T, score.alpha = 0.5) 
# The z column gives the Z-score of pathway over-dispersion relative to the genome-wide model (Z-score 
#     of 1.96 corresponds to P-value of 5%, etc.).
# "z.adj" column shows the Z-score adjusted for multiple hypothesis (using Benjamini-Hochberg correction).
# "score" gives observed/expected variance ratio
# "sh.z" and "adj.sh.z" columns give the raw and adjusted Z-scores of "pathway cohesion", which compares 
#     the observed PC1 magnitude to the magnitudes obtained when the observations for each gene are 
#     randomized with respect to cells. When such Z-score is high (e.g. for GO:0008009) then multiple 
#     genes within the pathway contribute to the coordinated pattern.

### Conclusion 
library(GO.db)
go.term <- Term(GOTERM)[df$name]
df$go.term <- go.term
df
#write.csv(df, file = "./Genomic Data/171020 PAGODA ASER-ASEL overdisp GO.csv")
# name npc  n   score        z    adj.z sh.z adj.sh.z
# 264 GO:0009975   1 14 1.92804 4.018418 2.361085   NA       NA
# 343 GO:0016849   1 14 1.92804 4.018418 2.361085   NA       NA
# 46  GO:0004383   1 14 1.92804 4.018418 2.361085   NA       NA

# GO:0009975  cyclase activity 
# GO:0016849  phosphorus-oxygen lyase activity
# GO:0004383  guanylate cyclase activity

##################################################################
# Step 5: Visualize significant aspects of heterogeneity

# translate group and sample source data from Pollen et al. into color codes. 
l2cols <- c("coral4", "slateblue")[as.integer(factor(gsub("\\d","", colnames(clean.data)), levels = c("ASER", "ASEL")))]

hc <- pagoda.cluster.cells(tam, varinfo)

pagoda.view.aspects(tam, cell.clustering = hc, box = T, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))

# Note: because the 3 GO terms include the same genes and therefore are the same gene set, this leads to a single 
#   GO. Therefore, the subsequent steps are useless.
tamr <- pagoda.reduce.loading.redundancy(tam, pwpca)
#tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = df$name, labCol = 1:60, box = TRUE, margins = c(10.5, 10.5), trim = 0)
tamr$df[tamr$df$adj.z>1.96,]

############################################################
# END


