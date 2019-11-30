setwd("/media/christophe/Windows/2nd stage - shortcut/1st semester/Integrated Bioinformatics Project/Project/HPC/1205-NeuronvsMusclevsHypodermis-1500cells/")
setwd("C:/Users/lechr/Documents/KUL/2nd stage/1st semester/Integrated Bioinformatics Project/Project/HPC/1205-NeuronvsMusclevsHypodermis-1500cells")
load("1205-fullResults.RData")

library(scde)

#-------------------------------------------------------------
### STEP5: output visualization

if(gene.set.annot == "GO"){
  require("GO.db")
  sel.go.id <- rownames(tamr2$xv)
  go.id <- gsub("#PC1# ", "", rownames(tamr2$xv))
  go.term <- Term(GOTERM)[go.id]
  go.size <- tamr2$df[tamr2$df$name %in% go.id,]$n
  rownames(tamr2$xv) <- paste("  ",go.term, " (",go.size,")", sep="")
}

go.term
nb.cells <- sapply(unique(gsub("-\\d*","",colnames(tamr2$xv))), function(x){sum(gsub("-\\d*","",colnames(tamr2$xv))==x)})

png("GO-HypMusNeur-1500-2.png", height = 20, width = 40, units = "in", res = 300)
l2cols <- c("coral", "lightblue", "bisque")[as.integer(factor(gsub("-\\d*","", hc$labels), levels = UMI.subset))]
pagoda.view.aspects(tamr2, cell.clustering = hc, box = T, labCol = NA, margins = c(0.5, 100), col.cols = l2cols, cexRow=5,xpd=T, row.clustering = F, show.row.var.colors =F)
legend(x = 1500, y = 80, legend = paste(UMI.subset, " (",nb.cells[UMI.subset],")",sep=""), fill =  c("coral", "lightblue", "bisque"), box.col = "white", cex = 6, xpd = T)
dev.off()

for (i in go.id){
  png(paste("GO-", Term(GOTERM)[i],"-HypMusNeur-1500.png"), height = 20, width = 20, res = 300, units = "in")
  pagoda.show.pathways(i, varinfo, gl.env, cell.clustering = hc, margins = c(1,1), show.cell.dendrogram = T, showRowLabels = F, showPC = F, colcols = l2cols, cexRow = 6, n.genes = Inf)
  dev.off()
}

i<- "GO:0007186" 
png(paste("GO-", Term(GOTERM)[i],"-test.png"), height = 40, width = 20, res = 300, units = "in")
pagoda.show.pathways(i, varinfo, gl.env, cell.clustering = hc, margins = c(1,3), show.cell.dendrogram = T, showRowLabels = T, showPC = F, colcols = l2cols, cexRow = 1, n.genes = Inf)
dev.off()

# alternative clustering 
cor.mat <- matrix(0, nrow = ncol(tamr2$xv), ncol = ncol(tamr2$xv))
for (i in 1:nrow(cor.mat)){
  for (j in i:nrow(cor.mat)){
    cor.mat[i,j] <- cor(tamr2$xv[,i], tamr2$xv[,j])
  }
}
d <- as.dist(1-cor.mat)
d <- dist(t(tamr2$xv))
hc <- hclust(d)
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label") 
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col= c("coral", "lightblue", "bisque")[as.integer(factor(gsub("-\\d*","", label), levels = UMI.subset))])
  }
  return(x)
}
## apply labelCol on all nodes of the dendrogram
d <- dendrapply(as.dendrogram(hc), labelCol)

pdf("test.pdf",   height = 20, width=100)
plot(d)
dev.off()
pagoda.view.aspects(tamr2, cell.clustering = hc, box = T, labCol = NA, margins = c(0.5, 10), col.cols = l2cols, cexRow=5,xpd=T, row.clustering = F, show.row.var.colors =F)


################### Keggg
load("1205-fullResults-kegg.RData")

library(scde)
library(KEGGREST)

sel.kegg.id <- rownames(tamr2$xv)
kegg.id <- gsub("#PC1# ", "", rownames(tamr2$xv))
rownames(tamr2$xv) <- sapply(kegg.id, function(x){ 
  name <- keggGet(x)[[1]]$NAME;
  name <- gsub(" - Caenorhabditis elegans \\(nematode\\)","",name)
  go.size <- tamr2$df[tamr2$df$name %in% x,]$n
  name <- paste(name, "(", go.size,")", sep="")
  return(gsub(" ","\n",name))})


nb.cells <- sapply(unique(gsub("-\\d*","",colnames(tamr2$xv))), function(x){sum(gsub("-\\d*","",colnames(tamr2$xv))==x)})

png("KEGG-HypMusNeur-1500-2.png", height = 20, width = 40, units = "in", res=300)
par(xpd=T)
l2cols <- c("coral", "lightblue", "bisque")[as.integer(factor(gsub("-\\d*","", hc$labels), levels = UMI.subset))]
pagoda.view.aspects(tamr2, cell.clustering = hc, box = T, labCol = NA, margins = c(0.5, 100), col.cols = l2cols, cexRow=5,xpd=T)
legend(x = 1500, y = 400, legend = paste(UMI.subset, " (",nb.cells[UMI.subset],")",sep=""), fill =  c("coral", "lightblue", "bisque"), box.col = "white", cex = 6, xpd = T)
dev.off()

for (i in kegg.id){
  png(paste("KEGG-", gsub(" - Caenorhabditis elegans \\(nematode\\)","",keggGet(i)[[1]]$NAME),"-HypMusNeur-1500.png"), height = 20, width = 20, units = "in", res = 300)
  pagoda.show.pathways(i, varinfo, gl.env, cell.clustering = hc, margins = c(1,1), show.cell.dendrogram = T, showRowLabels = F, showPC = F, colcols = l2cols, cexRow = 6, n.genes = Inf)
  dev.off()
}




####################
## Summarizing plot 

x <- t(tamr2$xv)
#x is the matrix to biplot, x is numeric, thus variables on ratio or interval scales
#x has dimnames(x)[[2]] defined (and eventually dimnames(x)[[1]] defined)
# PCA 
xm<-apply(x,2,mean)
y<-sweep(x,2,xm)
ss<-(t(y)%*%y)
s<-ss/(nrow(x)-1)
d<-(diag(ss))^(-1/2)
e<-diag(d,nrow=ncol(x),ncol=ncol(x))
z<-y%*%e
r<-t(z)%*%z
q<-svd(z)
gfd<-((q$d[1])+(q$d[2]))/sum(q$d)
gfz<-(((q$d[1])^2)+((q$d[2])^2))/sum((q$d)^2)
gfr<-(((q$d[1])^4)+((q$d[2])^4))/sum((q$d)^4)
l<-diag(q$d,nrow=ncol(x),ncol=ncol(x))
R.B<-q$u        #scores matrix
C.B<-q$v%*%l    #loadings
#possibility to stretch scores by a scale factor
#scalefactor<-3.5
#R.B<-q$u *scalefactor
library(plotrix)

png("GO-Biplot-1500.png", height = 20, width = 20, units = "in", res=300)
l2cols <- c("coral", "lightblue", "bisque")[as.integer(factor(gsub("-\\d*","", colnames(tamr2$xv)), levels = UMI.subset))]
plot(R.B[ ,1]/max(R.B[,1]),R.B[ ,2]/max(R.B[,2]),axes=F,xlim=c(-1.2,1.2),ylim=c(-1.1,1.1),xlab=' ',ylab=' ',cex=2, pch = 16, col = addTrans(l2cols, 100))
legend(-1.15,1.15, legend = paste(UMI.subset, " (",nb.cells[UMI.subset],")",sep=""), fill =  c("coral", "lightblue", "bisque"), box.col = "white", cex = 3, xpd = T)
draw.circle(0,0,1,border='grey', lwd=3, lty=2)
for (i in seq(1,nrow(C.B),by=1)){
  arrows(0,0,C.B[i,1],C.B[i,2], length = 0.75, angle = 10, lwd = 2)
}
text(C.B[,1]-.04,C.B[,2]+.04,as.character(dimnames(x)[[2]]),cex=3, xpd = T)
dev.off()




addTrans <- function(color,trans){
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}

