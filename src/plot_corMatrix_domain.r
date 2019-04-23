library(corrplot, quietly=T)
library(RColorBrewer, quietly=T)
library(gplots, quietly=T)
library(getopt, quietly=T)

# ======================================
# Options
# ======================================
spec <- matrix(c(
  "help",      "h", 0, "logical", "show this help", 
  "exprdir",   "e", 2, "character", "[required] expression directory", 
  "inputfile", "i", 2, "character", "[required] input file", 
  "corresp",   "t", 2, "character", "[required] correspondence table between samples ans domains", 
  "outfile",   "o", 2, "character", "[required] output file", 
  "cutoff",    "c", 1, "numeric",   "cutoff of success rates for heatmap plottings. (default: 0.6)"
), ncol=5, byrow=T)

opt <- getopt(spec)
if(!is.null(opt$help)){
  cat(getopt(spec, usage=T))
  quit(status=1)
}

if(is.null(opt$cutoff)){
  opt$cutoff = 0.6
}


#===============================
# Functions
#===============================
plot_corr_matrix <- function(expr){
  heatcols <- c(grey(1:20/20), rich.colors(20))
  corrplot(cor(expr), method='color', tl.col='black', cl.lim=c(0,1), col=heatcols)
}

#----------------------------
sub.plot_heatmap <-function(expr, main_title){
  par(cex.main=min(1, 40/nchar(main_title)))
  heatmap.2(expr, scale='none', dendrogram='column', Rowv=T, Colv=T, density.info='none', key.title='', 
            key.xlab = 'Correlation',  breaks=seq(0,1,0.05), 
            trace='none', main=main_title, col=rich.colors(20), cexCol=1.5, cexRow=1.5, keysize=1,
            hclustfun = hclust, offsetCol = 0.1, adjCol = c(NA,0.5))
}

#-----------------------------
plot_heatmap <- function(gos, idir, sample2domain, prefix){
  
  if(length(gos) > 0){
    for(go in gos){
      print(go)
      ifile <- paste(idir, 'exprs.log10.E1.', go, '.txt', sep='')
      expr <- read.table(ifile, row.names=1, header=T, sep='\t')
      
      #----- compute average expression level for each domain ------
      samples_d1 <- apply(expr[, rownames(subset(sample2domain, domain == 'D1'))], 1, mean)
      samples_d2 <- apply(expr[, rownames(subset(sample2domain, domain == 'D2'))], 1, mean)
      samples_d3 <- apply(expr[, rownames(subset(sample2domain, domain == 'D3'))], 1, mean)
      samples_d4 <- apply(expr[, rownames(subset(sample2domain, domain == 'D4'))], 1, mean)
      expr.domain <- matrix(c(samples_d1, samples_d2, samples_d3, samples_d4), ncol=4, byrow=F)
      colnames(expr.domain) <- c('D1', 'D2', 'D3', 'D4')
      rownames(expr.domain) <- rownames(expr)
      
      #----- plot figures -----
      main_title <- gsub("GO", "GO:", go)
      main_title <- gsub("comb_", "", go)
      sub.plot_heatmap(cor(expr.domain), main_title)
    }
  }else{
    plot(1:10, 1:10, type='n', axes=F, ann=F)
    text(5,5, 'No Data')
  }
}

#-----------------------------
plot_heatmap_allgenes <- function(ifile, sample2domain, title){
    expr <- read.table(ifile, row.names=1, header=T, sep='\t')
    
    #----- compute average expression level for each domain ------
    samples_d1 <- apply(expr[, rownames(subset(sample2domain, domain == 'D1'))], 1, mean)
    samples_d2 <- apply(expr[, rownames(subset(sample2domain, domain == 'D2'))], 1, mean)
    samples_d3 <- apply(expr[, rownames(subset(sample2domain, domain == 'D3'))], 1, mean)
    samples_d4 <- apply(expr[, rownames(subset(sample2domain, domain == 'D4'))], 1, mean)
    expr.domain <- matrix(c(samples_d1, samples_d2, samples_d3, samples_d4), ncol=4, byrow=F)
    colnames(expr.domain) <- c('D1', 'D2', 'D3', 'D4')
    rownames(expr.domain) <- rownames(expr)
    
    #----- plot figures -----
    sub.plot_heatmap(cor(expr.domain), title)
}

#===============================
# Main
#===============================
idir <- paste0(opt$exprdir, '/')
ifile <- opt$inputfile
ofile <- opt$outfile
cfile <- opt$corresp

tbl <- read.table(ifile, header=T, row.names=1, sep='\t')
tbl <- tbl[order(tbl$var),]
tbl <- tbl[order(tbl$success_rate, decreasing=T),]
tbl <- subset(tbl, success_rate >= as.numeric(opt$cutoff))
gos <- rownames(tbl)

sample2domain <- read.table(cfile, row.names=1, header=T, sep='\t')
pdf(file=ofile)
plot_heatmap(gos, idir, sample2domain, '')
dev.off()

