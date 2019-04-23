library(genefilter)
library(GO.db)
library(org.Mm.eg.db)
dir.create('exprs_go',showWarnings=F)
dir.create('exprs_peng', showWarnings=F)

#===========================
# class
#===========================
setClass(
  "gastrula", 
  representation(
    exprs = "matrix", 
    exprs.filtered = "matrix", 
    factors = "factor", 
    top_genes = "vector"
  ), 
  prototype = list(
    expr = NULL, 
    exprs.filtered = NULL, 
    factors = NULL, 
    top_genes = NULL
  )
)

#===============================
# functions
#===============================
filter_genes <- function(na.rm=T){
  function(x){
    (length(x[x > 1]) >= 2) && (var(log10(x+1)) > 0.05)
  }
}

#==============================
# generate_instances
#==============================
E1 <- new("gastrula")
E1@exprs <- as.matrix(read.table('data/GSE65924_E1.gene.expression.txt', header=T, row.names=1, sep='\t'))
E1@factors <- as.factor(
  c('D1', 'D4', 
    'D1', 'D4', 'D2', 'D2', 
    'D4', 'D2', 'D2', 
    'D1', 'D4', 'D2', 'D2', #  9~12
    'D1', 'D4', 'D2', 'D2', # 13~16
    'D1', 'D4', 'D2', 'D2', # 17~20
    'D1', 'D4', 'D3', 'D2', # 21~24
    'D1', 'D4', 'D3', 'D3', # 25~28
    'D1', 'D4', 'D3', 'D3', # 29~32
    'D1', 'D4', 'D3', 'D3', # 33~36
    'D3', 'D4', 'D3', 'D3') # 37~40
)

#==============================
# main
#==============================
#----- gene filtering -----
E1.good <- genefilter(E1@exprs, filterfun(filter_genes()))
E1@exprs.filtered <- E1@exprs[E1.good,]

#----- get expression table of all genes -----
all.exprs <- E1@exprs
ofpath <- 'exprs_peng/exprs.log10.E1.allgenes.txt'
write.table(log10(all.exprs + 1), ofpath, quote=F, sep='\t')

#----- get expression table of pca genes selected by peng et al. -----
peng158 <- sort(rownames(read.table('data/peng158.txt', header=F, row.names=1, sep='\t')))
peng158.exprs <- E1@exprs[match(peng158, rownames(all.exprs)), ]
ofpath <- 'exprs_peng/exprs.log10.E1.peng158.txt'
write.table(log10(peng158.exprs + 1), ofpath, quote=F, sep='\t')

#----- get expression table according to GOs -----
GO2gene <- read.table('output/go2gene.txt', header=T, sep='\t')
all_GOs <- read.table('output/go2term.txt', header=T, quote="", sep='\t')
go.db.BP <- as.list(GOBPOFFSPRING)
go.db.MF <- as.list(GOMFOFFSPRING)
go.db.CC <- as.list(GOCCOFFSPRING)
cnt <- 0
for(GO in all_GOs$go_id){
  GOs <- unique(c(GO, go.db.BP[[GO]], go.db.MF[[GO]], go.db.CC[[GO]]))
  genes <- subset(GO2gene, subset=go_id %in% GOs)$mgi_symbol
  genes <- genes[!is.na(genes)]
  genes_filtered <- intersect(rownames(E1@exprs.filtered), genes)
  feature_genes <- sort(unique(genes_filtered))
  ofile <- paste('exprs_go/exprs.log10.E1.', gsub(':', '', GO), '.txt', sep='')
  if(length(feature_genes) >= 3  && length(genes) <= 1000){
    print(GO)
    write.table(log10(E1@exprs.filtered + 1)[feature_genes, , drop=F], ofile, quote=F, sep='\t')
  }
}

