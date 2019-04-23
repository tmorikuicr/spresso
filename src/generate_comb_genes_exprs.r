library(getopt, quietly=T)

#======================================
# Options
#======================================
spec <- matrix(c(
  "help",     "h", 0, "logical", "show this help", 
  "inpexprs", "e", 1, "character", "[required] input gene expression table", 
  "ndel"  ,    "k", 2, "numeric",   "[required] specify the number of genes to be deleted", 
  "outdir",    "o", 2, "character", "[required] output directory"
), ncol=5, byrow=T)

opt <- getopt(spec)
if(!is.null(opt$help)){
  cat(getopt(spec, usage=T))
  quit(status=1)
}

if(is.null(opt$cutoff)){
  opt$cutoff = 0
}

#============================================
# Main
#============================================
dir.create(opt$outdir, showWarnings=F)
inpexprs <- read.table(opt$inpexprs, header=T, row.names=1, sep="\t")
candidate_genes <- as.vector(rownames(inpexprs))
combs <- combn(candidate_genes, length(candidate_genes) - opt$ndel)

for(i in 1:ncol(combs)){
  genes <- combs[,i]
  exprs <- read.table(opt$inpexprs, header=T, row.names=1, sep='\t')
  del_genes <- setdiff(candidate_genes, genes)
  exprs_del <- exprs[match(genes, rownames(exprs)), ]
  exprs_del <- exprs_del[order(rownames(exprs_del)),]
  comb_name <- paste(del_genes, collapse="-")
  ofname <- paste0('exprs.log10.E1.comb_del-', comb_name, '.txt')
  print(ofname)
  write.table(exprs_del, paste(opt$outdir, ofname, sep='/'), quote=F, sep='\t')
}

