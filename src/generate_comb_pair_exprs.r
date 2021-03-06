library(getopt, quietly=T)

#======================================
# Options
#======================================
spec <- matrix(c(
  "help",     "h", 0, "logical", "show this help", 
  "inptable", "t", 2, "character", "[required] input som result table", 
  "cutoff",   "c", 1, "numeric",   "[required] cutoff of the success rate", 
  "go"    ,   "g", 1, "character", "[required] specify a GO which is always included in combinations (e.g., GO0060412)", 
  "comb"  ,   "k", 2, "numeric",   "[required] specify k of n_C_k (choose k from n)", 
  "inpdir",   "i", 2, "character", "[required] input directory containing expression data", 
  "outdir",   "o", 2, "character", "[required] output directory"
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
inptable <- read.table(opt$inptable, header=T, sep="\t")
candidate_gos <- as.vector(subset(inptable, success_rate >= opt$cutoff)$sample)
combs <- combn(candidate_gos, opt$comb)
for(i in 1:ncol(combs)){
  gos <- combs[,i]
  exprs <- read.table(paste0(opt$inpdir, '/exprs.log10.E1.', opt$go, '.txt'), header=T, row.names=1, sep='\t')
  print(gos)
  for(go in gos){
    exprs.tmp <- read.table(paste0(opt$inpdir, '/exprs.log10.E1.', go, '.txt'), header=T, row.names=1, sep='\t')
    exprs <- rbind(exprs, exprs.tmp)
  }
  exprs <- unique(exprs)
  exprs <- exprs[order(rownames(exprs)),]
  ofname <- paste0('exprs.log10.E1.', opt$go, '-', go, '.txt')
  write.table(exprs, paste(opt$outdir, ofname, sep='/'), quote=F, sep='\t')
}

#----- create output directory -----
dir.create(opt$outdir, showWarnings=F)
odir <- gsub("\\/$", "", opt$outdir)


