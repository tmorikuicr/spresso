library(getopt, quietly=T)

red <- '#ee204d'
blue <- '#0d98ba'

#======================================
# Options
#======================================
spec <- matrix(c(
  "help",   "h", 0, "logical",   "show this help", 
  "term",   "t", 2, "character", "[required] input file showing the correspondence between GOs and their terms", 
  "size",   "s", 2, "character", "[required] input file showing the correspondence between GOs and their size", 
  "score",  "q", 2, "character", "[required] input file for the SOM result",
  "outdir", "o", 2, "character", "[required] output directory"
), ncol=5, byrow=T)

opt <- getopt(spec)
if(!is.null(opt$help)){
  cat(getopt(spec, usage=T))
  quit(status=1)
}

#============================================
# Main
#============================================
#----- load data --------------------------
go2term <- read.table(opt$term, header=T, quote="", sep="\t")
go2term <- rbind(go2term, data.frame(go_id="GO:NULL", name_1006='NA'))
rownames(go2term) <- sub(":", "", go2term$go_id)

#----- load SOM results ---------------------
tbl <- read.table(opt$score, header=T, sep='\t')
go2size <- read.table(opt$size, header=T, sep='\t')
tbl <- cbind(tbl, go2size) 
tbl <- subset(tbl, select=c("sample","num_cell","size","success_rate","var"))
rownames(tbl) <- sub("exprs.log10.E1.", "", tbl$sample)

#----- output -----
dir.create(opt$outdir, showWarnings=F)
odir <- gsub("\\/$", "", opt$outdir)
write.table(tbl, paste(odir, 'tbl.all.txt', sep='/'), quote=F, sep='\t')

#----- plot success rate vs. variance -----
pdf(file=paste(odir, 'distributions.pdf', sep='/'))
x <- tbl$success_rate
y <- tbl$var
plot(x, y, xlab='Success rate', ylab = 'Total variance', pch=16, col='royalblue3', cex.lab=1.5)

#----- plot the number of genes vs. success rate -----
#x <- tbl$size
#y <- tbl$success_rate
#plot(x, y, xlab='# features', ylab = 'Success rate', pch=16, col='royalblue3', cex.lab=1.5)
dev.off()
