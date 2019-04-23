library(biomaRt, quietly=T)
library(GO.db, quietly=T)

# ----- load ensembl data -----
dir.create("output")
db <- useMart("ensembl")
mm <- useDataset("mmusculus_gene_ensembl", mart = db)

# ----- get go2gene -----
res.go2gene <- getBM(attributes = c('go_id', 'mgi_symbol'), mart = mm)
res.go2gene <- res.go2gene[order(res.go2gene$go_id),]
res.go2gene <- subset(res.go2gene, subset=go_id != '')
rownames(res.go2gene) <- 1:nrow(res.go2gene)
write.table(res.go2gene, file='output/go2gene.txt', quote=F, row.names=F, sep='\t')

# ----- get go2term -----
res.go2term <- getBM(attributes = c('go_id', 'name_1006'), mart = mm)
res.go2term <- res.go2term[order(res.go2term$go_id),]
res.go2term <- subset(res.go2term, subset=go_id != '')
res.go2term <- unique(res.go2term)
rownames(res.go2term) <- 1:nrow(res.go2term)
write.table(res.go2term, file='output/go2term.txt', quote=F, row.names=F, sep='\t')

# ----- get go2offsprings -----
go2offsprings <- c()
for(i in 1:nrow(res.go2term)){
  print(paste(i, '/', nrow(res.go2term), ' GOs have been processed.'))
  go <- res.go2term[i,]$go_id
  offsprings <- unique(c(GOBPOFFSPRING[[go]], GOMFOFFSPRING[[go]], GOCCOFFSPRING[[go]]))
  go2offsprings <- c(go2offsprings, go, paste(offsprings, collapse=" "))
}
go2offsprings.dt <- as.data.frame(matrix(go2offsprings, ncol=2, byrow=T))
colnames(go2offsprings.dt) <- c('go_id', 'offsprings')
write.table(go2offsprings.dt, file="output/go2offsprings.txt", quote=F, row.names=F, sep="\t")

