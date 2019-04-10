# gastrula-reconst
Novel computational modeling and analysis of mouse mid-gastrula morphogenesis by self-organizing-map (SOM) clustering. 

## Requirement
* R version	3.5.1
* BiocInstaller	1.30.0
* GO.db			3.6.0
* biomaRt 		2.36.1 (Ensembl Release 94)
* genefilter		1.62.0
* goProfiles		1.42.0
* topGO				2.32.0
* org.Mm.eg.db			3.6.0
* corrplot			0.84
* gplots			3.0.1
* RColorBrewer			1.1.2
* getopt			1.20.2
* rgl					0.99.16

## Expression data
The input expression data ***GSE65924_E1.gene.expression.txt.gz*** can be downloaded 
from NCBI Gene Expression Omnibus [GSE65924](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE65924). 
The expression data was obtained from four sites (anterior, posterior, left, and right) of eleven frozen sections cut out 
from a mid-gastrula mouse embryo by [Peng et al.](https://www.ncbi.nlm.nih.gov/pubmed/27003939), 
where it is not single-cell data but it is composed of a small number of cells 
(about 20 cells per sample). The expression levels of 23,361 genes for 41 samples are stored as values normalized 
to FPKM (fragments per kilobase of exon per million reads mapped fragments).

## Usage

### 1. Get Gene Ontology Information
***get_GO.r*** creates ***output*** directory, 
and generates ***go2term.txt***, ***go2offsprings.txt***, 
and ***go2gene.txt*** in the directory using ***biomaRt*** package.
```
$ Rscript get_GO.r
```
**[Output]**
- output/go2term.txt
- output/go2offsprings.txt
- output/go2gene.txt


### 2. Generate sub gene expression tables
***generate_exprs.r*** generates sub gene expression tables according to PCA, GO, random variables, and DEG analysis.
```
$ Rscript generate_exprs.r
```
**[Output]**
- exprs_go/exprs.log10.E1.GOXXXXXX.txt
- exprs_pca_go/exprs.log10.E1.GOXXXXXX.txt
- exprs_peng/exprs.log10.E1.allgenes.txt
- exprs_peng/exprs.log10.E1.deg_peng.txt
- exprs_random/exprs.log10.E1.random.X.txt
- exprs_topGO/exprs.log10.E1.GOXXXXXX.txt
- output/topGO_KruskalWallisTest.txt


### 3. Count the number of genes
***count_genes.sh*** counts the number of genes for each expression table in 
***exprs_go***, ***exprs_pca_go***, ***exprs_peng***, ***exprs_topGO***, and ***exprs_random*** directories.
```
$ ./count_genes.sh
``` 
**[Output]**
- output/go2size.go.txt
- output/go2size.pca_go.txt
- output/go2size.peng.txt
- output/go2size.topGO.txt
- output/go2size.random.txt


### 4. Run Self-organizing map (SOM) algorithm (*submit jobs via `qsub`*)
***run_som.sh*** executes SOM clustering for all gene expression data in 
***exprs_go***, ***exprs_pca_go***, and other input directories. 
It will take several hours to finish SOM clustering for all expression data.
Please modify ***run_som.sh*** according to the environment of user's cluster machine if necessary.
```
$ ./run_som.sh exprs_go tiny
$ ./run_som.sh exprs_pca_go tiny
$ ./run_som.sh exprs_peng tiny
$ ./run_som.sh exprs_random tiny
$ ./run_som.sh exprs_topGO tiny
```
**[Output]**
- result_som.exprs_go
- result_som.exprs_pca_go
- result_som.exprs_peng
- result_som.exprs_random
- result_som.exprs_topGO


### 5. Evaluation (*submit jobs via `qsub`*)
***run_compScore.sh*** computes ***success rate*** and ***variance*** for all SOM results in \<input directory\>.
Please modify ***run_compScore.sh*** according to the environment of user's cluster machine if necessary.
```
$ mkdir result_score_tables
$ ./run_compScore.sh result_som.exprs_go result_score_tables/score_table_go.txt
$ ./run_compScore.sh result_som.exprs_pca_go result_score_tables/score_table_pca_go.txt
$ ./run_compScore.sh result_som.exprs_random result_score_tables/score_table_random.txt
$ ./run_compScore.sh result_som.exprs_peng result_score_tables/score_table_peng.txt
$ ./run_compScore.sh result_som.exprs_topGO result_score_tables/score_table_topGO.txt
```
**[Output]**
- result_score_tables/score_table_go.txt
- result_score_tables/score_table_pca_go.txt
- result_score_tables/score_table_random.txt
- result_score_tables/score_table_peng.txt
- result_score_tables/score_table_topGO.txt


### 6. Generate gene expression tables consisting of two GOs and run SOM
***generate_comb_pair_exprs.r*** generates gene expression tables based on combinations of two GOs including GO:0060412.
After running ***run_som.sh*** and ***run_compScore.sh***, the user will get ***score_table_go_comb2.txt***.
```
### 2 GOs combinations including GO:0060412 ###
$ Rscript generate_comb_pair_exprs.r --help
$ Rscript generate_comb_pair_exprs.r -t result_score_table/score_table_go.txt -g GO0060412 -k 1 -i exprs_go -o exprs_go_comb2
$ ./run_som.sh exprs_go_comb2 tiny
$ ./run_compScore.sh result_som.exprs_go_comb2 result_score_tables/score_table_go_comb2.txt
```
**[Output]**
- exprs_go_comb2/exprs.log10.E1.GOXXXXXX-GOXXXXXX.txt
- result_som.exprs_go_comb2
- result_score_tables/score_table_go_comb2.txt


### 7. Generate gene expression tables consisting of more than two GOs and run SOM
***generate_comb_multi_exprs.r*** generates gene expression tables based on combinations of three or more GOs including GO:0060412 with high success_rate.
After running ***run_som.sh*** and ***run_compScore.sh***, the user will get ***score_table_go_comb3.txt***, ***score_table_go_comb4.txt***, ***score_table_go_comb5.txt***, and ***score_table_go_comb6.txt***.
```
### 3 GOs combinations including GO:0060412 ###
$ Rscript generate_comb_multi_exprs.r --help
$ Rscript generate_comb_multi_exprs.r -t result_score_tables/score_table_go_comb2.txt -g GO0060412 -k 2 -c 0.85 -i exprs_go -o exprs_go_comb3
$ ./run_som.sh exprs_go_comb3 tiny
$ ./run_compScore.sh result_som.exprs_go_comb3 result_score_tables/score_table_go_comb3.txt

### 4 GOs combinations including GO:0060412 ###
$ Rscript generate_comb_multi_exprs.r -t result_score_tables/score_table_go_comb3.txt -g GO0060412 -k 3 -c 0.85 -i exprs_go -o exprs_go_comb4
$ ./run_som.sh exprs_go_comb4 tiny 
$ ./run_compScore.sh result_som.exprs_go_comb4 result_score_tables/score_table_go_comb4.txt

### 5 GOs combinations including GO:0060412 ###
$ Rscript generate_comb_multi_exprs.r -t result_score_tables/score_table_go_comb4.txt -g GO0060412 -k 4 -c 0.85 -i exprs_go -o exprs_go_comb5
$ ./run_som.sh exprs_go_comb5 tiny
$ ./run_compScore.sh result_som.exprs_go_comb5 result_score_tables/score_table_go_comb5.txt

### 6 GOs combinations including GO:0060412 ###
$ Rscript generate_comb_multi_exprs.r -t result_score_tables/score_table_go_comb5.txt -g GO0060412 -k 5 -c 0.85 -i exprs_go -o exprs_go_comb6
$ ./run_som.sh exprs_go_comb6 tiny
$ ./run_compScore.sh result_som.exprs_go_comb6 result_score_tables/score_table_go_comb6.txt
```
**[Output]**
- exprs_go_comb3/exprs.log10.E1.GOXXXXXX-GOXXXXXX-GOXXXXXX.txt
- exprs_go_comb4/exprs.log10.E1.GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX.txt
- exprs_go_comb5/exprs.log10.E1.GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX.txt
- exprs_go_comb6/exprs.log10.E1.GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX-GOXXXXXX.txt
- result_som.exprs_go_comb3
- result_som.exprs_go_comb4
- result_som.exprs_go_comb5
- result_som.exprs_go_comb6
- result_score_tables/score_table_go_comb3.txt
- result_score_tables/score_table_go_comb4.txt
- result_score_tables/score_table_go_comb5.txt
- result_score_tables/score_table_go_comb6.txt


### 8. Generate gene expression tables without one or more genes
***generate_comb_genes_exprs.r*** generates gene expression tables in which combinatorial genes were deleted. The number of genes to be deleted can be specified by `-k` option.
```
### Delete a gene ###
$ Rscript generate_comb_genes_exprs.r --help
$ Rscript generate_comb_genes_exprs.r -e exprs_go_comb5/exprs.log10.E1.GO0060412-GO0005021-GO2000392-GO0031994-GO0070986.txt -k 1 -o exprs_go_comb5_del1
$ ./run_som.sh exprs_go_comb5_del1 tiny
$ ./run_compScore.sh result_som.exprs_go_comb5_del1 result_score_tables/score_table_go_comb5_del1.txt

### Delete two genes ###
$ Rscript generate_comb_genes_exprs.r -e exprs_go_comb5/exprs.log10.E1.GO0060412-GO0005021-GO2000392-GO0031994-GO0070986.txt -k 2 -o exprs_go_comb5_del2
$ ./run_som.sh exprs_go_comb5_del2 tiny
$ ./run_compScore.sh result_som.exprs_go_comb5_del2 result_score_tables/score_table_go_comb5_del2.txt

### Delete three genes ###
$ Rscript generate_comb_genes_exprs.r -e exprs_go_comb5/exprs.log10.E1.GO0060412-GO0005021-GO2000392-GO0031994-GO0070986.txt -k 3 -o exprs_go_comb5_del3
$ ./run_som.sh exprs_go_comb5_del3 tiny
$ ./run_compScore.sh result_som.exprs_go_comb5_del3 result_score_tables/score_table_go_comb5_del3.txt
```
**[Output]**
- exprs_go_comb5_del1
- exprs_go_comb5_del2
- exprs_go_comb5_del3
- result_score_tables/score_table_go_comb5_del1.txt
- result_score_tables/score_table_go_comb5_del2.txt
- result_score_tables/score_table_go_comb5_del3.txt


### 9. Count genes of the expression tables consisting of multiple GOs
***count_genes_comb.sh*** counts the number of genes for each expression tables in 
***exprs_go_comb2***, ***exprs_go_comb3***, ***exprs_go_comb4***, ***exprs_go_comb5***, ***exprs_go_comb6***, 
***exprs_go_comb5_del1***, ***exprs_go_comb5_del12***, and ***exprs_go_comb5_del13***.
```
$ ./count_genes_comb.sh
``` 
**[Output]**
- output/go2size.go_comb2.txt
- output/go2size.go_comb3.txt
- output/go2size.go_comb4.txt
- output/go2size.go_comb5.txt
- output/go2size.go_comb6.txt
- output/go2size.go_comb5_del1.txt
- output/go2size.go_comb5_del2.txt
- output/go2size.go_comb6_del3.txt


### 10. Plot SOM results
***plot_som_results.r*** plots figures of "success rate vs variance", "the number of features vs success rate". 
```
$ Rscript plot_som_results.r --help 
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go.txt -o output_go
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.pca_go.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_pca_go.txt -o output_pca_go
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb2.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb2.txt -o output_go_comb2
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb3.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb3.txt -o output_go_comb3
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb4.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb4.txt -o output_go_comb4
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb5.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb5.txt -o output_go_comb5
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb6.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb6.txt -o output_go_comb6
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb5_del1.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb5_del1.txt -o output_go_comb5_del1
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb5_del2.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb5_del2.txt -o output_go_comb5_del2
$ Rscript plot_som_results.r -t output/go2term.txt -s output/go2size.go_comb5_del3.txt -r result_score_tables/score_table_random.txt -q result_score_tables/score_table_go_comb5_del3.txt -o output_go_comb5_del3
```
**[Output]**
- \<output directory\>/distributions.pdf
- \<output directory\>/go2term.sig.success_rate.txt	# listed GOs with FDR of success rate > 0.05
- \<output directory\>/go2term.sig.txt			    # listed GOs with FDR of success rate > 0.05 and FDR of variance > 0.05
- \<output directory\>/go2term.sig.variance.txt		      # listed GOs with FDR of variance > 0.05
- \<output directory\>/tbl.all.txt			      	         # listed all GOs
- \<output directory\>/tbl.sig.success_rate.txt				   # listed GOs with FDR of success rate > 0.05
- \<output directory\>/tbl.sig.txt					     	      # listed GOs with FDR of success rate > 0.05 and FDR of variance > 0.05
- \<output directory\>/tbl.sig.variance.txt					      	  # listed GOs with FDR of variance > 0.05


### 11. Plot 3D-models
***plot_3D-model.r*** generates png images of reconstructed 3d-model from the SOM results with the specified success rate by option `-c` in 
***result_som.exprs_go***, ***result_som.exprs_pca_go***, ***result_som.exprs_go_comb5***, and ***result_som.exprs_go_comb5_del2***. 
The cutoff of the success rate can be change by the `-c` option. If you want to generate movies, set `True` to the `-m` option. 
The type of 3D-model can be specified by `-t` option (`-t 1`: random model, `-t 2`: similarity-based model). 
```
$ Rscript plot_3D-model.r --help
$ Rscript plot_3D-model.r -i result_som.exprs_go -e exprs_go -s result_score_tables/score_table_go.txt -c 0.6 -m False -t 1 -o output_go/3d-model_exprs_go
$ Rscript plot_3D-model.r -i result_som.exprs_pca_go -e exprs_pca_go -s result_score_tables/score_table_pca_go.txt -c 0.6 -m False -t 1 -o output_pca_go/3d-model_exprs_pca_go
$ Rscript plot_3D-model.r -i result_som.exprs_go_comb5 -e exprs_go_comb5 -s result_score_tables/score_table_go_comb5.txt -c 0.95 -m False -t 1 -o output_go_comb5/3d-model_exprs_go_comb5
$ Rscript plot_3D-model.r -i result_som.exprs_go_comb5_del2 -e exprs_go_comb5_del2 -s result_score_tables/score_table_go_comb5_del2.txt -c 0.95 -m False -t 1 -o output_go_comb5_del2/3d-model_exprs_go_comb5_del2
```
**[Output]**
- output_go/3d-model_exprs_go
- output_pca_go/3d-model_exprs_pca_go
- output_go/3d-model_exprs_go_comb5
- output_go/3d-model_exprs_go_comb5_del2


### 12. Plot heatmaps
***plot_heatmap.r*** generates heat maps for all gene expression tables listed 
in ***\<directory name\>/tbl.all.txt***. 
The cutoff of the success rate can be change by the `-c` option.
```
$ Rscript plot_heatmap.r --help
$ Rscript plot_heatmap.r -e exprs_go -i output_go/tbl.all.txt -c 0.6 -o output_go/heatmap.pdf
$ Rscript plot_heatmap.r -e exprs_pca_go -i output_pca_go/tbl.all.txt -c 0.6 -o output_pca_go/heatmap.pdf
$ Rscript plot_heatmap.r -e exprs_go_comb5 -i output_go_comb5/tbl.all.txt -c 0.95 -o output_go_comb5/heatmap.pdf
$ Rscript plot_heatmap.r -e exprs_go_comb5_del2 -i output_go_comb5_del2/tbl.all.txt -c 0.95 -o output_go_comb5_del2/heatmap.pdf
```
**[Output]**  
- \<output directory\>/heatmap.pdf


### 13. Plot correlation matrix of domains
***plot_corMatrix_domain.r*** plots correlation matrix of domains as heat maps 
for all GOs listed in ***\<directory name\>/tbl.all.txt***.
The cutoff of the success rate can be change by the `-c` option.
```
$ Rscript plot_corMatrix_domain.r --help
$ Rscript plot_corMatrix_domain.r -e exprs_go -i output_go/tbl.all.txt -t data/sample2domain.txt -c 0.6 -o output_go/corr_domain.pdf
$ Rscript plot_corMatrix_domain.r -e exprs_pca_go -i output_pca_go/tbl.all.txt -t data/sample2domain.txt -c 0.6 -o output_pca_go/corr_domain.pdf
$ Rscript plot_corMatrix_domain.r -e exprs_go_comb5 -i output_go_comb5/tbl.all.txt -t data/sample2domain.txt -c 0.95 -o output_go_comb5/corr_domain.pdf
$ Rscript plot_corMatrix_domain.r -e exprs_go_comb5_del2 -i output_go_comb5_del2/tbl.all.txt -t data/sample2domain.txt -c 0.95 -o output_go_comb5_del2/corr_domain.pdf
```
**[Output]**
- \<output directory\>/corr_domain.pdf


### 14. Plot correlation matrix of cells
***plot_corMatrix_cell.r*** plots correlation matrix of cells as heat maps 
for all GOs listed in ***\<directory name\>/tbl.all.txt***.
The cutoff of the success rate can be change by the `-c` option.
```
$ Rscript plot_corMatrix_cell.r --help
$ Rscript plot_corMatrix_cell.r -e exprs_go -i output_go/tbl.all.txt -t data/sample2domain.txt -c 0.6 -o output_go/corr_cell.pdf
$ Rscript plot_corMatrix_cell.r -e exprs_pca_go -i output_pca_go/tbl.all.txt -t data/sample2domain.txt -c 0.6 -o output_pca_go/corr_cell.pdf
$ Rscript plot_corMatrix_cell.r -e exprs_go_comb5 -i output_go_comb5/tbl.all.txt -t data/sample2domain.txt -c 0.95 -o output_go_comb5/corr_cell.pdf
$ Rscript plot_corMatrix_cell.r -e exprs_go_comb5_del2 -i output_go_comb5_del2/tbl.all.txt -t data/sample2domain.txt -c 0.95 -o output_go_comb5_del2/corr_cell.pdf
```
**[Output]**
- \<output directory\>/corr_cell.pdf