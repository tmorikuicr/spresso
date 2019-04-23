library(rgl, quietly=T)
library(magick, quietly=T)
library(getopt, quietly=T)
pink <- '#ff5555'
orange <- '#ffae42'
blue2 <- '#5fd3bc'
purple <- '#aa87de'

domain_labels   <-   c('D1', 'D4', 
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

#======================================
# Options
#======================================
spec <- matrix(c(
  "help",     "h", 0, "logical", "show this help", 
  "inputdir", "i", 2, "character", "[required] input directory", 
  "exprdir",  "e", 2, "character", "[required] expression data directory", 
  "score",    "s", 2, "character", "[required] input score table", 
  "outdir",   "o", 2, "character", "[required] output directory", 
  "movie",    "m", 1, "logical", "generate movies (default: False)", 
  "cutoff",   "c", 1, "numeric", "cutoff of success rates for 3D plottings. (default: 0.6)"
), ncol=5, byrow=T)

opt <- getopt(spec)
if(!is.null(opt$help)){
  cat(getopt(spec, usage=T))
  quit(status=1)
}

if(is.null(opt$movie)){
  opt$movie = F
}

if(is.null(opt$cutoff)){
  opt$cutoff = 0.6
}

#================================================================
# functions
#================================================================
generate_colorVector_domain <- function(domain_labels){
  color_set <- c(pink, orange, blue2, purple)
  #color_set <- c('gray0', 'gray50', 'gray75', 'gray100')
  color_vec <- c()
  for(dom in domain_labels){
    if(dom == 'D1'){
      color_vec <- c(color_vec, color_set[1])
    }else if(dom == 'D4'){
      color_vec <- c(color_vec, color_set[2])
    }else if(dom == 'D2'){
      color_vec <- c(color_vec, color_set[3])
    }else if (dom == 'D3'){
      color_vec <- c(color_vec, color_set[4])
    }else{
      color_vec <- c(color_vec, 'black')
    }
  }
  return(color_vec)
}
color_vec <- generate_colorVector_domain(domain_labels)

#-----------------------------------------------------------------
corn_plot_3D_projection <- function(idx2class, map_vector, exprs){
  set.seed(1) #default
  pos.x <- c()
  pos.y <- c()
  pos.z <- c()
  
  n.cell <- nrow(idx2class)
  for(i in 1:n.cell){
    #-------------------------------------------------------------------
    # get x, y, z coordinates on cubic SOM (e.g.. 000 indicates x=0, y=0, z=0) 
    # * coord.c indicates the coordinates of centroid unit for cell i
    # * coord.c.adj_x (resp. ...adj_y, adj_z) indicates the coordinates 
    #   of the adjacents unit of the centroid unit
    #-------------------------------------------------------------------
    coord.c <- paste0(as.character(idx2class$x[i]), 
                        as.character(idx2class$y[i]), 
                        as.character(idx2class$z[i]))
    coord.c.adj_x <- paste0(as.character(as.integer(!idx2class$x[i])), 
                              as.character(idx2class$y[i]), 
                              as.character(idx2class$z[i]))
    coord.c.adj_y <- paste0(as.character(idx2class$x[i]), 
                              as.character(as.integer(!idx2class$y[i])), 
                              as.character(idx2class$z[i]))
    coord.c.adj_z <- paste0(as.character(idx2class$x[i]), 
                              as.character(idx2class$y[i]), 
                              as.character(as.integer(!idx2class$z[i])))
    
    #--------------------------------------------------------------------
    # get index of centroid and adjacent units (e.g., 000 -> 1, 001 -> 2, and so on)
    #--------------------------------------------------------------------
    index.c <- match(coord.c, sub("X", "", colnames(map_vector)))
    index.c.adj_x <- match(coord.c.adj_x, sub("X", "", colnames(map_vector)))
    index.c.adj_y <- match(coord.c.adj_y, sub("X", "", colnames(map_vector)))
    index.c.adj_z <- match(coord.c.adj_z, sub("X", "", colnames(map_vector)))
    
    #--------------------------------------------------------------------
    # get direction vectors from centroid unit to adjacent units
    #--------------------------------------------------------------------
    vec_cx <- map_vector[, index.c.adj_x] - map_vector[, index.c]
    vec_cy <- map_vector[, index.c.adj_y] - map_vector[, index.c]
    vec_cz <- map_vector[, index.c.adj_z] - map_vector[, index.c]
    
    #--------------------------------------------------------------------
    # get distance from centroid unit to adjacent units
    #--------------------------------------------------------------------
    dist_cx <- sqrt(sum(vec_cx * vec_cx))
    dist_cy <- sqrt(sum(vec_cy * vec_cy))    
    dist_cz <- sqrt(sum(vec_cz * vec_cz))
    
    #----------------------------------------------------------------------
    # get each distance from centroid unit to the projected point of cell i
    #----------------------------------------------------------------------
    dist_cp.x <- (sum(vec_cx * (exprs[, i] - map_vector[, index.c]))) / dist_cx
    dist_cp.y <- (sum(vec_cy * (exprs[, i] - map_vector[, index.c]))) / dist_cy
    dist_cp.z <- (sum(vec_cz * (exprs[, i] - map_vector[, index.c]))) / dist_cz
    
    #----------------------------------------------------------------------
    # get each scaled distance from centroid unit to the projected point of cell i
    #----------------------------------------------------------------------
    dist_cp.x.scaled <- dist_cp.x / dist_cx
    dist_cp.y.scaled <- dist_cp.y / dist_cy
    dist_cp.z.scaled <- dist_cp.z / dist_cz
    
    #----------------------------------------------------------------------
    # get the position of cell i
    #----------------------------------------------------------------------
    x <- (idx2class$x[i])
    y <- (idx2class$y[i])
    z <- (idx2class$z[i])
    if(x==0) {x <- x + dist_cp.x.scaled - 0.5} else {x <- x - dist_cp.x.scaled - 0.5}
    if(y==0) {y <- y + dist_cp.y.scaled - 0.5} else {y <- y - dist_cp.y.scaled - 0.5}
    if(z==0) {z <- z + dist_cp.z.scaled - 0.5} else {z <- z - dist_cp.z.scaled - 0.5}
    
    #----------------------------------------------------------------------
    # transform the cube structure to a embryo model
    #----------------------------------------------------------------------
    r <- sqrt(abs(z+1))
    theta <- atan2(y, x)
    pos.x <- c(pos.x, r*cos(theta))
    pos.y <- c(pos.y, r*sin(theta))
    pos.z <- c(pos.z, z)
  }
  
  open3d(windowRect=c(0, 0, 400, 400))
  plot3d(pos.x, pos.y, pos.z, 
         pch=21, col=color_vec, 
         type='s', size = 4,
         xlab='X', ylab='Y', zlab='Z',  lit=T)
  
  theta <- seq(0,2,pi/1000)*pi
  for(i in seq(-1,1,0.5)){
    r <- abs(polyroot(c(i+1, 0, 1))[1])
    lines3d(r*cos(theta), r*sin(theta), i, color='dimgrey')
  }
  r <- sqrt(abs(1+1))
  x <- seq(-r, r, 0.001)
  lines3d(x, 0, x^2-1, color='dimgrey')
  lines3d(0, x, x^2-1, color='dimgrey')
  lines3d(x, 0, x^2-1, color='dimgrey')
}

#-------------------------------------------------
make_example_list <- function(ipath, idirs){
  example_list <- c()
  for(idir in idirs){
    ifile <- paste(ipath, idir, 'scores.txt', sep='/')
    tbl <- read.table(ifile, header=F, '\t')
    colnames(tbl) <- c('method', 'seed', 'radii', 'itr', 'score')
    subtbl <- subset(tbl, subset = score == 1.0)
    if(nrow(subtbl) > 0){
      for(i in 1:nrow(subtbl)){
        example_list <- c(example_list, idir, subtbl[i,]$seed)
      }
    }
  }
  example_list <- as.data.frame(matrix(example_list, ncol=2, byrow=T))
  colnames(example_list) <- c('dirname', 'seed')
  return(example_list)
}

#-------------------------------------------------
run_corn_plot_3D <- function(ipath, epath, opath, example_list, movie){
  for(i in 1:nrow(example_list)){
    dirname <- example_list[i,1]
    seed <- example_list[i,2]
    print(dirname)
    ifile <- paste(ipath, '/', dirname, '/idx2class/idx2class.', seed, '.txt', sep='')
    idx2class <- read.table(ifile, row.names=1, header=T, sep='\t')
    mfile <- paste(ipath, '/', dirname, '/map_vector/map_vector.', seed, '.txt', sep='')
    map_vector <- read.table(mfile, row.names=1, header=T, sep='\t')
    #corn_plot_3D(idx2class)
    efile <- paste(epath, '/exprs.log10.E1.', dirname, '.txt', sep='')
    exprs <- read.table(efile, row.names=1, header=T, sep='\t')
    corn_plot_3D_projection(idx2class, map_vector, exprs)
    view_mat <- matrix(
                    c(0.8662382, -0.4996285, -0.001548417, 0, 
                      0.2102456, 0.3617015, 0.908277810, 0, 
                      -0.4532415, -0.7871108, 0.418364525, 0, 
                      0, 0, 0, 1), 
                    ncol=4, byrow=T)
    rgl.viewpoint(userMatrix = view_mat, zoom=0.9)
    go <- dirname
    ofile <- paste(opath, '/', go,'_3d_', seed, '.png', sep='')
    rgl.snapshot(ofile)
    if(movie == T){
      movie3d(spin3d(axis=c(0,0,1), rpm=5), duration=10, dir=opath, movie=paste0(go, '_3d_', seed, '_movie'))
    }
    rgl.close()
  }
}

#--------------------------------------------------
run_corn_plot_3D_single <- function(ifile){
  idx2class <- read.table(ifile, row.names=1, header=T, sep='\t')
  corn_plot_3D(idx2class)
}

#=================================================================
# main
#=================================================================
movie <- opt$movie
ipath <- opt$inputdir
epath <- opt$exprdir
spath <- opt$score
opath <- opt$outdir
cutoff <- opt$cutoff

dir.create(opath, recursive=T, showWarnings=F)
idirs <- list.files(ipath)
scoretbl <- read.table(spath, header=T, sep="\t")
topsamples <- as.vector(subset(scoretbl, success_rate >= cutoff)$sample)
idirs <- intersect(idirs, topsamples)

if(length(idirs) > 0){
  example_list <- make_example_list(ipath, idirs)
  print(example_list)
  run_corn_plot_3D(ipath, epath, opath, example_list, movie)
}else{
  print(paste0('No valid data in ', ipath))
}

#run_corn_plot_3D_single('result_som.exprs_pca_go/GO0044389/idx2class/idx2class.51.txt')
#rgl.snapshot('output/3d_model_pca_go_GO0044389_51.png')
#movie3d(spin3d(axis=c(0,0,1), rpm=5), duration=10, dir='output', movie='3d_movie_GO0044389_51')

