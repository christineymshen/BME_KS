library(doMC)
library(Matrix) # create block diagonal matrix
library(tidyverse)
library(abind)
library(wesanderson)
library(bde)
library(amen)
library(matrixStats)
library(rstan)

outfolder <- ""

pr_ge <- function(x,y,tol=1.5e-8){
  
  # probability of greater or equal to 0 (ge from the latex greater than or equal to)
  # x: matrix, m x n
  # y: matrix, m x n
  
  diff <- x-y
  diff[abs(diff)<tol] <- 0
  colMeans(diff >= 0)
  
}

combine <- function(labels,...){
  res <- list(...)
  nr <- length(res)
  nd <- length(res[[1]]); np <- length(res[[1]][[1]])
  
  for (i in 1:nd){
    for (j in 1:np){
      for (k in 2:nr){
        res[[1]][[i]][[j]] <- cbind(res[[1]][[i]][[j]], res[[k]][[i]][[j]])
      }
      colnames(res[[1]][[i]][[j]]) <- labels
    }
  }
  res[[1]]
}

# compile stan file, this is to be run in the terminal
cp_stan <- function(filename){
  
  op_filename <- paste0(str_remove(filename,".stan"),"_mod.rds")
  mod <- stan_model(filename)
  saveRDS(mod,op_filename)
  
}

# this version is used for runs in the 20240804 folder, with run time
agg <- function(runnum,folder,p_label){
  
  fdir <- paste0(outfolder,folder,"/")
  
  # Files to aggregate 
  f <- list.files(fdir)
  f <- f[stringr::str_detect(f, paste0(runnum,"_",p_label))]
  
  res <- readRDS(paste0(fdir,f[1]))
  res <- res$output
  
  nd <- dim(res)[2]
  labels <- dimnames(res)
  discrepancy <- labels[[2]]
  priors <- labels[[1]]
  
  # Read into memory
  output <- vector(mode="list", length=length(nd))
  runtime <- vector(mode="list",length=length(f))
  for (d in 1:nd){
    
    out <- NULL
    for (i in seq_along(f)){
      
      idx <- as.integer(str_extract(f[i],"(?<=_)\\d+(?=\\.rds)"))
      
      res <- readRDS(paste0(fdir,f[i]))
      runtime[[idx]] <- res$runtime
      out <- rbind(out, (res$output[,d]))
      
    }
    
    colnames(out) <- priors
    output[[d]] <- out
  }
  names(output) <- discrepancy
  res <- list(output=output,runtime=runtime)
  saveRDS(res, paste0("../res/",folder, "_",runnum,"_",p_label,".rds"))
  
}

wanted <- function(runnum,folder,p_label,nrun=1000){
  
  fdir <- paste0(outfolder,folder,"/")
  check <- logical(nrun)
  
  f <- list.files(fdir)
  f <- f[stringr::str_detect(f, paste0(runnum,"_",p_label,"_"))]
  for (i in seq_along(f)){
    idx <- as.integer(str_extract(f[i],"(?<=_)\\d+(?=\\.rds)"))
    check[idx] <- T
  }
  list(total=sum(check==F), list=paste(which(check==F), collapse = ","))
  
}

getrange <- function(res, op_density=1, na.rm=F){
  p <- dim(res)[2]
  upper <- numeric(p)
  
  for (i in 1:p){
    
    if (na.rm){
      res_i <- res[,i][!is.na(res[,i])]
    } else {
      res_i <- res[,i]
    }
    
    if (op_density==1){
      upper[i] <- max(density(res_i)$y)
    } else if (op_density==2) {
      upper[i] <- max(bde(res_i, estimator="boundarykernel")@densityCache )
    } else if (op_density==3) {
      upper[i] <- max(bd_density(res_i[[1]][,prior])$y)
    }
  }
  
  # if the largest is too large compared to the others, return the second largest
  if (max(max(upper) / upper) >3 && max(upper) > 5){
    max(upper[-which.max(upper)])
  } else {
    max(upper)
  }
}

# adapted from plotp_density3
# made changes so that ymax is calculated for each plot in a row before applying
# in the past there are cases where the two pictures are of different scales but only the scale for the plot on the left was showing up
plotp_density4 <- function(res, spec=NULL){
  
  # things can be specified in spec:
  # op_density: 1 for usual density; 2 (default) for bde; 3 for the density codes online
  # - d_idx: numeric vector, discrepancies to be plotted
  # - fs: font size, default at 0.8
  # - prior_idx: numeric vector, priors to be plotted
  # - p_idx: numeric vector, p-values to be plotted
  # - prior_label: character vector, labels for the priors, the order of this should follow the original order in the data rather than prior_idx
  # - d_label: character vector, labels for the discrepancy measures
  # - ymax: to control the plot ymax
  # - na.rm: if T, remove all na's. Default to be F so that I can get warnings for NAs
  
  # these labels are in the order I always want to be presenting (based on the structural comparison figure in the draft)
  p_labels_all <- c("ppost","pdpost","psplit","pdsplit","ptm","pepd","pptpost","pdptpost")
  p_labels_all_latex <- c(expression(p[post]),
                          expression(p[dpost]),
                          expression(p[split]),
                          expression(p[dsplit]),
                          expression(p[tm]),
                          expression(p[epd]),
                          expression(p[ptpost]),
                          expression(p[dptpost]))
  colors_all <- c(1,1,4,4,3,2,6,6)
  color_mapping <- data.frame(p_label=p_labels_all,color=colors_all)
  default_col_order <- c(1,4,5,6,3,2)
  
  # note cannot change to use ifelse when the TRUE/FALSE values are not same length
  if (is.null(spec$d_idx)){
    d_idx <- c(1:length(res))
  } else {
    d_idx <- spec$d_idx
  }
  
  if (is.null(spec$fs)){
    fs <- 0.8
  } else {
    fs <- spec$fs
  }
  
  if (is.null(spec$prior_idx)){
    prior_idx <- c(1:dim(res[[1]][[1]])[2])
  } else {
    prior_idx <- spec$prior_idx
  }
  
  if (is.null(spec$p_idx)){
    p_idx <- c(1:length(res[[1]]))
  } else {
    p_idx <- spec$p_idx
  }
  
  if (is.null(spec$prior_label)){
    prior_label <- colnames(res[[1]][[1]])
  } else {
    prior_label <- spec$prior_label
  }
  
  if (is.null(spec$p_label)){
    p_label <- names(res[[1]])
  } else {
    p_label <- spec$p_label
  }
  
  if (is.null(spec$d_label)){
    d_label <- names(res)
  } else {
    d_label <- spec$d_label
  }
  
  if (is.null(spec$op_density)) {
    op_density <- 2
  } else {
    op_density <- spec$op_density
  }
  
  na.rm <- ifelse(is.null(spec$na.rm), FALSE, spec$na.rm)
  
  nprior <- length(prior_idx); np <- length(p_idx); nd<-length(d_idx)
  
  mypalette <- wes_palette("Zissou1", 6, type = "continuous")
  
  p_label_idx <- which(p_labels_all %in% p_label)
  if (length(p_label_idx)==length(p_label)){
    # map colors -- this is for plotting
    mycols_idx <- data.frame(p_label=p_label) %>%
      left_join(color_mapping, join_by(p_label)) %>%
      select(color) %>% pull()
    mycols <- mypalette[mycols_idx]
    
    # map latex legend labels
    p_labels_latex <- p_labels_all_latex[p_label_idx]
    
    p_label_ordered <- p_labels_all[p_labels_all %in% p_label]
    mycols_idx_ordered <- data.frame(p_label=p_label_ordered) %>%
      left_join(color_mapping, join_by(p_label)) %>%
      select(color) %>% pull()
    mycols_ordered <- mypalette[mycols_idx_ordered]
  } else {
    p_labels_latex <- p_label
    mycols_ordered <- mypalette[default_col_order]
    mycols <- mypalette[default_col_order]
  }
  
  # plotting
  
  # get ymax based on all results so that the scale is suitable for all plots
  res_ymax <- numeric(nd)
  ymax_tmp <- numeric(nprior)
  if (is.null(spec$ymax)){
    for (i in 1:nd){
      dm <- d_idx[i]
      for (j in 1:nprior){
        prior <- prior_idx[j]
        tmp <- sapply(res[[dm]], function(x) x[,prior])
        ymax_tmp[j] <- round(getrange(tmp,op_density,na.rm),2)
      }
      res_ymax[i] <- max(ymax_tmp)
    }
  }

  par(mfrow=c(nd,nprior),mar=c(0.5,0,1.5,1), mgp=c(1.75,0.75,0), oma=c(4,4,1,1))
  for (i in 1:nd){
    dm <- d_idx[i]
    if (is.null(spec$ymax)){
      ymax <- res_ymax[i]
    } else {
      ymax <- spec$ymax
    }
    for (j in 1:nprior){
      prior <- prior_idx[j]

      for (k in seq_along(p_idx)){
        res_i <- res[[i]][[p_idx[k]]][,prior]
        if (na.rm) res_i <- res_i[!is.na(res_i)]
        
        if (op_density==1){
          est_density <- density(res_i)
        } else if (op_density==2) {
          # tested different options of mu under the boundary kernel
          # the default (mu=1) seems to be the most smooth
          est_density <- bde(res_i, estimator="boundarykernel", options = list(mu=1)) 
        } else if (op_density==3){
          est_density <- bd_density(res_i)
        }
        
        if (k==1){
          plot(est_density, xlim=c(0,1), type="l",xlab="",lwd=2, main="",
               col=mycols[1],ylim=c(0,ymax),ylab="",
               xaxt="n",yaxt="n",cex.axis=fs, cex.lab=fs)
        } else {
          lines(est_density, col=mycols[k],lwd=2)
        }
      }
      
      #if (i==1 & j==nprior)
      #  legend(legend_loc, legend = p_labels_latex, y.intersp=0.75,
      #         col = mycols_ordered, lwd=rep(2,np), bty="n", cex = fs)
      if (i==1) title(main=list(prior_label[prior],cex=fs*1.1),line=0.25,font.main=1)
      if (j==1) mtext(side=2, text=d_label[i], line=1.75, cex=fs*0.8)
      if (i==nd) axis(side = 1, at = seq(0,1,0.5), labels = T, cex.axis=fs)
      if (j==1) axis(side=2, labels=T, cex.axis=fs)
    }
  }
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = 'l', bty = 'n', xaxt = 'n', yaxt = 'n')
  legend('bottom', legend=p_labels_latex[p_idx], col=mycols_ordered[c(1:length(p_idx))], lwd=rep(2,np), 
         xpd = TRUE, horiz = TRUE, cex = fs*1.1, seg.len=2, bty = 'n')
  
}

