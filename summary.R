source("functions.R")

# to get summary statistics from run results
for (i in 1:8){
  agg(i,"gamma","ptm1")
  agg(i,"gamma","ptm2")
}

# to get plots
runname <- "gamma"
labels_out <- c("n=10","n=20","n=100","n=500")

currfolder <- ""
res_ptm1_level <- read_res(currfolder,runname,"ptm1",c(1:4),labels_out)
res_ptm2_level <- read_res(currfolder,runname,"ptm2",c(1:4),labels_out)

res_ptm_level <- combine(labels=c("good prior w posterior mean", "bad prior w posterior mean",
                                  "good prior w MLE", "bad prior w MLE"),
                         res_ptm1_level$output,res_ptm2_level$output)
spec <- vector("list")
spec$fs <- 1.3

pdf(file=paste0(currfolder,"/fig/gm_ptm_level.pdf"), width=10, height=2.5)
plotp_density4(res_ptm_level,spec)
dev.off()

res_ptm1_5 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",5,"_","ptm1",".rds"))$output)
res_ptm2_5 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",5,"_","ptm2",".rds"))$output)
res_ptm1_6 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",6,"_","ptm1",".rds"))$output)
res_ptm2_6 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",6,"_","ptm2",".rds"))$output)
res_ptm1_7 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",7,"_","ptm1",".rds"))$output)
res_ptm2_7 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",7,"_","ptm2",".rds"))$output)
res_ptm1_8 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",8,"_","ptm1",".rds"))$output)
res_ptm2_8 <- list(readRDS(paste0(currfolder,"/res/", runname,"_",8,"_","ptm2",".rds"))$output)

res_ptm1_power <- combine(labels=c("null", "gamma GLM", "weibull", "lognormal"),
                         res_ptm1_5, res_ptm1_7, res_ptm1_6, res_ptm1_8)
res_ptm2_power <- combine(labels=c("null", "gamma GLM", "weibull", "lognormal"),
                          res_ptm2_5, res_ptm2_7, res_ptm2_6, res_ptm2_8)
spec <- vector("list")
spec$fs <- 1.3
spec$ymax <- 5

pdf(file=paste0(currfolder,"/fig/gm_ptm1_power.pdf"), width=10, height=2.5)
plotp_density4(res_ptm1_power,spec)
dev.off()

pdf(file=paste0(currfolder,"/fig/gm_ptm2_power.pdf"), width=10, height=2.5)
plotp_density4(res_ptm2_power,spec)
dev.off()
