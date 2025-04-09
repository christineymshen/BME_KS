slurm_id <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
n_cpus <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))

outfolder <- ""

simnum <- 2
folder <- "gamma" # gamma model

# p-value to be calculated
p_label <- "ptm1"

if (!file.exists(paste0(outfolder, folder))) dir.create(paste0(outfolder, folder))

output_name <- paste0(outfolder,folder,"/sim", simnum, "_", p_label,"_",slurm_id, ".rds")

if (!file.exists(output_name)){
  
  library(rstan)
  library(tidyverse)
  library(doMC)
  
  registerDoMC(cores = n_cpus)
  source("functions.R")
  
  # fixed parameters
  burnin <- 1000
  thin <- 5; nchains <- 1
  nsim <- 1000 # number of MCMC posterior samples/ posterior predictive rep datasets
  nsim1 <- rep(1,nsim)
  
  # number of test statistics/ discrepancy measures
  d_labels <- c("KS")
  nd <- length(d_labels)
  
  # general parameters
  n <- 20
  
  # true parameters
  a <- 2
  b <- 0.5
  
  # data.frame(x = 0:10) %>%
  #   ggplot(aes(x = x)) +
  #   stat_function(fun = dgamma, args=list(shape=4,rate=1)) +
  #   theme_bw()
  
  # simulate observed data
  set.seed(slurm_id)
  y <- rgamma(n,a,b)
  
  # prior parameter
  prior_labels <- c("good","bad")
  nprior <- length(prior_labels)
  
  # prior for a
  mus <- c(2.5,1)
  sigmas <- c(4,0.5)
  a0s <- c(1,3)
  b0s <- c(1,1.25)

  mod <- readRDS("gamma_mod.rds")
  output <- matrix(0,nprior,nd,dimnames = list(prior_labels,d_labels))
  
  runtime <- system.time({
    for (q in 1:nprior){
    
      data <- list(n=n,y=y,mu=mus[q],sigma=sigmas[q],a0=a0s[q],b0=b0s[q])
      stan_mod <- sampling(mod, data = data,
                           chains = nchains, warmup = burnin,iter = (burnin+thin*nsim/nchains),
                           cores = nchains, thin = thin, seed = slurm_id, refresh=0)
      p_a <- as.matrix(stan_mod, pars="a")
      p_b <- as.matrix(stan_mod, pars="b")
      
      ptm_obs <- ks.test(y,pgamma,mean(p_a),mean(p_b))$statistic
      
      y_rep <- matrix(rgamma(nsim*n,p_a,p_b),nrow=nsim)
      
      ptm_rep <- foreach(i=1:nsim, .combine="rbind") %dopar% {
        data$y <- y_rep[i,]
        stan_mod_i <- sampling(mod, data = data,
                               chains = nchains, warmup = burnin,iter = (burnin+thin*nsim/nchains),
                               cores = nchains, thin = thin, seed = slurm_id, refresh=0)
        p_a_i_m <- mean(as.matrix(stan_mod_i, pars="a"))
        p_b_i_m <- mean(as.matrix(stan_mod_i, pars="b"))
        return(ks.test(y_rep[i,],pgamma,p_a_i_m,p_b_i_m)$statistic)
      }
      output[q,] <- pr_ge(ptm_rep,outer(nsim1, ptm_obs))
    }
  })
  
  res <- list(runtime=runtime,output=output)
  saveRDS(res, output_name)
}

