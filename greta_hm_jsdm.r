
# Negative binomial in greta
# With a hierarchical structure on the rows (alpha)
# Can I put true variances on the hierarchy and estimate them correctly? 

library(greta)
library(geiger)
#library(extraDist)
library(mvtnorm)
library(doMC)
library(foreach)

set.seed(666)

N_cores <- N_chains <- 5
registerDoMC(cores=N_cores)
N_dataSet <- 20

n_samples = 20000
warmup = 15000
thin = 40

## Data & deterministic parameters - these are fixed across simulations!
n.species <- 5 # number of OTUs/ASVs
n.hosts <- 20 # number of host species
n.latent <- 2 # number of latent factors
n.samples <- n.hosts*5 # number of samples
hostID <- rep(1:n.hosts,5)

### Simulate the host phylogeny
birth.rate = 1
death.rate = 0
phy <- geiger::drop.extinct(geiger::sim.bdtree(b = birth.rate, d = death.rate, stop = "taxa", n = n.hosts, extinct = FALSE))
phy$tip.label = 1:n.hosts
vphy <- ape::vcv(phy) # var-cov matrix
vphy <- vphy/(det(vphy)^(1/n.hosts)) # standardize the matrix
C <- cov2cor(vphy) # corr matrix

### Simulate factor loadings

## Sample lvl
Loading_true <- matrix(0,nrow=n.species,ncol=n.latent)
for(k in 1:n.species) {Loading_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))} 

## Host-species lvl
Loading_H_true <- matrix(0,nrow=n.species,ncol=n.latent)
for(k in 1:n.species) {Loading_H_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))} 

### Fix variance parameters for the different row effects

sigma.sample_true <- 6 # quantifies variance that can be attributed to the sample level
sigma.host_true <- 2 # quantifies variance that can be attributed to the host-species level
scale.phylo_true <- rexp(1, rate = 0.1) # quantifies variance that can be attributed to the phylogenetic effect

## Total variance of row effects
tot.var.alpha_true <- sigma.sample_true^2 + sigma.host_true^2 + scale.phylo_true^2

## Function to check for overdispersion (i.e., var > mean)
nb.overdispersion <- function(eta, phi, returnVals=FALSE){ 
  nb.mean <- apply(eta, 2, mean)
  nb.var <- nb.mean + phi*nb.mean^2 
  bool <- nb.var > nb.mean
  if(returnVals==TRUE){
    return(rbind(mean=nb.mean, var=nb.var))
  }else{
    return(bool)
  }
}

dosim <- function(seed_num) {
  ############################
  ###### Generate data ######
  ############################
  set.seed(seed_num)
  
  eta_true <- z_true <- y <- matrix(0, nrow = n.samples, ncol = n.species)
  
  gamma_true <- runif(n.species, -1, 1) # smaller range of species-specific intercepts; otherwise you can too rich and too poor prevalent species
  
  ## Simulate the hierarchical structure on the row effects alpha

  host_phylo.eff_true <- mvtnorm::rmvnorm(1, mean = rep(0,n.hosts), sigma = C) # the phylogenetic host effect
  
  mu.host_true <- rnorm(n.hosts, 0, 1)
  host.eff_true <-  rnorm(n.hosts, mean = mu.host_true, sd = sigma.host_true)
  host.mean_true <- host.eff_true + host_phylo.eff_true*scale.phylo_true
  alpha_true <-  rnorm(n.samples, mean = host.mean_true[hostID], sd = sigma.sample_true)
  
  ### Simulate latent factors
  
  ## Factor loadings are defined above and fixed across simulations
  
  ## Sample lvl
  LV_true <- matrix(0,nrow=n.samples,ncol=n.latent)
  for(k in 1:n.samples) {LV_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))}
  
  ## Host-species lvl
  LV_H_true <- matrix(0,nrow=n.hosts,ncol=n.latent)
  for(k in 1:n.hosts) {LV_H_true[k,] <- mvtnorm::rmvnorm(1, mean = rep(0,n.latent))}
  
  epsilon_true <- as.matrix(LV_true)%*%t(as.matrix(Loading_true))
  epsilon_H_true <- as.matrix(LV_H_true)%*%t(as.matrix(Loading_H_true))

  ## Simulate compositional count data under the NB distribution 
  
  phi_true <- rep(0.1, n.species) # overdispersion parameter
  
  for(i in 1:n.samples) {
    for(j in 1:n.species) {
      eta_true[i,j] <- alpha_true[i] + gamma_true[j] + epsilon_true[i,j] + epsilon_H_true[hostID[i],j]
      #y[i,j] <- rnbinom(n = 1, prob = 1/phi_true[j]/(exp(eta_true[i,j])+1/phi_true[j]), size = 1/phi_true[j])
      y[i,j] <- rnbinom(n = 1, mu = exp(eta_true[i,j]), size = 1/phi_true[j]) # alternative parametrization via mean
    } 
  }
  
  y[is.na(y)] <- 0 # if NA creeps in; sometimes happens when the mean gets too large
  
  nb_overdisp1 <- nb.overdispersion(eta_true,phi_true,returnVals=T) # check overdispersion

  ## Put true params in list to save with the mcmc samples 
  true_params <- list( nb_overdisp1=nb_overdisp1 
                      ,eta_true=eta_true
                      ,alpha_true=alpha_true
                      ,gamma_true=gamma_true
                      ,host.mean_true=host.mean_true
                      ,host.eff_true=host.eff_true
                      ,host_phylo.eff_true=host_phylo.eff_true
                      ,LV_true=LV_true
                      ,LV_H_true=LV_H_true
                      ,Loading_true=Loading_true
                      ,Loading_H_true=Loading_H_true
                      ,tot.var.alpha_true=tot.var.alpha_true
                      ,sigma.sample_true=sigma.sample_true
                      ,sigma.host_true=sigma.host_true
                      ,scale.phylo_true=scale.phylo_true )
  
  ############################
  #### Model & estimation ####
  ############################

  z_host_phylo.eff = greta::normal(0, 1, dim = n.hosts)
  host_phylo.eff <- chol(C) %*% z_host_phylo.eff
  scale.phylo = greta::exponential(0.1)
  
  # non-centered parameteriztion of mu.host
  mu.host = greta::normal(0, 5, dim = n.hosts)
  z_host.eff = greta::normal(0, 1)
  sigma.host = greta::cauchy(0, 1, truncation = c(0,Inf))
  host.eff = mu.host + sigma.host*z_host.eff
  
  # non-centered parameteriztion of alpha0
  host.mean <- host.eff + host_phylo.eff*scale.phylo # the row effect's linear predictor
  z_alpha = greta::normal(0, 1, dim = n.samples)
  sigma.sample = greta::cauchy(0, 1, truncation = c(0,Inf))
  alpha0 = host.mean[hostID] + sigma.sample*z_alpha
  
  ## Total variance of row effects
  tot.var.alpha <- sigma.sample^2 + sigma.host^2 + scale.phylo^2
  
  var.host_phylo <- scale.phylo^2
  var.host <- sigma.host^2
  var.sample <- sigma.sample^2
  
  # latent factors
  LV = greta::normal(0, 1, dim = c(n.samples, n.latent))
  LV_H = greta::normal(0, 1, dim = c(n.hosts, n.latent))
  
  # corner constraints on the loading matrices
  Loading <- Loading_H <- greta::zeros(n.species, n.latent)
  diag(Loading) = greta::normal(0, 1,dim = c(n.latent), truncation=c(0,Inf))
  Loading[lower.tri(Loading)] = greta::normal(0, 1,dim = length(Loading[lower.tri(Loading, diag = F)]))
  diag(Loading_H) = greta::normal(0, 1,dim = c(n.latent),truncation = c(0, Inf))
  Loading_H[lower.tri(Loading_H)] = greta::normal(0, 1, dim=length(Loading_H[lower.tri(Loading_H, diag = F)]))
  
  gamma0 = greta::normal(0, 1, dim = c(n.species))
  
  eta = LV%*%t(Loading) + LV_H[hostID,] %*% t(Loading_H)
  eta = greta::sweep(eta, 1, alpha0, "+")
  eta = greta::sweep(eta, 2, gamma0, "+")
  
  # overdispersion parameter distributed as a half-cauchy (following Polson & Scott 2012)
  phi = greta::cauchy(0, 2.5, truncation=c(0,Inf))
  
  #direct negbin parameterization
  expeta <- exp(eta)
  p <- 1/phi/(expeta+1/phi)
  greta::distribution(y) = greta::negative_binomial(size = 1/phi, prob = p)
  
  ## Build model
  m_fit <- model(  eta
                  ,expeta
                  ,gamma0
                  ,alpha0
                  ,host.mean
                  ,host.eff
                  ,host_phylo.eff
                  ,sigma.sample
                  ,sigma.host
                  ,scale.phylo
                  ,tot.var.alpha
                  ,var.sample
                  ,var.host
                  ,var.host_phylo
                  ,LV, LV_H
                  ,Loading, Loading_H
                  ,precision = "double" ) #, n_cores=1 ) #this n_cores neq N_cores above!
  
  ## Sampling
  draws <- mcmc(m_fit, n_samples=n_samples, warmup=warmup, thin=thin, verbose=T, chains=N_chains)

  out <- list( draws=draws, true_params=true_params, C=C )
  save(out,file=paste0("host_microbiota_run1.",seed_num,".RData"))
  return(out)      
}

start_time <- Sys.time()
out <- foreach(i=1:N_dataSet) %dopar% dosim(seed_num = i)
end_time <- Sys.time()
time_taken <- end_time - start_time

q('no')

###############################################################################################
###################################### ANALYZE OUTPUT #########################################
###############################################################################################

eta_draws <- draws[[1]][,grep("^eta\\d+",colnames(draws[[1]]))]
alpha0_draws <- draws[[1]][,grep("alpha0\\d+",colnames(draws[[1]]))]
host.mean_draws <- draws[[1]][,grep("host.mean\\d+",colnames(draws[[1]]))]
host.eff_draws <- draws[[1]][,grep("host.eff\\d+",colnames(draws[[1]]))]
host_phylo.eff_draws <- draws[[1]][,grep("host_phylo.eff\\d+",colnames(draws[[1]]))]

LV_draws <- draws[[1]][,grep("LV\\d+",colnames(draws[[1]]))]
LV_H_draws <- draws[[1]][,grep("LV_H\\d+",colnames(draws[[1]]))]
Loading_draws <- draws[[1]][,grep("Loading\\d+",colnames(draws[[1]]))]
Loading_H_draws <- draws[[1]][,grep("Loading_H\\d+",colnames(draws[[1]]))]

var.host_draws <- draws[[1]][,grep("var.host$",colnames(draws[[1]]))]
var.sample_draws <- draws[[1]][,grep("var.sample",colnames(draws[[1]]))]
tot.var.alpha_draws <- draws[[1]][,grep("tot.var.alpha",colnames(draws[[1]]))]
var.host_phylo_draws <- draws[[1]][,grep("var.host_phylo",colnames(draws[[1]]))]

mcmcSetting <- "(mcmc:20000/15000/40)"

plot(coda::effectiveSize(eta_draws));title(paste("eta effectSize",mcmcSetting))
plot(coda::effectiveSize(alpha0_draws));title(paste("alpha0 effectSize",mcmcSetting))
plot(coda::effectiveSize(host.mean_draws));title(paste("host.mean effectSize",mcmcSetting))
plot(coda::effectiveSize(host_phylo.eff_draws));title(paste("host_phylo.eff effectSize",mcmcSetting))
plot(coda::effectiveSize(LV_draws));title(paste("LV effectSize",mcmcSetting))
plot(coda::effectiveSize(LV_H_draws));title(paste("LV_H effectSize",mcmcSetting))
plot(coda::effectiveSize(Loading_draws));title(paste("Loading effectSize",mcmcSetting))
plot(coda::effectiveSize(Loading_H_draws));title(paste("Loading_H effectSize",mcmcSetting))

post.eta <- vector("list",500)
post.alpha0 <- vector("list",500)
post.host.eff <- vector("list",500)
post.host.mean <- vector("list",500)
post.host_phylo.eff <- vector("list",500)

post.LV <- vector("list",500)
post.LV_H <- vector("list",500)
post.Loading <- vector("list",500)
post.Loading_H <- vector("list",500)

post.var.host_draws <- vector("list",500)
post.var.sample_draws <- vector("list",500)
post.tot.var.alpha_draws <- vector("list",500)
post.var.host_phylo_draws <- vector("list",500)

for(i in 1:500){
  post.eta[[i]] <- matrix(unlist(split(eta_draws[i,], ceiling(seq_along(eta_draws[i,])/n.samples))),ncol=n.samples)
  post.alpha0[[i]] <- matrix(unlist(split(alpha0_draws[i,], ceiling(seq_along(alpha0_draws[i,])/n.samples))),ncol=n.samples)
  post.host.mean[[i]] <- matrix(unlist(split(host.mean_draws[i,], ceiling(seq_along(host.mean_draws[i,])/n.hosts))),ncol=n.hosts)
  post.host.eff[[i]] <- matrix(unlist(split(host.eff_draws[i,], ceiling(seq_along(host.eff_draws[i,])/n.hosts))),ncol=n.hosts)
  
  post.host_phylo.eff[[i]] <- matrix(unlist(split(host_phylo.eff_draws[i,], ceiling(seq_along(host_phylo.eff_draws[i,])/n.hosts))),ncol=n.hosts)
  
  post.LV[[i]] <- matrix(unlist(split(LV_draws[i,], ceiling(seq_along(LV_draws[i,])/n.samples))),ncol=n.samples)
  post.LV_H[[i]] <- matrix(unlist(split(LV_H_draws[i,], ceiling(seq_along(LV_H_draws[i,])/n.hosts))),ncol=n.hosts)
  post.Loading[[i]] <- matrix(unlist(split(Loading_draws[i,], ceiling(seq_along(Loading_draws[i,])/n.species))),ncol=n.species)
  post.Loading_H[[i]] <- matrix(unlist(split(Loading_H_draws[i,], ceiling(seq_along(Loading_H_draws[i,])/n.species))),ncol=n.species)
}

post.eta.mean <- Reduce("+", post.eta) / length(post.eta)
post.alpha.mean <- Reduce("+", post.alpha0) / length(post.alpha0)
post.host.mean <- Reduce("+", post.host.mean) / length(post.host.mean)
post.host_phylo.eff.mean <- Reduce("+", post.host_phylo.eff) / length(post.host_phylo.eff)
post.host.eff.mean <- Reduce("+", post.host.eff) / length(post.host.eff)

post.LV.mean <- Reduce("+", post.LV) / length(post.LV)
post.LV_H.mean <- Reduce("+", post.LV_H) / length(post.LV_H)
post.Loading.mean <- Reduce("+", post.Loading) / length(post.Loading)
post.Loading_H.mean <- Reduce("+", post.Loading_H) / length(post.Loading_H)

plot(c(post.eta.mean),c(eta_true));abline(0,1) 
plot(c(post.alpha.mean),c(alpha_true));abline(0,1) 
plot(c(post.host.mean),c(host.mean_true));abline(0,1)
plot(c(post.host.eff.mean),c(host.eff_true));abline(0,1) 
plot(c(post.host_phylo.eff.mean),c(host_phylo.eff_true));abline(0,1) 

plot(c(post.LV.mean),c(LV_true));abline(0,1) 
plot(c(post.LV_H.mean),c(LV_H_true));abline(0,1) 
plot(c(post.Loading.mean),c(Loading_true));abline(0,1) 
plot(c(post.Loading_H.mean),c(Loading_H_true));abline(0,1) 

tot.var.alpha_draws <- draws[[1]][,grep("tot.var.alpha",colnames(draws[[1]]))]
var.sample_draws <- draws[[1]][,grep("var.sample",colnames(draws[[1]]))]
var.host_draws <- draws[[1]][,grep("var.host$",colnames(draws[[1]]))]
var.host_phylo_draws <- draws[[1]][,grep("var.host_phylo",colnames(draws[[1]]))]

boxplot(data.frame(x=tot.var.alpha_draws));abline(h=tot.var.alpha_true,col="red")
boxplot(data.frame(x=var.sample_draws));abline(h=sigma.sample_true^2,col="red")
boxplot(data.frame(x=var.host_draws),ylim=c(0,5));abline(h=sigma.host_true^2,col="red")
boxplot(data.frame(x=var.host_phylo_draws));abline(h=scale.phylo_true^2,col="red")

  
  