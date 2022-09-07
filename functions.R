# helper functions for combined evidence imbalance experiment #
require(phytools)
require(RPANDA)
require(CollessLike)
require(parallel)
require(TreeSimGM)
require(ggplot2)
require(tidyverse)
require(reshape2)
require(TreeSim)
require(phangorn)
require(knitr)
require(stringr)
require(tidyverse)
require(BAMMtools)
require(R.devices)
require(caper)
require(DDD)

setwd("~/Dropbox/imbalance_combined_evidence")

# calculate colless index for multiphylo object and return list of scores
ci.norm <- function(x){ colless.like.index(tree = x, norm = TRUE)}
#calculate Sakin index for multiphylo object and return list of scores
si.norm <- function(x){ sackin.index(x, norm = T)}

#calculate the distance of a tip to the current root

tree.length <- function(tree) {
  max(diag(vcv(tree)))
}

# count the number of extinct lineages on the tree
count_dead <- function(phy){
  count <- length(getExtinct(phy, tol = 0.001))
  return(count)
}

# count the number of living lineages on the tree
count_living <- function(phy){
  count <- length(getExtant(phy, tol = 0.001))
  return(count)
}

# prouces a list of a trees tips with their corresponding node labels
tips_w_nodes <- function(phy){
  t <- ladderize(phy, right = F)
  tips <- t$edge[,2] <= length(t$tip.label)
  ordered_tips <- t$edge[tips, 2]
  names(ordered_tips) <- t$tip.label[ordered_tips]
  return(ordered_tips)
}

# give it a list of tip names and get the heights back
tip_heights <- function(list, phy) {
  tips <- tips_w_nodes(phy = phy)
  tips <- tips[list]
  heights <- vector()
  for(i in 1:length(tips)){
    heights <- c(heights, nodeheight(phy, tips[[i]]))
  }
  names(heights) <- names(tips)
  return(heights)
}

#and we need something to fit a bd model to out tree and estimate lambda and mu


#just a general function for building trees so I can test that these functions work
tree.build <- function(n, nsim, par) {sim.bd.taxa(n = n, numbsim = nsim, lambda = par[1], mu = par[2])}
par <- c(0.8, 0.3)
tree <- tree.build(n = 100, nsim = 1, par = par)[[1]]

# simulates set of trees where the number of extinct species are always equal to or greater than the empirical tree
# calls simulation function
make.trees <- function(FUN, phy, nsim, par, e = count_dead(phy), n = count_living(phy))
  {
    #t <- FUN(n = count_living(phy), nsim = nsim, par = par)
    t <- FUN(n = n, nsim = nsim, par = par)
    targets <- which(unlist(mclapply(t, count_dead, mc.cores = 3)) < e)
    while(length(targets) > 0){
        new.t <- FUN(n = count_living(phy), nsim = length(targets), par = par)
        t[targets] <- new.t
        targets <- which(unlist(mclapply(t, count_dead, mc.cores = 3)) < e)
    }
    return(t)
  }

#trees <- make.trees(tree.build, tree = tree, nsim = 1000, par = par)
#which(unlist(mclapply(trees, count_dead, mc.cores = 3)) < count_dead(tree))

#ci. <- unlist(mclapply(trees, ci.norm, mc.cores = 3))
#hist(ci)
#cool... so now we can make this a function

#just here to help check the rest of the function
#ci <- unlist(mclapply(trees, ci.norm, mc.cores = 3))



#sample fossils randomly
fossil.random.sampling <- function(phy, e){ #t is the list tree and e is the number of extinct lineages to be sampled
  #randomly sample fossil tips
  s <- sample(getExtinct(phy, tol = 0.001), e)
  s <- c(s, getExtant(phy, tol = 0.001))
  tr <- keep.tip(phy, s)
  return(tr)
}



#sample fossil according to an evolved (BM) trait  
fossil.bm.sampling <- function(phy, e) {
  s <- fastBM(phy)
  s <- s[getExtinct(phy, tol = 0.001)]
  x <- runif(1, range(s)[1], range(s)[2])
  s <- sort(abs(s - x), decreasing = F) #sample a random point in range of BM values
  s <- names(s)[1:e]
  s <- c(s, getExtant(phy, tol = 0.001))
  tb <- keep.tip(phy, s)
  tb$trait_value <- x
  return(tb)
}


# sample fossils to maximise phylogenetic distance
fossil.pd.sampling <- function(phy, e, method = "bootstrap", reps = 1000) {
  extinct.tree <- keep.tip(phy, getExtinct(phy, tol = 0.001))
  if(method == "bootstrap"){
  s <- pd_bootstrap_sampler(tree = extinct.tree, ntips = round(e), reps = reps, method = "TBL", tip.weights = NULL)
  s <- s[sample(which(s[,1] == max(s[,1])), 1),]
  #s <- as.character(s[-1])
  tp <- keep.tip(phy, tip = c(getExtant(phy, tol = 0.001),as.character(s[-1])))
  tp$pd_scores <- s[[1]]
  return(tp)
  } else if (method == "pda"){
  file <- tempfile(pattern = "tree", tmpdir = "pd")
  print(file)
  write.tree(extinct.tree, file = paste0(file, ".tree"))
  execute <- combine_words(c("./pd/pda", paste0(file, ".tree"), "-k", e, "-g", paste0(file, ".txt")), sep = " ")
  system(execute)
  #read in the output of pda
  s <- try(read.delim(file = paste0(file, ".txt")))
 # if(class(s) == 'try-error'){
#  unlink(paste0(file, ".tree"))
#    unlink(paste0(file, ".tree.log"))
#    unlink(paste0(file, ".txt"))
#    } else {
  s <- s[[1]] #and some processing
  x <- as.character(s[[11]])
  x <- as.numeric(str_split(x, " ")[[1]])[10]
  x <- x/sum(extinct.tree$edge.length)
  s <- as.character(s[13:(e+12)]) #to extract just the names of our sample tips
  s <- c(s, getExtant(phy, tol = 0.001))
  tp <- keep.tip(phy, s)
  tp$pd_scores <- x
  unlink(paste0(file, ".tree"))
  unlink(paste0(file, ".tree.log"))
  unlink(paste0(file, ".txt"))
  return(tp)
  }
    }
  #}



# sample fossils according to their proximity in a random point in time
fossil.time.sampling <- function(phy, e) {
  s <- tip_heights(list = getExtinct(phy, tol = 0.001), phy = phy)
  x <- runif(1, range(nodeHeights(phy)[,1])[1], range(nodeHeights(phy)[,1])[2]) #sample a random point in time across our tree
  s <- sort(abs(s - x)) #order our fossil tips based on their distance in time from the sample 
  s <- names(s)[1:e]
  s <- c(s, getExtant(phy, tol = 0.001))
  tt <- keep.tip(phy, s)
  tt$time_point <- x
  #print(x)
  return(tt)
}

#some functions for extracting variables of interest
get_pd <- function(phy) {feat <- phy$pd_scores; return(feat)}
get_tp <- function(phy) {feat <- phy$time_point; return(feat)}
get_tv <- function(phy) {feat <- phy$trait_value; return(feat)}

# function that takes multiphylo list, applies sampling function
multiphylo_sampler <- function(mphylo, sampling_function, e)
{
  #we need to create a little worker function so we can run this in parallel
  sampler <- function(phy){sampling_function(phy, e = e)}
  trees <- (mclapply(mphylo, FUN = sampler, mc.cores = 3))
  return(trees)
}

#a function for bootstrapping pd values from a tree
pd_bootstrap_sampler <- function (tree, ntips, reps = 1000, method = "TBL", tip.weights = NULL) 
{
  method <- match.arg(method, c("TBL", "MST", "UEH", "SBL", 
                                "TIP"))
  if (class(tree) != "clade.matrix") {
    if (class(tree) == "phylo") {
      warning("Converting phylo object to clade.matrix object")
      cm <- clade.matrix(tree)
    }
    else {
      stop("pd.calc requires a phylogeny")
    }
  }
  total.nb.tips <- dim(cm$clade.matrix)[2]
  if (!(ntips %in% 1:(total.nb.tips - 1))) {
    stop("'sample' must be a positive integer lower than the number of tips")
  }
  pd.store <- numeric(reps)
  tips <- 1:total.nb.tips
  if (!is.null(tip.weights)) {
    if (!is.null(names(tip.weights))) {
      wght.match <- match(cm$tip.label, names(tip.weights))
      if (any(is.na(wght.match))) {
        warning("The returned tip labels have no matching named element in tip.weights")
        return(cm$tip.label[is.na(wght.match)])
      }
      tip.weights <- tip.weights[wght.match]
    }
    else {
      stop("'weights' must be a vector of weights, named to match the tip labels")
    }
  }
  pd_samples <- matrix(nrow = 0, ncol = ntips+1)
  for (rep in seq(along = pd.store)) {
    which.tips <- sample(tips, ntips, prob = tip.weights)
    pd <- pd.calc(cm, tip.subset = which.tips, 
                  method = method)
    pd_samples <- rbind(pd_samples, c(pd[1], tree$tip.label[which.tips]))
  }
  
  pd_samples <- (as.tibble(pd_samples))
  colnames(pd_samples) <- c("PD", paste0("tip_", seq(1:ntips)))
  pd_samples$PD <- as.numeric(pd_samples$PD)
  return(pd_samples)
}



#now a function that applies all four version of the above, calculates CIs, makes a tibble

get_data <- function(mphylo, e){
  #apply sampling procedures to multiphylo list and calculate CI for each
  random <- multiphylo_sampler(mphylo = mphylo, sampling_function = fossil.random.sampling, e = e)
    random <- unlist(mclapply(random, ci.norm, mc.cores = 3))
    print("I've finsihed the random sample")
  brownian <- multiphylo_sampler(mphylo = mphylo, sampling_function = fossil.bm.sampling, e = e)
    trait_value <- unlist(mclapply(brownian, get_tv, mc.cores = 3))
    brownian <- unlist(mclapply(brownian, ci.norm, mc.cores = 3))
    print("and now the brownian")
  phylogenetic_diversity <- multiphylo_sampler(mphylo = mphylo, sampling_function = fossil.pd.sampling, e = e)
    pd_scores <- unlist(mclapply(phylogenetic_diversity, get_pd, mc.cores = 3)) 
    phylogenetic_diversity <- unlist(mclapply(phylogenetic_diversity, ci.norm, mc.cores = 3))
    print("pd is done too")
  time_dependent <- multiphylo_sampler(mphylo = mphylo, sampling_function = fossil.time.sampling, e = e)
    time_point <- unlist(mclapply(time_dependent, get_tp, mc.cores = 3))
    time_dependent <- unlist(mclapply(time_dependent, ci.norm, mc.cores = 3))
    print("and finally time dependent")
  #get some other variables we may be interested in
  n_fossils <- unlist(mclapply(mphylo, count_dead, mc.cores = 3))
  tree_lengths <- unlist(mclapply(mphylo, tree.length, mc.cores = 3))
  
  shift_no <- vector()
    for(i in 1:length(mphylo)){
      shift_no <- c(shift_no, mphylo[[i]]$shifts)
    }
  
  df <- tibble(tree = seq(1, length(mphylo), 1), random, brownian, phylogenetic_diversity, time_dependent, n_fossils, 
               tree_lengths, shift_no, trait_value, time_point, pd_scores)
#you could add PD, BM and time point scores to this too and maybe length of the tree
  df <- df %>% pivot_longer(-c(tree, n_fossils, tree_lengths, shift_no, trait_value, time_point, pd_scores), names_to = "sampling", values_to = "CI")
  #df$sampling <- as.factor(df$sampling)
  }

#takes the output of get_data and calculates the number of runs where CI exceeds a given threshold by sampling regime
per_past_threshold <- function(data, threshold, stat, models){
#  models <- as.factor(data[,models])
  levels <- levels(data$models)
  excess <- vector()
  for(i in 1:length(levels)){
  s <- which(data[,models] == levels[i])
  s <- length(which(data[s, stat] >= threshold))
  excess <- c(excess, s)
  }
  names(excess) <- levels
  return(excess)
}


get_params <- function(edata, rtt, nsamples = 1, replace = F, names = T){
  #positive_div <- which(rtt$lambda[,1] > rtt$mu[,1]) #we only want samples were net diversification is positive intially, otehrwise TREESIM fails.. IS THIS STILL TRUE
  #positive_div <- seq(1,901,1) 
  positive_div <- which(rtt$lambda[,1] > rtt$mu[,1] & rtt$lambda[,1]/4 < rtt$mu[,1]) #we only want samples were net diversification is positive intially, otehrwise TREESIM fails.. IS THIS STILL TRUE
  samples <- sample(positive_div,nsamples, replace = replace)
  params <- vector()
  for(i in 1:length(samples)){
    lam_init <- rtt$lambda[samples[i],1]
    mu_init <- rtt$mu[samples[i],1]
    shifts <- edata$numberEvents[samples[i]]-1
    if(shifts == 0){
      shift_prob <- 0
      lam_range <- c(0,0)
      mu_range <- c(0,0)
    } else {
      shift_prob <- shifts/(edata$Nnode) #double check that's it
      lam_rates <- unique(edata$tipLambda[[samples[i]]])/rtt$lambda[samples[i],1]
      lam_rates <- lam_rates[-which(lam_rates == 1)]
      if(length(lam_rates) == 1){ 
          lam_range <- c(0, lam_rates)
        }else{
          lam_range <- range(lam_rates)
        }
      mu_rates <- unique(edata$tipMu[[samples[i]]])/rtt$mu[samples[i],1]
      mu_rates <- mu_rates[-which(mu_rates == 1)]
      if(length(mu_rates) == 1){ 
        mu_range <- c(0, mu_rates)
      }else{
        mu_range <- range(mu_rates)
      }   
    }
    params <- rbind(params, c(lam_init, mu_init, shifts, shift_prob, lam_range, mu_range, samples[i]))
  }
  if(names == T){
    colnames(params) <- c("lam_init", "mu_init", "shift_no", "shift_prob", "min_lam_shift", "max_lam_shift", "min_mu_shift", "max_mu_shift", "samples")
  }
  #params <- as_tibble(params)
  return(params)
}



draw_params <- function(nsim, hyper_lambda = runif(nsim, 3, 5), hyper_mu = runif(nsim, .1, .5), mu_as_proportion = F, 
                        shift_no = sample(seq(1,5,1), nsim, replace = T), lambda_scalers = c(0.5, 10), 
                        mu_scalers = c(0.1, 2), N){
  lambda_init <- hyper_lambda #draw the speciation rates from specified distribution
  
  if(mu_as_proportion == T) {
    mu_init <- hyper_lambda*hyper_mu
  } else {
      mu_init <- hyper_mu #draw extinction rates from specified distribution
  }
  shift_prob <- shift_no/(N-1) #the number of desired shifts and the nodes in the desired tree - if trees of different sizes provide list of sizes as argument N
  min_lam_shift <- rep(eval(sort(lambda_scalers)[1]), nsim)
  max_lam_shift <- rep(eval(sort(lambda_scalers)[2]), nsim)
  min_mu_shift <- rep(eval(sort(mu_scalers)[1]), nsim)
  max_mu_shift <- rep(eval(sort(mu_scalers)[2]), nsim)
  df <- data.frame(lambda_init, mu_init, shift_no, shift_prob, min_lam_shift, max_lam_shift,
                min_mu_shift, max_mu_shift)
  return(df)
}


#parameters <- draw_params(nsim = 10, N = 100)

### Simulate trees and measure summary stats for any treebuilder and fossil sampler combination
# params needs to be a dataframe in the format of the output of the get_params() or draw_params() functions



#want to add gsa to this!!!
get_sum_stats <- function(params, N, e, treebuilder, fossilsampler = "none", gsa = FALSE, m, age, extinction_date = NA, survivorship = NA, ddmodel = NA,
                          KM = NA, tinn = NA, shape = NA, e_is_prop = F, sampling_fraction = NULL){
  #t <- tree.build.symmetric(N, 1, par = par)[[1]]
  
  if(is.null(sampling_fraction) == F){
    N <- N/sampling_fraction
  }
  
  t <- pbtree(n = 10)
  if(gsa == FALSE) {
  while(length(getExtinct(t, tol = 0.001)) < e){
    if(treebuilder == "sym_con"){
      t <- tree.build.symmetric(N, 1, par =  params)[[1]] 
    } else if(treebuilder == "asym_con") {
      t <- tree.build.asymmetric(N, 1, par =  params)[[1]] 
    } else if (treebuilder == "sym_shiftsp"){
      t <- tree.build.shiftsp(N, 1, par = params)[[1]] 
    } else if (treebuilder == "asym_shiftsp"){
      t <- tree.build.asymmetric_shiftsp(N, 1, par = params)[[1]] 
    } else if (treebuilder == "sym_shiftex"){
      t <- tree.build.shiftex(N, 1, par = params)[[1]]
    } else if (treebuilder == "asym_shiftex") {
      t <- tree.build.asymmetric_shiftex(N, 1, par = params)[[1]]
    } else if (treebuilder == "sym_shiftsp_shiftex"){
      t <- tree.build.shiftsp_shiftext(N, 1, par = params)[[1]] 
    } else if (treebuilder == "asym_shiftsp_shiftex"){
      t <- tree.build.asymmetric_shiftsp_shiftext(N, 1, par = params)[[1]] 
    } else if (treebuilder == "con_mass_extinction"){
      t <- tree.build.massex(N, 1, par = params, age = age, extinction_date = extinction_date, survivorship = survivorship)
    } else if (treebuilder == "shift_mass_extinction"){
      t <- tree.build.massex(N, 1, par = params, age = age, extinction_date = extinction_date, survivorship = survivorship)
    } else if (treebuilder == "diversity_dependent"){
      t <- tree.build.dd(N, 1, par = params, age = age, ddmodel = ddmodel)
    } else if (treebuilder == "diversity_dependent_KI") {
      t <- tree.build.ddKI(N, par = params, age = age, ddmodel = ddmodel, KM = KM, tinn = tinn)
    } else if (treebuilder == "diversity_dependent_KI_shiftsp"){
      t <- tree.build.ddKI.spshift(N, par = params, age = age, ddmodel = ddmodel, KM = KM, tinn = tinn)
    } else if (treebuilder == "sym_age_dependent") {
      t <- tree.build.AD.symmetric(N, 1, par = params, shape = shape)[[1]]
    } else if (treebuilder == "asym_age_dependent") {
      t <- tree.build.AD.asymmetric(N, 1, par = params, shape = shape)[[1]] 
    } else if (treebuilder == "sym_uni_shiftsp"){
      t <- tree.build.uni_sp_shift(N, 1, par = params, age = age, shift_date = tinn)
    }
  }
  } else  {
  while(length(getExtinct(t, tol = 0.001)) < e){
    if(treebuilder == "sym_con"){
      t <- tree.build.symmetric(N, 1, par =  params, m = m, gsa = TRUE)[[1]] 
    } else if(treebuilder == "asym_con") {
      t <- tree.build.asymmetric(N, 1, par =  params, m = m, gsa = TRUE)[[1]] 
    } else if (treebuilder == "sym_shiftsp"){
      t <- tree.build.shiftsp(N, 1, par = params, m = m, gsa = TRUE)[[1]] 
    } else if (treebuilder == "asym_shiftsp"){
      t <- tree.build.asymmetric_shiftsp(N, 1, par = params, m = m, gsa = TRUE)[[1]] 
    } else if (treebuilder == "sym_shiftex"){
      t <- tree.build.shiftex(N, 1, par = params, m = m, gsa = TRUE)[[1]]
    } else if (treebuilder == "asym_shiftex") {
      t <- tree.build.asymmetric_shiftex(N, 1, par = params, m = m, gsa = TRUE)[[1]]
    } else if (treebuilder == "sym_shiftsp_shiftex"){
      t <- tree.build.shiftsp_shiftext(N, 1, par = params, m = m, gsa = TRUE)[[1]] 
    } else if (treebuilder == "asym_shiftsp_shiftex"){
      t <- tree.build.asymmetric_shiftsp_shiftext(N, 1, par = params, m = m, gsa = TRUE)[[1]] 
    } else if (treebuilder == "sym_age_dependent") {
      t <- tree.build.AD.symmetric(N, 1, par = params, m = m, gsa = TRUE, shape = shape)[[1]]
    } else if (treebuilder == "asym_age_dependent") {
      t <- tree.build.AD.asymmetric(N, 1, par = params, m = m, gsa = TRUE, shape = shape)[[1]] 
    } else if (treebuilder == "sym_uni_shiftsp"){
      t <- tree.build.uni_sp_shift(N, 1, par = params, age = age, shift_date = tinn)
    }
  }
  }  
  foss_no <- length(getExtinct(t, tol = 0.001))
  
  if(is.null(sampling_fraction) == F){
    extant <- getExtant(t, tol = 0.001)
    N <- N*sampling_fraction
    t <- keep.tip(t, c(sample(extant, N), getExtinct(t, tol = 0.001)))
  }
  
  
  if(e_is_prop == T){
  e <- round(e*foss_no)
  } else {e <- e}
  #you could put in another if else statement here for a version were it applies all sampling methods to same tree
  if(fossilsampler == "random"){
    t <- fossil.random.sampling(t, e)
    t$random <- TRUE
  } else if (fossilsampler == "diversified"){
    t <- fossil.pd.sampling(t, e)  
  } else if (fossilsampler == "time_dependent"){
    t <- fossil.time.sampling(t,e)
  } else if (fossilsampler == "trait_dependent"){
    t <- fossil.bm.sampling(t, e)  
  } else if (fossilsampler == "none") {
    t <- t
    print("no fossilsampler specified, all extint lineages retained")
  } else if (fossilsampler == "all"){
    t <- c(fossil.random.sampling(t, e), fossil.pd.sampling(t,e), fossil.time.sampling(t,e),  fossil.bm.sampling(t, e), t)
  }
  
  if(fossilsampler == 'all'){
    output <- c(ci.norm(t[[1]]), ci.norm(t[[2]]), ci.norm(t[[3]]),ci.norm(t[[4]]),ci.norm(t[[5]]),
                si.norm(t[[1]]), si.norm(t[[2]]), si.norm(t[[3]]),si.norm(t[[4]]),si.norm(t[[5]]),
                gammaStat(t[[1]]), gammaStat(t[[2]]), gammaStat(t[[3]]), gammaStat(t[[4]]), gammaStat(t[[5]]),
                foss_no, e, tree.length(t[[5]]))
    names <- c("CI_random", "CI_diversified", "CI_time_dependent", "CI_trait_dependent", "CI_all",
                       "SI_random", "SI_diversifed", "SI_time_dependent", "SI_trait_dependent", "SI_all",
                       "Gamma_random", "Gamma_diversified", "Gamma_time_dependent", "Gamma_trait_dependent", "Gamma_all", 
                       "Fossil_no", "Preservation_no", "Tree_length")
    if(treebuilder == "con_mass_extinction" | treebuilder == "shift_mass_extinction"){
      output <- c(output, extinction_date, survivorship)
      names(output) <- c(names, "extinction_date", "survivorship")
   # } else if (treebuilder == "diversity_dependent"){
  #  output <- unlist(c(output, t[[5]][5]))
  #  names(output) <- c(names, names(t[[5]][5]))
    # } else if (treebuilder == "diversity_dependent_KI"){
    # output <- unlist(c(output, t[[5]][5],t[[5]][6], t[[5]][7]))##########
    # names(output) <- c(names, names(t[[5]][5]),names(t[[5]][6]), names(t[[5]][7]))
    # } else if (treebuilder == "diversity_dependent_KI_shiftsp"){
      # output <- c(output, t[[5]][5],t[[5]][6], t[[5]][7], t[[5]][8])
      # names(output) <- c(names, names(t[[5]][5]),names(t[[5]][6]), names(t[[5]][7]),  names(t[[5]][8]))
    } else if (treebuilder == "sym_uni_shiftsp"){
      output <- c(output, t[[5]][5],t[[5]][6], t[[5]][7])
      names(output) <- c(names, names(t[[5]][5]),names(t[[5]][6]), names(t[[5]][7]))
    } else {
      names(output) <- names
    }
  } else {  
  if(treebuilder == "con_mass_extinction" | treebuilder == "shift_mass_extinction"){
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), extinction_date, survivorship)
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length", "extinction_date", "survivorship")
  } else if (treebuilder == "diversity_dependent"){
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), t[[5]])
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length", names(t[5]))
  } else if (treebuilder == "diversity_dependent_KI") {
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), t[[5]], t[[6]], t[[7]])
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length",  names(t[5]), names(t[6]), names(t[7]))
  } else if (treebuilder == "diversity_dependent_KI_shiftsp") {
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), t[[5]], t[[6]], t[[7]], t[[8]])
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length",  names(t[5]), names(t[6]), names(t[7]), names(t[8]))
  } else if (treebuilder == "sym_uni_shiftsp") {
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), t[[5]], t[[6]], t[[7]])
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length",  names(t[5]), names(t[6]), names(t[7]))
  }  else {
    output <- c(ci.norm(t), si.norm(t), gammaStat(t), foss_no, e, tree.length(t), t[[12]])
    names(output) <- c("CI", "SI", "Gamma", "Fossil_no", "Preservation_no", "Tree_length", names(t[12]))
  }
  }
  
  if(is.null(sampling_fraction) == F){
    output <- setNames(c(output, sampling_fraction), nm = c(names(output), 'sampling_fraction'))
  }
  
  return(output)
}    
#now you need to add gsa to all of the tree builders

#helper function for adding model names to output of parallel experiments
add_model <- function(df, model){
  df <- as.tibble(df)
  df$model <- paste(model)
  return(df)
}


#below is setting up a parallel run of simulated data under symmertic with speciation shifts and random fossil sampling
sym.ci.sim <- function(params){
  get_sum_stats(params = params, N = 100, e = 10, gsa = TRUE, m = 500,  treebuilder = "sym_con", fossilsampler = "diversified")
}


add_starting_params <- function(ep, params){
  names <- c(names(ep), "lambda_init", "mu_init", "shift_no")
  ep <- c(ep, params[[1]], params[[2]], params[[3]])
  names(ep) <- names
  return(ep)
}

#system.time(stats <- mclapply(X = split(parameters, seq(nrow(parameters))), FUN = sym.ci.sim, mc.cores = detectCores()))

#you need to specify the number of extant (N) and extinct (e) tips, the tree builder and fossil sampler in a secondary function for paralellisation
#sym.ci.shift.macro <- function(params){
#  get_sum_stats(params = params, N = 60, e = 9, treebuilder = "sym_shiftsp", fossilsampler = "random")
#}

#example- run for a single set of parameter values
#tt <- get_sum_stats(params = parameters[1,], treebuilder = "sym_shiftsp", fossilsampler = "random")

#run single core full dataframe of parameter values
#system.time(stats <- apply(X = parameters, FUN = sym.ci.shift.macro, MARGIN = 1))

#multicore full table of parameter values
#system.time(stats <- mclapply(X = split(parameters, seq(nrow(parameters))), 
#                              FUN = sym.ci.shift.macro, mc.cores = detectCores()))

#quickly convert list of lists to dataframe
#stats <- do.call(rbind, stats)
#colnames(stats) <- c("CI", "Gamma")

####################
### TreeBuilders ###
####################
### constant_rates ###
tree.build <- function(n, nsim, par) {sim.bd.taxa(n = n, numbsim = nsim, lambda = par[1], mu = par[2])}


### varying rates - symmetric ###
tree.build.symmetric <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = T, gsa = FALSE)
  } else {
  sim.taxa(n = n, m = m, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = T, gsa = TRUE)
  }
}


### varying rates - asymmetric ###
tree.build.asymmetric <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  #}
  if(gsa == FALSE){
    sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = F)
  } else {
    sim.taxa(n = n, m = m, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = F, gsa = TRUE)
  }
}


### shifting speciation rates ###
tree.build.shiftsp <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftsp <- list(prob = par[4], strength = paste0("runif(", par[5], ",", par[6], ")"))
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, tiplabel = c("sp.", "ext.","", ""))
  } else {
  sim.taxa(n = n, m = m, gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, tiplabel = c("sp.", "ext.","", ""))
  }
}


### shifting extinction rates ###
tree.build.shiftex <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")  #I dont love this because it makes the speciation rate increadibly small and unrealistic for a fixed system... I think I should take the midpoint scaler and multiply by initial rate
  #waitsp <- paste0("rexp(", mean(c(par[[5]], par[[6]]))*par[[1]], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftex <- list(prob = par[4], strength = paste0("runif(", par[7], ",", par[8], ")"))
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftex = shiftex, tiplabel = c("sp.", "ext.","", ""))
  } else {
    sim.taxa(n = n, m = m, gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftex = shiftex, tiplabel = c("sp.", "ext.","", ""))
  }
}


### asymmetric shifting speciation rates ###
tree.build.asymmetric_shiftsp <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftsp <- list(prob = par[4], strength = paste0("runif(", par[5], ",", par[6], ")"))
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, tiplabel = c("sp.", "ext.","", ""), symmetric = F)
  } else {
  sim.taxa(n = n, m = m,  gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, tiplabel = c("sp.", "ext.","", ""), symmetric = F)
  }
}

### assymetic shifting extinction rates ###
tree.build.asymmetric_shiftex <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")  #I dont love this because it makes the speciation rate increadibly small and unrealistic for a fixed system... I think I should take the midpoint scaler and multiply by initial rate
  #waitsp <- paste0("rexp(", mean(c(par[[5]], par[[6]]))*par[[1]], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftex <- list(prob = par[4], strength = paste0("runif(", par[7], ",", par[8], ")"))
  #}
  if(gsa == FALSE){
    sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftex = shiftex, tiplabel = c("sp.", "ext.","", ""), symmetric = F)
  } else {
    sim.taxa(n = n, m = m, gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftex = shiftex, tiplabel = c("sp.", "ext.","", ""))
  }
}


### shifting speciation and extinction rates ###
tree.build.shiftsp_shiftext <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftsp <- list(prob = par[4]/2, strength = paste0("runif(", par[5], ",", par[6], ")")) #divide probability of shift by 2 because extinction and speciation shifts independnet
  shiftext <- list(prob = par[4]/2, strength = paste0("runif(", par[7], ",", par[8], ")"))
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, shiftext = shiftext, tiplabel = c("sp.", "ext.","", ""))
  } else {
  sim.taxa(n = n, m = m, gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, shiftext = shiftext, tiplabel = c("sp.", "ext.","", ""))
  }
}


### asymmetric shifting speciation and extinction rates ###
tree.build.asymmetric_shiftsp_shiftext <- function(n, nsim, par, gsa = FALSE, m){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rexp(", par[1], ")")
  waitext <-  paste0("rexp(", par[2], ")")
  shiftsp <- list(prob = par[4]/2, strength = paste0("runif(", par[5], ",", par[6], ")")) #divide probability of shift by 2 because extinction and speciation shifts independnet
  shiftext <- list(prob = par[4]/2, strength = paste0("runif(", par[7], ",", par[8], ")"))
  #}
  if(gsa == FALSE){
  sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, shiftext = shiftext, tiplabel = c("sp.", "ext.","", ""), symmetric = F)
  } else {
    sim.taxa(n = n, m = m, gsa = TRUE, numbsim = nsim, waitsp = waitsp, waitext = waitext, shiftsp = shiftsp, shiftext = shiftext, tiplabel = c("sp.", "ext.","", ""), symmetric = F)
  }
}


### constant rates but with mass extinction ###
tree.build.massex <- function(n, nsim, par, age, extinction_date =  runif(1, min = 0, max = age), 
                              survivorship = runif(1, min = 0.05, max =0.5)){
  t <- sim.rateshift.taxa(n = n, numbsim = nsim, lambda = c(par[[1]], par[[1]]), mu = c(par[[2]], par[[2]]), 
                          frac = c(1, survivorship), times = c(0, extinction_date), complete = T)[[1]]
  t$extinction_date <- extinction_date
  t$survivorship <- survivorship
  return(t)
}

### mass extinction and subsequent speciation shift ###
tree.build.massex.shift <- function(n, numbsim, par, age, extinction_date =  runif(1, min = 0, max = age), 
                                    survivorship = runif(1, min = 0.05, max =0.5)){
  shift_scaler <- runif(1, par[[5]], par[[6]]) #this is the scaler for speciation rate, we draw a scaler from the range reported from the BAMM analysis
  t <- sim.rateshift.taxa(n = n, numbsim = numbsim, lambda = c(par[[1]]*shift_scaler, par[[1]]), mu = c(par[[2]], par[[2]]), 
                          frac = c(1, survivorship), times = c(0, extinction_date), complete = T)[[1]]
  t$extinction_date <- extinction_date
  t$survivorship <- survivorship
  t$shift_scaler <- shift_scaler
  return(t)
}

### universal speciation shift ###
tree.build.uni_sp_shift <- function(n, numbsim, par, age, shift_date =  runif(1, min = 0, max = age))
  {
  shift_scaler <- runif(1, par[[5]], par[[6]]) #this is the scaler for speciation rate, we draw a scaler from the range reported from the BAMM analysis
  t <- sim.rateshift.taxa(n = n, numbsim = numbsim, lambda = c(par[[1]]*shift_scaler, par[[1]]), mu = c(par[[2]], par[[2]]), 
                          frac = c(1, 1), times = c(0, shift_date), complete = T)[[1]]
  t$shift_date <- shift_date
  t$shift_scaler <- shift_scaler
  return(t)
}


### diversity dependent ###
#here K is assumed to be equal to n
tree.build.dd <- function(n, nsim, par, age, ddmodel){
  pars <- c(par[[1]], par[[2]], n)
  t <- dd_sim(pars = pars, age = age, ddmodel = ddmodel)
  t <- t[[2]]
  return(t)
}

t <- tree.build.dd(20, 1, c(2, .5), 10, ddmodel = 1)

### diversity dependent key innovation *constant* ###
tree.build.ddKI <- function(n, nsim, par, age, ddmodel, KM = sample(seq(round(n*.1), n, 1),1), tinn = runif(1, min = age*0.1, max = age)){
  #KM <- sample(seq(round(n*.1), n, 1),1) #carrying capacity of main clade
  KS <- n - KM#carrying capacity of subclade
  #tinn <- runif(1, min = age*0.1, max = age)
  pars <- c(par[[1]], par[[2]], KM, par[[1]], par[[2]], KS, tinn)
  suppressGraphics(
    {t <- dd_KI_sim(pars = pars, age = age, ddmodel = ddmodel)}
  )
  t <- t[[2]]
  t$n_sub <- KS
  t$tinn <- tinn
  return(t)
}

### diversity dependent key innovation *speciation shift* ###
tree.build.ddKI.spshift <- function(n, nsim, par, age, ddmodel, KM = sample(seq(round(n*.1), n, 1),1), tinn = runif(1, min = age*0.1, max = age)){
  #KM <- sample(seq(round(n*.1), n, 1),1) #carrying capacity of main clade
  shift_scaler <- runif(1, par[[5]], par[[6]])
  KS <- n - KM#carrying capacity of subclade
  #tinn <- runif(1, min = age*0.1, max = age)
  #shift_scaler <- runif(1, min = par[[5]], max= par[[6]])
  pars <- c(par[[1]], par[[2]], KM, (par[[1]]*shift_scaler), par[[2]], KS, tinn)
  suppressGraphics(
    {t <- dd_KI_sim(pars = pars, age = age, ddmodel = ddmodel)}
  )
  t <- t[[2]]
  t$n_sub <- KS
  t$tinn <- tinn
  t$shift_scaler <- shift_scaler
  return(t)
} #test this

### age dependent symmetric
tree.build.AD.symmetric <- function(n, nsim, par, gsa = FALSE, m, shape){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rweibull(", shape[1],", ",  (1/par[1]), ")") #how do you select the shape parameter
  waitext <-  paste0("rweibull(", shape[2],", ",  (1/par[2]), ")")
  #}
  if(gsa == FALSE){
    sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = T, gsa = FALSE)
  } else {
    sim.taxa(n = n, m = m, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = T, gsa = TRUE)
  }
}

### age dependent asymmetric
tree.build.AD.asymmetric <- function(n, nsim, par, gsa = FALSE, m, shape){#, species_wait = 'expo', extinct_wait = 'expo'){
  #if(species_wait == 'expo'){
  waitsp <- paste0("rweibull(", shape[1],", ",  (1/par[1]), ")") #how do you select the shape parameter
  waitext <-  paste0("rweibull(", shape[2],", ",  (1/par[2]), ")")
  #}
  if(gsa == FALSE){
    sim.taxa(n = n, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = F, gsa = FALSE)
  } else {
    sim.taxa(n = n, m = m, numbsim = nsim, waitsp = waitsp, waitext = waitext, symmetric = F, gsa = TRUE)
  }
}


tree.build.ClaDS <- function(n, par){
  t <- sim_ClaDS( lambda_0 = par[[1]], mu_0 = par[[2]],      
                  sigma_lamb=runif(1,0,1),         
                  alpha_lamb=runif(1,0,1),     
                  condition="taxa",    
                  taxa_stop = n,    
                  prune_extinct = F)
  t <- t$tree
}
#plot(ladderize(t))

#ttt <- tree.build.ddKI.spshift(100, 1, par = params[1, ], age = tree.length(tree), ddmodel = 1)
#get_sum_stats(params[1,], n, e, treebuilder = "diversity_dependent_KI_shiftsp", fossilsampler = 'random', gsa = F, ddmodel = 1, KM = 5, tinn = 10, age = 20, shift_scaler = 10)
