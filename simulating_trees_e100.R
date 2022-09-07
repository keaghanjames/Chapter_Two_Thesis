#file for simulating phylogenies using diversification models/preservation models and uninformative priors
setwd("/Users/keaghanyaxley/Dropbox/imbalance_combined_evidence")
source("functions.R")
setwd("BAMM_analyses/simulated_diversification")
require(gtools)
require(parallel)
require(opencpu)



eval_fork <- function(..., timeout = 60) {
  myfork <- parallel::mcparallel({
    eval(...)
  }, silent = FALSE)
  
  # wait max n seconds for a result.
  myresult <- parallel::mccollect(myfork, wait = FALSE, timeout = timeout)
  # kill fork after collect has returned
  tools::pskill(myfork$pid, tools::SIGKILL)
  tools::pskill(-1 * myfork$pid, tools::SIGKILL)
  
  # clean up:
  parallel::mccollect(myfork, wait = FALSE)
  # timeout?
  if (is.null(myresult))
    stop("reached elapsed time limit")
  
  # move this to distinguish between timeout and NULL returns
  myresult <- myresult[[1]]
  
  # send the buffered response
  return(myresult)
}

time_limit <- 120



#tree_size <- c(25, 50, 100, 200, 400, 800)
#tree_size <- c(400, 800)
tree_size <- 200
n_replicates <- 1000
tree_size <- rep(tree_size, n_replicates)
age = 3

#preservation_no <- (0.5)
preservation_no <- 100

  names <- c("CI", "SI","Gamma", "Fossil_no", "Preservation_no", "Tree_length", "random") 




some_params <- function(nsim, size) {
  draw_params(nsim = nsim, hyper_lambda = runif(nsim, 1, 5), hyper_mu = runif(nsim, 0.1, 0.5), mu_as_proportion = T, shift_no = rpois(nsim,2)+1, mu_scalers = c(0.1, 3), lambda_scalers =  c(0.1, 3), N = size)
}

some_params(1, 20)


#constant
sym_con <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_con", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_con", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  }
  ep <- c(ep, params)
  ep$n_extant <- size
  return(unlist(ep))
}

#sym_con(100)

system.time(sym_con <- mclapply(X = tree_size, FUN = sym_con, mc.cores = detectCores()-2))
sym_con <- do.call(rbind, sym_con)
sym_con <- add_model(sym_con, "constant")



#shifting speciation rate
shift_sp <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftsp", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftsp", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
  return(unlist(ep))
  }
system.time(shift_sp <- mclapply(X = tree_size, FUN = shift_sp, mc.cores = detectCores()-2))
shift_sp <- do.call(rbind, shift_sp)
shift_sp <- add_model(shift_sp, "shifting_specation")
df <- bind_rows(sym_con, shift_sp)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))


#shifting_extinction rates
shift_ex <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftex", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftex", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
} 
system.time(shift_ex <- mclapply(X = tree_size, FUN = shift_ex, mc.cores = detectCores()-2))
shift_ex <- do.call(rbind, shift_ex)
shift_ex <- add_model(shift_ex, "shifting_extinction")
df <- bind_rows(df, shift_ex)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))

#shifting speciation and extinction 
shift_ex_sp <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftsp_shiftex", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
     ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_shiftsp_shiftex", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
} 
system.time(shift_ex_sp <- mclapply(X = tree_size, FUN = shift_ex_sp, mc.cores = detectCores()-2))
shift_ex_sp <- do.call(rbind, shift_ex_sp)
shift_ex_sp <- add_model(shift_ex_sp, "shifting_speciation_extinction")
df <- bind_rows(df, shift_ex_sp)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))


#constant except mass extinction
con_mass_extinction <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "con_mass_extinction", fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2)))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "con_mass_extinction", fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2)))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
}
system.time(con_mass_extinction <- mclapply(X = tree_size, FUN = con_mass_extinction, mc.cores = detectCores()-2))
#system.time(con_mass_extinction <- lapply(X = tree_size, con_mass_extinction))

con_mass_extinction <- do.call(rbind, con_mass_extinction)
con_mass_extinction <- add_model(con_mass_extinction, "mass_extinction_constant")
df <- bind_rows(df, con_mass_extinction)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))



#shifting with mass extinction
shift_mass_extinction <- function(size){
  setTimeLimit(time_limit)
    params <- some_params(1, size = size)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "shift_mass_extinction", fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2)))
  while(class(ep) == 'try-error') {
    setTimeLimit(time_limit)
    params <- some_params(1, size = size)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "shift_mass_extinction", fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2)))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
}
system.time(shift_mass_extinction <- mclapply(X = tree_size, FUN = shift_mass_extinction, mc.cores = detectCores()-2))
shift_mass_extinction <- do.call(rbind, shift_mass_extinction)
shift_mass_extinction <- add_model(shift_mass_extinction, "mass_extinction_shifting")
df <- bind_rows(df, shift_mass_extinction)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))

# age-dependent
age_dependent <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_age_dependent", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, shape = c(.7,1)))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "sym_age_dependent", fossilsampler = 'all', age = age, gsa = T, m = size*5, e = preservation_no, ddmodel = 1, e_is_prop = F, shape = c(.7,1)))
  }
  ep <- c(ep, params)
  ep$n_extant <- size
  return(unlist(ep))
}
system.time(age_dependent <- mclapply(X = tree_size, FUN = age_dependent, mc.cores = detectCores()-2))
age_dependent <- do.call(rbind, age_dependent)
age_dependent <- add_model(age_dependent, "age_dependent")
df <- bind_rows(df, age_dependent)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))


## universal speciation shift ###
uni_sp_shift <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, e = preservation_no, treebuilder = "sym_uni_shiftsp", fossilsampler = "all",  age = age, tinn = runif(1, min = 0, max = age), e_is_prop = F))
  while(class(ep) == "try-error"){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, e = preservation_no, treebuilder = "sym_uni_shiftsp", fossilsampler = "all",  age = age, tinn = runif(1, min = 0, max = age), e_is_prop = F))
  }
  ep <- c(ep, params)
  ep$n_extant <- size
  return(unlist(ep))
}

system.time(uni_sp_shift <- mclapply(X = tree_size, FUN = uni_sp_shift, mc.cores = detectCores()-2))
uni_sp_shift <- do.call(rbind, uni_sp_shift)
uni_sp_shift <- add_model(uni_sp_shift, "universal_speciation_shift")
df <- bind_rows(df, uni_sp_shift)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))



#diversity dependent
diversity_dependent <- function(size){
  #setTimeLimit(14400, transient = T)
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = size, treebuilder = "diversity_dependent", fossilsampler = 'all', age = age, gsa = F, e = preservation_no, ddmodel = 1, e_is_prop = F))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "diversity_dependent", fossilsampler = 'all', age = age, gsa = F, e = preservation_no, ddmodel = 1, e_is_prop = F))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
}
#system.time(diversity_dependent <- lapply(tree_size, diversity_dependent))
system.time(diversity_dependent <- mclapply(X = tree_size, FUN = diversity_dependent, mc.cores = detectCores()-2))
diversity_dependent <- do.call(rbind, diversity_dependent)
diversity_dependent <- add_model(diversity_dependent, "diversity_dependent")
df <- bind_rows(df, diversity_dependent)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))


#diversity dependent with key innovation


diversity_dependent_KI <- function(size){
  params <- some_params(1, size = size)
  setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "diversity_dependent_KI", 
                          fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, 
                          ddmodel = 1, e_is_prop = F, survivorship = runif(1, min = 0.2, max = 0.2),
                          KM = round(runif(1, min = size*0.1, max = size*0.5)), 
                          tinn = runif(1, min = age*0.1, max = age)))
  while(class(ep) == 'try-error') {
    params <- some_params(1, size = size)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = size, treebuilder = "diversity_dependent_KI", 
                            fossilsampler = 'all', age = age, gsa = F, m = size*5, e = preservation_no, 
                            ddmodel = 1, e_is_prop = F, survivorship = runif(1, min = 0.2, max = 0.2), 
                            KM = round(runif(1, min = size*0.1, max = size*0.5)), 
                            tinn = runif(1, min = age*0.1, max = age)))
  }
    ep <- c(ep, params)
    ep$n_extant <- size
    return(unlist(ep))
}
system.time(diversity_dependent_KI <- mclapply(X = tree_size, 
                                               FUN = diversity_dependent_KI, 
                                               mc.cores = detectCores()-1))
diversity_dependent_KI <- do.call(rbind, diversity_dependent_KI)
diversity_dependent_KI <- add_model(diversity_dependent_KI, "diversity_dependent_KI")
df <- bind_rows(df, diversity_dependent_KI)
write.csv(df, file = paste0("simulation_summary_simulated_e100.csv"))
save(df, file = paste0("simulation_summary_simulated_e100.Rdata"))





