require(RPANDA)
require(beepr)
require(pbmcapply)
setwd("/Users/keaghanyaxley/Dropbox/imbalance_combined_evidence/")
source("functions.R")

set.seed(888)
counter <- function(start, inc, target){
  x <- start
  while(x < target){
    x <- x + inc
    print(x)
  }
  return(x)
}

beep(2, counter(0, 5, 10000))

#CE_study <- "Brennan_in_press"
#CE_study <- "Dornburgh_et_al_2015"
CE_study <- "Gustafosn_et_al_2017"
#CE_study <- "Kealy_and_Beck_2017"
#CE_study <- "Slater_et_al_2017"
#CE_study <- "Koch_et_al_2020"
#CE_study <- "Paterson_et_al_2020"
#CE_study <- 'Aze_et_al_2011'
require(gtools)
require(fitdistrplus)
require(R.devices)
require(beepr)
setwd(paste0("BAMM_analyses/", CE_study))
tree_file <- 'trim.tree'
#taxa <- "Macropodinae"
#taxa <- "Holocentridae"
#taxa <- "Foraminifera"
#taxa <- "lemur"
#taxa <- "Dasyuromorphia"
#taxa <- 'Mysticeti'
taxa <- 'Gyrinidae'
#taxa <- "echinoidea"
#taxa <- "Pinnipedia"
#taxa <- 'Foraminifera'
dev.off()
#tree <- read.tree(file = "Macropodinae.tree")
tree <- read.tree(file = tree_file)
#tree <- read.tree('Foraminifera.tree')
plot(tree, show.tip.label = F)


#lets get our fixed parameters
age <- tree.length(tree)
n <- length(getExtant(tree, tol = age*0.01))
n
e <- length(tree$tip.label) - n
e
time_limit <- 60
nsim <- 1000
n
e


#sample_fraction <- 56/60 #macropodidae
#sample_fraction <- 0.1 #Gyrinindae
#sample_fraction <- n/33 #pinipedia
sample_fraction <- n/87 #Holocentridae
#sample_fraction <- n/73 #Dasyuromorphia
#sample_fraction <- n/950 #echinoidea
#sample_fraction <- n/16 #mysticetti
sample_fraction

#clad0 <- fit_ClaDS0(tree = tree, iteration = 1000, thin = 10, nCPU = 3)
#beep(3, clad <- fit_ClaDS(tree, iterations = 100000, thin = 100, it_save = 100, file_name = 'clad', nCPU = 3, sample_fraction = sample_fraction))

#save(clad, file = "cladDS.R")
load('cladDS.R')

MAPS <- getMAPS_ClaDS(clad)
lambda <- MAPS[4]
mu <- MAPS[3]*lambda
plot_ClaDS_chains(clad)
plot_ClaDS_phylo(tree, MAPS[-(1:4)])

params <- draw_params(1, N = n, shift_no = 10)
params[1] <- lambda
params[2] <- mu
scaler <- range(MAPS[-c(1:4)]/lambda)
params[5:8] <- c(scaler, scaler)
params
#system.time(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_shiftsp", fossilsampler = "all", sampling_fraction = sample_fraction))
time_limit <- 20

m <- n*5

df <- as.tibble(read.csv("ClaDS_summary.csv"))
summary(as.factor(df$model))
### Shifting Speciation ####

runner <- function(run){
  params[[3]] <- sample(seq(1,10,1), 1)
  params[[4]] <- params[[3]]/(n+e-1)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = (n/sample_fraction)*5, gsa = F, e = e, treebuilder = "sym_shiftsp", fossilsampler = "all", sampling_fraction = sample_fraction))
   while(class(ep) == 'try-error'){
     params[[3]] <- sample(seq(1,10,1), 1)
     params[[4]] <- params[[3]]/(n+e-1)
     setTimeLimit(time_limit)
     ep <- try(get_sum_stats(params = params, N = n, m = (n/sample_fraction)*5, gsa = F, e = e, treebuilder = "sym_shiftsp", fossilsampler = "all", sampling_fraction = sample_fraction))
   }
  ep <- add_starting_params(ep, params)
  return(ep)
}
system.time(runner(1))
system.time(sym_shiftsp <- pbmclapply(X = seq(1,nsim,1), FUN = runner, mc.preschedule = T , mc.cores = detectCores()-2))
sym_shiftsp[grep("Error", sym_shiftsp)] <- NULL
sym_shiftsp <- as.data.frame(split(unlist(sym_shiftsp), names(unlist(sym_shiftsp))))
sym_shiftsp <- add_model(sym_shiftsp, "shifting_sp")
#df <- sym_shiftsp
df <- bind_rows(df, sym_shiftsp)
df

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)


### Shifting Extinction ###

runner <- function(run){
  params[[3]] <- sample(seq(1,10,1), 1)
  params[[4]] <- params[[3]]/(n+e-1)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_shiftex", fossilsampler = "all", sampling_fraction = sample_fraction))
  
   while(class(ep) == 'try-error'){
     params[[3]] <- sample(seq(1,10,1), 1)
     params[[4]] <- params[[3]]/(n+e-1)
     setTimeLimit(time_limit)
     ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_shiftex", fossilsampler = "all", sampling_fraction = sample_fraction))
   }
  ep <- add_starting_params(ep, params)
  return(ep)
}
#runner(1)

system.time(sym_shiftex <- pbmclapply(X = seq(1,nsim,1), FUN = runner, mc.preschedule = T, mc.cores = detectCores()-2))
sym_shiftex[grep("Error", sym_shiftex)] <- NULL
sym_shiftex <- as.data.frame(split(unlist(sym_shiftex), names(unlist(sym_shiftex))))
sym_shiftex <- add_model(sym_shiftex, "shifting_ex")
df <- bind_rows(df, sym_shiftex)
df

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
summary(as.factor(df$model))

### Shifting Speciation and Extinction ###

runner <- function(run){
  params[[3]] <- sample(seq(1,10,1), 1)
  params[[4]] <- params[[3]]/(n+e-1)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_shiftsp_shiftex", fossilsampler = "all", sampling_fraction = sample_fraction))
  while(class(ep) == 'try-error'){
    params[[3]] <- sample(seq(1,10,1), 1)
    params[[4]] <- params[[3]]/(n+e-1)
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_shiftsp_shiftex", fossilsampler = "all", sampling_fraction = sample_fraction))
  }
  ep <- add_starting_params(ep, params)
  return(ep)
}

system.time(sym_shiftsp_shiftex <- pbmclapply(X = seq(1, nsim,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
sym_shiftsp_shiftex[grep("Error", sym_shiftsp_shiftex)] <- NULL
sym_shiftsp_shiftex <- as.data.frame(split(unlist(sym_shiftsp_shiftex), names(unlist(sym_shiftsp_shiftex))))
sym_shiftsp_shiftex <- add_model(sym_shiftsp_shiftex, "shifting_sp_ex")
df <- bind_rows(df, sym_shiftsp_shiftex)
summary(as.factor(df$model))

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)


### Constant w Mass Extinction ###

runner <- function(run){
  ep <- try(get_sum_stats(params = params, N = n, gsa = FALSE, e = e, treebuilder = "con_mass_extinction", fossilsampler = "all",
                          age = age, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2), sampling_fraction = sample_fraction))
  ep <- add_starting_params(ep, params)
  return(ep)
}
system.time(con_mass_extinction <- pbmclapply(X = seq(1,50,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
con_mass_extinction[grep("Error", con_mass_extinction)] <- NULL
con_mass_extinction <- as.data.frame(split(unlist(con_mass_extinction), names(unlist(con_mass_extinction))))
con_mass_extinction <- add_model(con_mass_extinction, "con_ME")
df <- bind_rows(df, con_mass_extinction)

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
summary(as.factor(df$model))


### Age Dependent ###

runner <- function(run){
 setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_age_dependent", fossilsampler = "all", shape = c(.7,1), sampling_fraction = sample_fraction))
  while(class(ep) == 'try-error'){
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_age_dependent", fossilsampler = "all", shape = c(.7,1), sampling_fraction = sample_fraction))
  }
  ep <- add_starting_params(ep, params)
  return(ep)
}
#runner(1)
system.time(sym_age_dependent <- pbmclapply(X = seq(1,nsim,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
sym_age_dependent[grep("Error", sym_age_dependent)] <- NULL
sym_age_dependent <- as.data.frame(split(unlist(sym_age_dependent), names(unlist(sym_age_dependent))))
sym_age_dependent <- add_model(sym_age_dependent, "age_dependent")
# df <- bind_rows(df, sym_age_dependent)
# 
# write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
# df

# ### Diversity Dependent ###
# 
 runner <- function(run){
   setTimeLimit(time_limit)
   ep <- try(get_sum_stats(params = params, N = n, gsa = FALSE, e = e, treebuilder = "diversity_dependent", fossilsampler = "all", age = age*2, ddmodel = 1, sampling_fraction = sample_fraction))
   # while(class(ep) == 'try-error'){
   #   setTimeLimit(60)
   #   ep <- try(get_sum_stats(params = params, N = n, gsa = FALSE, e = e, treebuilder = "diversity_dependent", fossilsampler = "all", age = age, ddmodel = 1, sampling_fraction = sample_fraction))
   # }
   ep <- add_starting_params(ep, params)
   return(ep)
 }
 system.time(runner(1))
 
 # sim_DD <- list()
 # for(i in 1:20){
 # #  setTimeLimit(60)
 #   sim_DD[[i]] <- runner(1)
 #   print(i)
 #  }
 # 
 system.time(sim_DD <- pbmclapply(X = seq(1,600,1), FUN = runner, mc.preschedule = T, mc.cores = detectCores()-1))
 sim_DD[grep("Error", sim_DD)] <- NULL
 sim_DD <- as.data.frame(split(unlist(sim_DD), names(unlist(sim_DD))))
 sim_DD <- add_model(sim_DD, "diversity dependent")
 df <- bind_rows(df, sim_DD)
 df$model <- as.factor(df$model)
 summary(df$model)
 
 write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
 
# 
### Diversity Dependent w KI ###
 
 runner <- function(run){
   setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = n, gsa = FALSE, e = e, treebuilder = "diversity_dependent_KI", fossilsampler = "all", age = age*2, ddmodel = 1,
                           KM = round(runif(1, min = (n/sample_fraction)*0.1, max = (n/sample_fraction)*0.5)), tinn = runif(1, min = 1, max = age*2), sampling_fraction = sample_fraction))
   # while(class(ep) == 'try-error'){
   #   ep <- try(get_sum_stats(params = params, N = n, gsa = FALSE, e = e, treebuilder = "diversity_dependent_KI", fossilsampler = "all", age = age, ddmodel = 1,
   #                           KM = round(runif(1, min = (n/sample_fraction)*0.1, max = (n/sample_fraction)*0.5)), tinn = runif(1, min = 1, max = age), sampling_fraction = sample_fraction))
   # }
   ep <- add_starting_params(ep, params)
   return(ep)
 }
 #runner(1)
 system.time(DD_KI <- pbmclapply(X = seq(1,2000,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
 DD_KI[grep("Error", DD_KI)] <- NULL
 DD_KI <- as.data.frame(split(unlist(DD_KI), names(unlist(DD_KI))))
 DD_KI <- add_model(DD_KI, "diversity dependent w KI")
 df <- bind_rows(df, DD_KI)
 
 write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
 summary(as.factor(df$model))
# 

### Shifting w Mass Extinction ###

runner <- function(run){
  params[[3]] <- sample(seq(1,10,1), 1)
  params[[4]] <- params[[3]]/(n+e-1)
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = FALSE, e = e, treebuilder = "shift_mass_extinction", fossilsampler = "all",
                          age = age, extinction_date =  runif(1, min = 0, max = age), survivorship = runif(1, min = 0.2, max = 0.2), sampling_fraction = sample_fraction))
  ep <- add_starting_params(ep, params)
  return(ep)
}
runner(1)
system.time(shift_mass_extinction <- pbmclapply(X = seq(1,100,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
shift_mass_extinction[grep("Error", shift_mass_extinction)] <- NULL
shift_mass_extinction <- as.data.frame(split(unlist(shift_mass_extinction), names(unlist(shift_mass_extinction))))
shift_mass_extinction <- add_model(shift_mass_extinction, "shiftsp_ME")
df <- bind_rows(df, shift_mass_extinction)

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
summary(as.factor(df$model))

### Universal Speciation Shift ###

runner <- function(run){
  setTimeLimit(20)
  ep <- try(get_sum_stats(params = params, N = n, e = e, treebuilder = "sym_uni_shiftsp", fossilsampler = "all",  age = age, tinn = runif(1, min = 0, max = age), sampling_fraction = sample_fraction))
  ep <- add_starting_params(ep, params)
  return(ep)
}
runner(1)
system.time(uni_shift_sp <- pbmclapply(X = seq(1,150,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
uni_shift_sp[grep("Error", uni_shift_sp)] <- NULL
uni_shift_sp <- as.data.frame(split(unlist(uni_shift_sp), names(unlist(uni_shift_sp))))
uni_shift_sp <- add_model(uni_shift_sp, "universal_sp_shift")
df <- bind_rows(df, uni_shift_sp)

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
beepr(print('hey'))

### Constant ###

runner <- function(run){
  setTimeLimit(time_limit)
  ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_con", fossilsampler = "all", sampling_fraction = sample_fraction))
  while(class(ep) == 'try-error'){
    setTimeLimit(time_limit)
    ep <- try(get_sum_stats(params = params, N = n, m = n*5, gsa = F, e = e, treebuilder = "sym_con", fossilsampler = "all", sampling_fraction = sample_fraction))
  }
  ep <- add_starting_params(ep, params)
  return(ep)
}

system.time(x) <- runner(1)
system.time(sym_con <- pbmclapply(X = seq(1,nsim,1), FUN = runner, mc.preschedule = F, mc.cores = detectCores()-1))
sym_con[grep("Error", sym_con)] <- NULL
sym_con <- as.data.frame(split(unlist(sym_con), names(unlist(sym_con))))
sym_con <- add_model(sym_con, "con")
df <- bind_rows(df, sym_con)

df$model <- as.factor(df$model)
summary(df$model)

write.csv(df, file = 'ClaDS_summary.csv', row.names = F)
save(df, file = 'ClaDS_summary.Rdata')
##### GO BACK AND DO DD ### 
