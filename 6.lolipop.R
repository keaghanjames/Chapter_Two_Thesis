setwd("/Users/keaghanyaxley/Dropbox/imbalance_combined_evidence/")
source("functions.R")

study <- 'Brennan_in_press'
taxa <- 'Macropodinae'

setwd(paste0('BAMM_analyses/', study))
tree <- read.tree("trim.tree")

imbalance <- setNames(ci.norm(tree), nm = taxa)

n <- length(tree$tip.label)
extant <- length(getExtant(tree, tol = 0.001))

setwd("../..")

###################################
  
study <- 'Gustafosn_et_al_2017'
taxa <- 'Gyrinidae'

setwd(paste0('BAMM_analyses/', study))
tree <- read.tree("trim.tree")

n <- c(n, length(tree$tip.label))
extant <- c(extant, length(getExtant(tree, tol = 0.001)))

imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
setwd("../..")

##################################
#study <- 'Herrera_and_davalos'
# taxa <- 'lemur'

study <- 'Koch_et_al_2020'
taxa <- 'Echinoidea'

setwd(paste0('BAMM_analyses/', study))
tree <- read.tree("trim.tree")

n <- c(n, length(tree$tip.label))
extant <- c(extant, length(getExtant(tree, tol = 0.001)))

imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
setwd("../..")

################################

study <- 'Dornburgh_et_al_2015'
taxa <- 'Holocentridae'

setwd(paste0('BAMM_analyses/', study))
tree <- read.tree("trim.tree")

n <- c(n, length(tree$tip.label))
extant <- c(extant, length(getExtant(tree, tol = 0.001)))

imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
setwd("../..")


################################

study <- 'Slater_et_al_2017'
taxa <- 'Mysticeti'

setwd(paste0('BAMM_analyses/', study))
tree <- read.tree("trim.tree")

n <- c(n, length(tree$tip.label))
extant <- c(extant, length(getExtant(tree, tol = 0.001)))

imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
setwd("../..")


# ################################
# 
# study <- 'Paterson_et_al_2020'
# taxa <- 'Pinnipedia'
# 
# setwd(paste0('BAMM_analyses/', study))
# tree <- read.tree("trim.tree")
# 
# n <- c(n, length(tree$tip.label))
# extant <- c(extant, length(getExtant(tree, tol = 0.001)))
# 
# 
# imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
# setwd("../..")
# 
# ###########################
# 
# study <- 'Kealy_and_Beck_2017'
# taxa <- 'Dasyuromorphia'
# 
# setwd(paste0('BAMM_analyses/', study))
# tree <- read.tree("trim.tree")
# 
# n <- c(n, length(tree$tip.label))
# extant <- c(extant, length(getExtant(tree, tol = 0.8)))
# 
# imbalance <- c(imbalance, setNames(ci.norm(tree), nm = taxa))
# setwd("../..")
# 


#imbalance <- sort(imbalance)
imbalance <- tibble(imbalance = imbalance, n = n, extant = extant, clade = names(imbalance))
#imbalance <- imbalance %>% arrange(imbalance)
imbalance$clade <- factor(imbalance$clade, levels = imbalance$clade, ordered = T)
imbalance <- sort(imbalance)
imbalance$size_n <- imbalance$n/8
imbalance$size_e <- imbalance$extant/8
imbalance

require(ggplot2)

p1 <- ggplot(imbalance, aes(x = clade, y = imbalance))+
  geom_segment( aes(x=clade, xend=clade, y=0, yend=imbalance), color="grey") +
  geom_point(color = "grey", size= imbalance$size_n, alpha = .4) +
  geom_point(color = "darkgrey", size= imbalance$size_e, alpha = 1) +
  ylim(0,0.6) +
  theme_light() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.border = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  xlab("") +
  scale_x_discrete(limits = as.character(imbalance$clade)) +
  ylab("Normalised Colless Index") 
p1

ggsave(device = 'pdf', filename = "lollipop.pdf", width = 19, height = 12)

