require(tidyverse)
require(broom)
require(lme4)
require(emmeans)
require(dplyr)
require(usdm)
require(knitr)
require(reshape2)
require(fastDummies)
require(lmtest)
require(hexbin)
require(lmerTest)

setwd("~/Dropbox/imbalance_combined_evidence/BAMM_analyses/simulated_diversification/")

############################################################
### Importing, Formating and Subsampling Simulation Data ###
############################################################

df <- as_tibble(read.csv("CI_simulation_output_ne100.csv"))

#now we need to make sure our variables are coded correctly 

#we add a tree ID
df$tree_ID <- as.factor(paste0("tree_",rep(seq(1, (nrow(df)/4), 1), 4)))

#make sure our categorical data is factor
df$preservation <- as.factor(df$preservation)
df$model <- factor(df$model, levels = c("constant", "universal_speciation_shift", "shifting_specation", "shifting_extinction",
                                        "shifting_speciation_extinction", "mass_extinction_constant",
                                        "mass_extinction_shifting", "age_dependent", "diversity_dependent", "diversity_dependent_KI"))
levels(df$model) <- str_to_title(str_replace_all(levels(df$model), "_", " "))

df$Preservation_no <- as.factor(df$Preservation_no)
#calculate our weights
df$weights <- df$Fossil_no+df$n_extant

df$preservation <- factor(df$preservation, levels = c("random", "diversified", "trait_dependent", "time_dependent"))

#creat a tidy object so we can run correlations between treatments
by_model <-  group_by(df[which(df$preservation == 'random'),], model)


#lets look at how our initial lambda and mu values are distributed 
#lets look more closely at those relationships
ggplot(df, aes(x = lambda_init, y = mu_init, color = model))+
  geom_point(alpha=.1) +
  geom_smooth(method = 'lm') +
  theme_classic() +
  xlab("Initial speciation rate") +
  ylab("Initial extinction rate")



p0 <- ggplot(df, aes(x = index, y = model, fill = preservation))+
  stat_summary(geom = "boxplot", 
               fun.data = function(x) setNames(quantile(x, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("ymin", "lower", "middle", "upper", "ymax")), 
               position = "dodge")  + 
  xlab('Normalised Colless Index') +
  theme_classic() +
  theme(legend.title = element_blank())+
  scale_fill_manual(labels = c('Random', 'Diversified', 'Trait Dependent', 'Time Dependent'), values = c("#264653", "#2a9d8f", "#f4a261", "#e76f51")) +
  ylab("")
p0
ggsave("figures/colless_by_model.pdf", device = 'pdf', scale = 0.7)

p1 <- ggplot(df[which(df$preservation == 'random'),], aes(lambda_init, y = mu_init))+
      geom_hex(bins = 20) +
#      geom_smooth(method = 'lm', color = 'black', lwd = 0.5) +
#      stat_density2d(color = 'black', alpha = 0.3, linejoin = 'bevel') +
 # stat_smooth(method = 'lm', se = T, fullrange = T) +
  
      theme_classic()+ 
      xlab('Initial speciation rate')+
      ylab('Initial extinction rate')+
      scale_fill_distiller(palette=8, direction=1) +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      facet_wrap( ~ model, nrow = 5)
p1

ggsave("figures/initial_lamb_mu.pdf", device = 'pdf', scale = .7)
save(p1, file = 'initial_lamb_mu.Rdata')


results <- do(by_model, tidy(cor.test(.$lambda_init, .$mu_init)))
filter(results, estimate >= 0.1)
write.csv(results, file = 'cor_lambda_mu.csv')


starting_params <- tibble(score = c(df[which(df$preservation == 'random'),]$lambda_init, df[which(df$preservation == 'random'),]$mu_init),
       index = c(rep("lambda", length(which(df$preservation == 'random'))), rep("mu", length(which(df$preservation == 'random')))),
       model = rep((df$model[which(df$preservation == 'random')]),2),
       fossil_no = rep(df$Fossil_no[which(df$preservation == 'random')], 2),
       tree_length = rep(df$Tree_length[which(df$preservation == 'random')], 2),
       imbalance = rep(df$index[which(df$preservation == 'random')],2)
       )


p2 <- ggplot(starting_params, aes(x = score, y = log(fossil_no), color = index))+
  geom_point(alpha = 0.3)+
  facet_wrap( ~ model, nrow = 5) +
  theme_classic() +
  xlab('Initial rate') +
  ylab("ln(Fossil number)") +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  labs(color = 'Parameter')
p2
ggsave("figures/initial_paramters_fossil_no.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test(.$lambda_init, .$Fossil_no)))
filter(results, estimate >= 0.1)
write.csv(results, file = 'cor_lambda_fossilno.csv')

results <- do(by_model, tidy(cor.test(.$mu_init, .$Fossil_no)))
filter(results, estimate >= 0.1)
write.csv(results, file = 'cor_mu_fossilno.csv')

p3 <- ggplot(starting_params, aes(x = score, y = log(tree_length), color = index))+
  geom_point(alpha = 0.1)+
  facet_wrap( ~ model, nrow = 5) +
  theme_classic() +
  xlab('Initial rate') +
  ylab("ln(Tree length)") +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  labs(color = 'Parameter')
p3
ggsave("figures/initial_paramters_tree_length.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test(.$lambda_init, .$Tree_length)))
write.csv(results, file = 'cor_lambda_treelength.csv')

results <- do(by_model, tidy(cor.test(.$mu_init, .$Tree_length)))
write.csv(results, file = 'cor_mu_treelength.csv')


p4 <- ggplot(starting_params, aes(x = score, y = imbalance, color = index))+
  geom_point(alpha = 0.1)+
  facet_wrap( ~ model, nrow = 5) +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  theme_classic() +
  xlab('Initial rate') +
  ylab("Normalised Colless index")+
  labs(color = 'Parameter')
p4
ggsave("figures/initial_paramters_normailised_colless_index.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test(.$lambda_init, .$index)))
write.csv(results, file = 'cor_lambda_colless.csv')

results <- do(by_model, tidy(cor.test(.$mu_init, .$index)))
write.csv(results, file = 'cor_mu_coless.csv')



p5 <- ggplot(df[which(df$preservation == 'random'),], aes(x = index, y = lambda_init-mu_init)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  facet_wrap(~ model, nrow = 5) +
  theme_classic() +
  xlab('Normalised Colless Index') +
  ylab('Initial diversification rate')
p5
ggsave("figures/initial_diversification_normailised_colless_index.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test((.$lambda_init - .$mu_init), .$index)))
results
write.csv(results, file = 'cor_diversification_coless.csv')

p6 <- ggplot(df[which(df$preservation == 'random'),], aes(x = index, y = mu_init/lambda_init)) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  facet_wrap(~ model, nrow = 5) +
  theme_classic() +
  xlab('Normalised Colless Index') +
  ylab('Initial turnover rate')
p6
ggsave("figures/initial_turnover_normailised_colless_index.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test((.$mu_init/.$lambda_init), .$index)))
results
write.csv(results, file = 'cor_turnover_colless.csv')


p7 <- ggplot(df[which(df$preservation == 'random'),], aes(x = index, y = log(Fossil_no))) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  facet_wrap(~ model, nrow = 5) +
  theme_classic() +
  xlab('Normalised Colless Index') +
  ylab('ln(Fossil number)')
p7
ggsave("figures/initial_log_fossil_no_normailised_colless_index.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test((.$Fossil_no), .$index)))
results
write.csv(results, file = 'cor_fossilno_colless.csv')


p8 <- ggplot(df[which(df$preservation == 'random'),], aes(x = index, y = log(Tree_length))) +
  geom_point(alpha = 0.1) +
  stat_smooth(method = 'lm', se = F, fullrange = T) +
  facet_wrap(~ model, nrow = 5) +
  theme_classic() +
  xlab('Normalised Colless Index') +
  ylab('ln(Tree length)')
p8
ggsave("figures/initial_log_tree_length_normailised_colless_index.pdf", device = 'pdf', scale = .7)

results <- do(by_model, tidy(cor.test((.$Tree_length), .$index)))
results
write.csv(results, file = 'cor_treelength_colless.csv')

#and scale continuous variables 
df$Fossil_no <- scale(df$Fossil_no)
df$Tree_length <- scale(df$Tree_length)
df$mu_init <-scale(df$mu_init)
df$lambda_init <-scale(df$lambda_init)

#we are now going to reduce our dataset, at present it contains 40,000 observations whcih is too large for this computer
#we will reduce it to 4,000 wile ensuring equal representation among treatments

# models <- unique(df$model)
# preservation <- unique(df$preservation)

# new_df <- tibble()
# for(i in 1:length(models)){
#   for(j in 1:length(preservation)){
#     new_df <- rbind(new_df, df[which(df$model == models[i] & df$preservation == preservation[j]),][1:100,])
#   }
# }
# unique(new_df$preservation)
# 
# df <- new_df
# df

####################
### Co-linearity ###
####################

#before we start building our models we want to check for co-linearity between our continous predictor variables
colin <- as.data.frame(df[,1:5])
head(colin)
colin$Preservation_no <- NULL
colin$model <- df$model
#we are also going to for now ignore our different preservation treatments and subsample only those with random fossil records
colin$preservation <- df$preservation
colin$turnover_init <- df$mu_init/df$lambda_init
colin$diversification_init <- df$lambda_init - df$mu_init
colin <- colin[which(colin$preservation == 'random'),]
colin['preservation'] <- NULL
head(colin)

get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

cor_by_model <- function(colin, name){
  m <- colin[which(colin$model == name),]
  m$model <- NULL
  cor <- round(cor(m), 2)
  cor <- get_upper_tri(cor)
  cor <- melt(cor, na.rm = T)
  ggplot(data = cor, aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    ggtitle(name)+
    scale_fill_gradient2(low = colours[1], high = colours[4], mid = colours[2], 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 5, hjust = 1), 
          axis.text.y = element_text(angle = 45, vjust = 1, 
                                     size = 5, hjust = 1))+
    coord_fixed()+ 
    geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
    scale_x_discrete(labels = c("N extinct", "Age", expression(paste(lambda)), expression(paste(mu)), expression(paste(mu,"/",lambda)),
                                expression(paste(lambda,"-",mu))))+
    scale_y_discrete(labels = c("N extinct", "Age", expression(paste(lambda)), expression(paste(mu)), expression(paste(mu,"/",lambda)),
                                expression(paste(lambda,"-",mu))))+
        theme(
          plot.title = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = 'none',
      legend.direction = "horizontal")+
    guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                                 title.position = "top", title.hjust = 0.5))
  
}

cor_by_model_simple <- function(name) {cor_by_model(colin = colin, name = name)}

models <- levels(colin$model)

cor_tables <- lapply(models, cor_by_model_simple)

#cor_tables
n <- length(cor_tables)
nCol <- floor(sqrt(n))
grid <- do.call("grid.arrange", c(cor_tables, ncol=2))
ggsave(plot = grid, "figures/correlation_heatmaps.pdf", device = 'pdf', scale = .7)
              
              
grid.arrange(unlist(cor_tables), nrow = 5)

cor_by_model(colin, name = 'Constant')



#we need to create dummr variables for the tree buidling models
colin <- dummy_cols(colin)

colin <- as.data.frame(df[which(df$preservation == 'random'),1:5])
head(colin)
colin$Preservation_no <- NULL
colin$imbalance <- df$index[which(df$preservation == 'random')]
head(colin)


vif(colin)
vifcor(colin, th = .9)
#the VIF test would suggest that we don't have a major (>90%) co-linearity probelm, our biggest correlation is 
#between starting lambda and mu. This only really tells us about our continous predictors though
#to examine correlations among our models we will construct a correlation heatmap

#vif(as.data.frcolin[,1:2])
cormat <- round(cor(colin), 2)

# Get lower triangle of the correlation matrix


upper_tri <- get_upper_tri(cormat)
upper_tri

melted_cormat <- melt(upper_tri, na.rm = TRUE)

colours <- c('#264653', "#2a9d8f", "#f4a261", "#e76f51")

# Heatmap
library(ggplot2)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = colours[1], high = colours[4], mid = colours[2], 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  coord_fixed()+ 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 5) +
  scale_x_discrete(labels = c("N extinct", "Age", expression(paste('Initial ', lambda)), expression(paste('Initial ',mu)), 
                              "CI"))+
  scale_y_discrete(labels = c("N extinct", "Age", expression(paste('Initial ',lambda)), expression(paste('Initial ',mu)), 
                              'CI'))+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text = element_text(face = 'bold.italic', size = 14),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.35,.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))

ggsave("figures/predictor_correlation_heatmap.pdf", device = 'pdf', scale = .7)
#okay, we can see that there is a very strong relationship between lambda and mu and between those two and Tree length.
#we can also see a correlation between Tree length and fossil number and a number of correlations between all of the above and 
#the tree buidling models








#therefore, believe the best fiting model 

best_model <- lmer(log(index) ~ model*preservation + mu_init + model + lambda_init + (1|tree_ID), data = df)
best_model <- lmer(log(index) ~ model*preservation*lambda_init*model*lambda_init*mu_init + (1|tree_ID), data = df)

AICc(best_model)
plot(best_model)
plot(resid(best_model))
qqnorm(resid(best_model, type = 'pearson'))
qqline(resid(best_model,  type = 'pearson'))
#not great

null <- lmer(log(index) ~ 1 + (1|tree_ID), data = df)
lrtest(best_model, null)
1-(logLik(best_model)/logLik(null))
summary(best_model)

#work out how much variance our random effect is capturing
0.05602/(0.05602 + 0.01594)
#okay, so this tells us that the Tree itself explains 78% of the variation in our dataset is explained by the Tree_ID

require(stargazer)




anova(best_model)

#emmeans struggles with large datasets so we are going to reduce the outpuit of the simulation by an order of magnitude. 
# models <- unique(df$model)
# preservation <- unique(df$preservation)
# new_df <- tibble()
# for(i in 1:length(models)){
#   for(j in 1:length(preservation)){
#     new_df <- rbind(new_df, df[which(df$model == models[i] & df$preservation == preservation[j]),][1:100,])
#   }
# }
# unique(new_df$preservation)
# 
# df <- new_df
# df


#lets subsample df

(head(df))


df$model <- as.character(df$model)
df$tree_ID <- as.character(df$tree_ID)

models <- unique(df$model)
models
ids <- vector()
for(i in 1:length(models)){  
  s <- df$tree_ID[which(df$model == models[i])]
  s <- unique(s)
  s <- s[1:(length(s)/2)]
  ids <- c(ids, s)
}

df_sub <- df[which(df$tree_ID %in% ids),]
df_sub$model <- factor(df_sub$model, levels = c("Constant", "Universal Speciation Shift", "Shifting Specation", "Shifting Extinction",
                                                "Shifting Speciation Extinction", "Mass Extinction Constant",
                                                "Mass Extinction Shifting", "Age Dependent", "Diversity Dependent", "Diversity Dependent Ki"))
df_sub$tree_ID <- as.factor(df_sub$tree_ID)


#we have to fit our model again
best_model <- lmer(log(index) ~ model*preservation + model*lambda_init + model*mu_init + lambda_init*mu_init + lambda_init*Tree_length + Fossil_no*model + (1|tree_ID), data = df_sub, na.action = 'na.fail')


require(lmerTest)

AICc(best_model)
summary(best_model)
anova(best_model)
# class(best_model) <- "lmerMod"
# stargazer(best_model, type = 'text', digits = 3, star.cutoffs = c(0.05,0.01,0.001), digit.separator = "")
# stargazer(best_model, digits = 3, star.cutoffs = c(0.05,0.01,0.001), digit.separator = "")

summary(best_model)
an_mod <- anova(best_model, test = 'ChiSq')
require(xtable)
an_mod

xtable(an_mod)

VarCorr(best_model)
confidence_intervals <- confint(best_model)
save(confidence_intervals, file = 'confidence_int.Rdata')
plot(best_model)

require(visreg)

pdf('figures/lambda_by_CI_fixed.pdf', width = 11, height = 11)
visreg(best_model, 'lambda_init', 'model', xlab = 'Initial lambda', ylab = 'log(Normalised Colless Index)', type = 'contrast', rug = 2, partial = F)
dev.off()

pdf('figures/mu_by_CI_fixed.pdf', width = 11, height = 11)
visreg(best_model, 'mu_init', 'model', xlab = 'Initial mu', ylab = 'log(Normalised Colless Index)', points=list(cex=.8, pch=1))
dev.off()

pdf('figures/extinct_by_CI_fixed.pdf', width = 11, height = 11)
visreg(best_model, 'Fossil_no', 'model', xtrans = log, xlab = 'log(N Extinct)', ylab = 'log(Normalised Colless Index)', type = 'contrast')
dev.off()

pdf('figures/labda_by_tree_age.pdf', width = 11, height = 11)
visreg(best_model, 'Tree_length', 'model', xtrans = log, xlab = 'log(Tree Age)', ylab = 'log(Normalised Colless Index)', layout = c(5,2))
dev.off()

visreg2d(best_model, 'lambda_init', 'mu_init', 'model')


best_model <- lmer(log(index) ~ model*preservation + model*lambda_init + model*mu_init + lambda_init*mu_init + lambda_init*Tree_length + Fossil_no*model + (1|tree_ID), data = df_sub, na.action = 'na.fail')

best_model1 <- lmer(log(index) ~ preservation + lambda_init + mu_init + lambda_init*mu_init + lambda_init*Tree_length + Fossil_no + (1|tree_ID), data = df_sub, na.action = 'na.fail')
best_model2 <- lmer(log(index) ~ model + model*lambda_init + model*mu_init + lambda_init*mu_init + lambda_init*Tree_length + Fossil_no*model + (1|tree_ID), data = df_sub, na.action = 'na.fail')
best_model3 <- lmer(log(index) ~ lambda_init + mu_init + lambda_init*mu_init + lambda_init*Tree_length + Fossil_no + (1|tree_ID), data = df_sub, na.action = 'na.fail')


anova(best_model, best_model1)
anova(best_model, best_model2)
anova(best_model, best_model3)
anova(best_model1, best_model3)
anova(best_model2, best_model3)




v1 + theme_classic()



visreg2d(best_model, 'lambda_init', 'mu_init', plot.type="rgl")



xtable(an_mod)
null_model <- lmer(log(index) ~ 1 + (1|tree_ID), data = df_sub, na.action = 'na.fail')

lrtest(best_model, null_model)

stargazer(anova(best_model),digit.separator = "", type = 'text') 

require(MuMIn)
dredger <- dredge(best_model, trace = T)
dredger
save(dredger, file = 'dredger_model_interactions.Rdata')

require(pbkrtest)
require(afex)
?mixed

#mixed(best_model, data = df)


setwd("../..")
source("functions.R")
setwd("BAMM_analyses/simulated_diversification/")

starting_random <- data_frame(draw_params(nsim = 1000, N = 100)[,1:2])
#starting_random$lambda_init <- starting_random$lambda_init[1]
#starting_random$mu_init <- starting_random$mu_init[1]
starting_random$Fossil_no <- 100
starting_random$Tree_length <- 10
starting_random$model <- as.factor(rep(unique(df$model), nrow(starting_random)/10))
starting_random$preservation <- as.factor('random')
starting_random$tree_ID <- as.factor(paste0('tree_', (seq(1, nrow(starting_random), 1))))
starting_random

starting_random$index <- predict(best_model, starting_random)

ggplot(starting_random, aes(x = mu_init/lambda_init, y = index))+
  geom_point()+
  facet_wrap(~ model, nrow = 5)+
  stat_smooth(method = 'lm') +
  theme_classic()





?predict




emm_options(pbkrtest.limit = 40001)
emm_options(lmerTest.limit = 40001)



system.time(emm3<- emmeans(best_model, specs = pairwise ~ model + preservation, type = 'response'))
emm3$contrasts
emm3$emmeans

save(emm3, file = "emm_means.Rdata")
load('emm_means.Rdata')


adjusted_means <- emm3$emmeans %>%as.data.frame()
adjusted_means$preservation <- factor(adjusted_means$preservation, levels = c("random", "diversified", "trait_dependent", "time_dependent"))
adjusted_means


adjusted_means$model <- factor(adjusted_means$model, levels = c("Constant", "Universal Speciation Shift", "Shifting Specation", "Shifting Extinction",
                                        "Shifting Speciation Extinction", "Mass Extinction Constant",
                                        "Mass Extinction Shifting", "Age Dependent", "Diversity Dependent", "Diversity Dependent Ki"))

ggplot(adjusted_means, aes(x = response, y = model, color = preservation))+
  scale_color_manual(labels = c("Random", "Diversified", "Trait Dependent", "Time Dependent"), values = c("#264653", "#2a9d8f", "#f4a261", "#e76f51")) +
  geom_pointrange(aes(xmin = upper.CL, xmax = lower.CL), position = position_dodge(width = .5)) +
  ylab(element_blank())+
  xlab("Adjusted normalised Colless Index") +
  labs(color = "Preservation model") +
  theme_classic()
ggsave('figures/emmeans_preservation_model.pdf', device = 'pdf', scale = .7)



emmip(best_model, model ~ Fossil_no | preservation,  cov.reduce = range)




#######################
### Defining Models ###
#######################

hist(log(df$index))

lmer.null <- lmer(log(index) ~ 1 + (1|tree_ID), data = df)
summary(lmer.null)

lmer.simple <- lmer(log(index) ~ model*preservation + (1|tree_ID), data = df)
summary(lmer.simple)

lmer.lammu <- lmer(log(index) ~ model*preservation + lambda_init + mu_init +(1|tree_ID), data = df)
summary(lmer.lammu)

lmer.div_turnover <- lmer(log(index) ~ model*preservation + intial_turnover + intial_diversification +(1|tree_ID), data = df)
summary(lmer.div_turnover)

anova(lmer.lammu, lmer.div_turnover)
anova(lmer.div_turnover, lmer.lammu)
AIC(lmer.lammu)
AIC(lmer.div_turnover)

#okay, so I think this all supports that the model with the initial lambda and mu is better 

lmer.lammu_model <- lmer(log(index) ~ model*preservation + model*lambda_init + model*mu_init + (1|tree_ID), data = df)
AIC(lmer.lammu_model)

lmer.lammu_model_interaction <- lmer(log(index) ~ model*preservation + model*lambda_init*mu_init + (1|tree_ID), data = df)
AIC(lmer.lammu_model_interaction)

lmer.lammu_model_fossil_no <- lmer(log(index) ~ model*preservation + lambda_init + mu_init + Fossil_no +(1|tree_ID), data = df)
AIC(lmer.lammu_model_fossil_no)

lmer.no_int <- lmer(log(index) ~ model + preservation + lambda_init + mu_init + Fossil_no +(1|tree_ID), data = df)
AIC(lmer.no_int)

lmer.lammu_model_fossil_no_tree_length <- lmer(log(index) ~ model*preservation + lambda_init + mu_init + Fossil_no + Tree_length + (1|tree_ID), data = df)
AIC(lmer.lammu_model_fossil_no_tree_length)

lmer.lammu_model_tree_length <- lmer(log(index) ~ model*preservation + lambda_init + mu_init + Tree_length + (1|tree_ID), data = df)
AIC(lmer.lammu_model_tree_length)

lmer.lammu_model_fossil_no_tree_length_complex <- lmer(log(index) ~ model*preservation + lambda_init + mu_init + Fossil_no*model + Tree_length*model + (1|tree_ID), data = df)
AIC(lmer.lammu_model_fossil_no_tree_length_complex)

lmer.lammu_model_fossil_no_tree_length_extra_complex <- lmer(log(index) ~ model*preservation + lambda_init*Tree_length + mu_init*Fossil_no + Fossil_no*model + Tree_length*model + (1|tree_ID), data = df)
AIC(lmer.lammu_model_fossil_no_tree_length_extra_complex)

lmer.lammu_model_fossil_no_tree_length_extra_complex <- lmer(log(index) ~ model*preservation + lambda_init*Tree_length*mu_init + mu_init*Fossil_no + Fossil_no*model + Tree_length*model + (1|tree_ID), data = df)


lmer.everything_interacts <- lmer(log(index) ~ model*preservation*lambda_init*Tree_length*mu_init*Fossil_no*Fossil_no*model*Tree_length*model + (1|tree_ID), data = df)
AIC(lmer.everything_interacts)

lmer.lammu_model_fossil_no_tree_length_extra_complex2 <- lmer(log(index) ~ model+preservation + lambda_init*Tree_length + mu_init*Fossil_no + Fossil_no*model + Tree_length*model + (1|tree_ID), data = df)
lmer.lammu_model_fossil_no_tree_length_extra_complex3 <- lmer(log(index) ~ model*preservation + lambda_init*Tree_length*mu_init*Fossil_no*model + (1|tree_ID), data = df)


anova(lmer.best, lmer.null, lmer.no_int, lmer.lammu_model_fossil_no, lmer.lammu_model_interaction, lmer.lammu_model, lmer.simple, lmer.lammu_model_fossil_no_tree_length, lmer.lammu_model_fossil_no_tree_length_complex, 
      lmer.lammu_model_fossil_no_tree_length_extra_complex, lmer.everything_interacts, lmer.lammu_model_fossil_no_tree_length_extra_complex2, lmer.lammu_model_fossil_no_tree_length_extra_complex3)


anova(lmer.best, lmer.lammu_model_fossil_no_tree_length_extra_complex, lmer.lammu_model_fossil_no_tree_length_extra_complex2)
anova(lmer.lammu_model_fossil_no_tree_length_extra_complex2, lmer.lammu_model_fossil_no_tree_length_extra_complex)

#okay so it seems that the best fitting model is the one that follows our intuition, one that includes all predictors but places an interaction terms on model and preservation,
#tree_length and fossil_no

#I think this bets model is based on my thinking 

lmer.best <- lmer(log(index) ~ model*preservation + lambda_init*mu_init*model + (1|tree_ID), data = df)
lmer.other <- lmer(log(index) ~ model*preservation + lambda_init*model + mu_init*model + (1|tree_ID), data = df)
lmer.worst <- lmer(log(index) ~ model*preservation + lambda_init + mu_init + (1|tree_ID), data = df)
lmer.terrible <- lmer(log(index) ~ model*preservation + lambda_init*model + (1|tree_ID), data = df)
lmer.terrible2 <- lmer(log(index) ~ model*preservation + mu_init*model + (1|tree_ID), data = df)

anova(lmer.best, lmer.other, lmer.worst, lmer.terrible, lmer.terrible2)

#and lmer.other comes out as the best, suggetsing we don't need an interaction between lambda and mu 
#but worth metioning that it is basically the equivelant of the first

ggplot(data = df, aes(x = Tree_length, y = model, color = preservation)) +
  geom_boxplot()
ggplot(data = df, aes(x = Fossil_no, y = model, color = preservation)) +
  geom_boxplot()

ggplot(data = df, aes(x = Tree_length, y = lambda_init, color = model, alpha = .6))+
  geom_point()

ggplot(data = df, aes(x = Fossil_no, y = lambda_init, color = model, alpha = .6))+
  geom_point()

#I think the inclusion of Fossil no and Tree Length is spurious, I think these are just picking up 
#on the differences between the models, rather than actually having an effect on imbalance

#so although the 

