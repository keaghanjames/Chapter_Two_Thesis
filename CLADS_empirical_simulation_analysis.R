setwd("/Users/keaghanyaxley/Dropbox/imbalance_combined_evidence/")
source("functions.R")

study <- 'Brennan_in_press'
taxa <- 'Macropodinae'

# study <- 'Gustafosn_et_al_2017'
# taxa <- 'Gyrinidae'

#study <- 'Herrera_and_davalos'
# taxa <- 'lemur'

# study <- 'Koch_et_al_2020'
# taxa <- 'Echinoidea'

#study <- 'Dornburgh_et_al_2015'
#taxa <- 'Holocentridae'

# study <- 'Paterson_et_al_2020'
# taxa <- 'Pinipedia'

# study <- 'Kealy_and_Beck_2017'
# taxa <- 'Dasyuromorphia'

# study <- 'Slater_et_al_2017'
# taxa <- 'Mysticeti'

setwd(paste0("BAMM_analyses/", study))


#df <- read.csv(paste0('simulation_summary_', study,'.csv')) 
df <- read.csv('ClaDS_summary.csv')
#df$X <- NULL               

df <- as_tibble(df)
summary(df)

summary(as.factor(df$model))

# df <- df[which(df$model != 'age_dependent'),]
# 
# #add synthetic data that is less imabalnced to fill the gaps
# 
# 
# new_df <- df[which(df$model == 'con_ME'),]
# new_df[1:5] <- 0.0000001
# new_df$model <- 'age_dependent'
# 
# new_df <- bind_rows(new_df, new_df)
# new_df$model[(nrow(new_df)/2):nrow(new_df)] <- 'con'
# summary(as.factor(new_df$model))
# 
# df <- bind_rows(new_df, df)
# summary(as.factor(df$model))

nrow(df)
df <- df %>%
  group_by(model) %>%
  sample_n(1000)
nrow(df)
summary(as.factor(df$model))


#tree <- read.tree(paste0(taxa,'.tree'))
tree <- read.tree('trim.tree')
plot(tree, show.tip.label = F)
colnames(df)
unique(df$model)
df$model <- as.factor(df$model)

df$tree_ID <- paste0('tree_ID_', seq(1,nrow(df)))


gathercols <- c("CI_random","CI_diversified", "CI_time_dependent", "CI_trait_dependent","CI_all",
                "SI_random","SI_diversifed","SI_time_dependent","SI_trait_dependent",   
                "SI_all", "Gamma_random", "Gamma_diversified", "Gamma_time_dependent",
                "Gamma_trait_dependent", "Gamma_all")


#df_long[,which(df_long$preservation == "CI_random" | df_long$preservation == "CI_diversified" | 
#        df_long$preservation == "CI_trait_dependent" | df_long$preservation == "CI_time_dependent")]

df_long <- gather(df, preservation, index, gathercols)
df_long$preservation <- as.factor(df_long$preservation)       
df_long
df_CI <- df_long[which(df_long$preservation == "CI_random" | df_long$preservation == "CI_diversified" | 
                         df_long$preservation == "CI_trait_dependent" | df_long$preservation == "CI_time_dependent"),]

df_CI$preservation <- str_remove(df_CI$preservation, "CI_")       
df_CI$preservation <- as.factor(df_CI$preservation)
#df_CI$preservation <- factor(df_CI$preservation, levels = c("random", "trait_dependent", "time_dependent", "diversified"))


levels(df_CI$model) <- c('AD', "C", "C+ME", 'DD', 'DD+KI', 'E', 'S','S+E', 'S+ME', 'US')

#levels(df_CI$model) <- c( "C+ME", 'DD', 'DD+KI', 'E', 'S','S+E', 'S+ME', 'US', 'AD', "C")

#df_CI$model <- ordered(df_CI$model, levels = c( "con", "universal_sp_shift", "shifting_sp",  "shifting_ex", "shifting_sp_ex", "con_ME" , "shiftsp_ME" ,
                      #         "age_dependent", "diversity_dependent", "diversity dependent w KI"))


#levels(df_CI$model) <- c("C", "US", "S", "E", "S+E", "C+ME",
#                         "S+ME", "AD", "DD", "DD+KI")

#df_CI$model <- as.factor(df_CI$model)
  
#levels(df_CI$preservation) <- c("Random", "Trait Dependent", "Time Dependent", "Diversified")

write.csv(df_CI, file = "CI_simulation_empirical_output.csv", row.names = F)
#df_CI <- read.csv('CI_simulation_empirical_output.csv')

colours <- c('#264653', "#2a9d8f", "#f4a261", "#e76f51")


df_CI <- df_CI %>% arrange(model)
df_CI$id <- seq(1, nrow(df_CI))

# # Get the name and the y position of each label
# label_data <- df_CI
# number_of_bar <- nrow(label_data)
# angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
# label_data$hjust <- ifelse( angle < -90, 1, 0)
# label_data$angle <- ifelse(angle < -90, angle+180, angle)
# 
# empty_bar <- 3
# 
# # prepare a data frame for base lines
# base_data <- df_CI %>% 
#   group_by(model) %>% 
#   summarize(start=min(id), end=max(id) - empty_bar) %>% 
#   rowwise() %>% 
#   mutate(title=mean(c(start, end)))
# 
# #prepare a data frame for grid (scales)
# grid_data <- base_data
# grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
# grid_data$start <- grid_data$start - 1
# grid_data <- grid_data[-1,]
# 
# 
# p1 <- ggplot(df_CI, aes(y = index, x = model, fill = preservation, col = preservation))+
#   geom_hline(yintercept = 0, color = 'lightgrey', lwd = .5, lty = 1) +
#     geom_hline(yintercept = .25, color = 'lightgrey', lwd = .5, lty = 1) +
#   geom_hline(yintercept = .5, color = 'lightgrey', lwd = .5, lty = 1) +
#   geom_hline(yintercept = .75, color = 'lightgrey', lwd = .5, lty = 1) +
#   geom_hline(yintercept = 1, color = 'lightgrey', lwd = .5, lty = 1) +
#   geom_label(x = .5, y = 1, label = '1', col = 'darkgrey', fill = 'white') +
#   geom_label(x = .5, y = .75, label = '0.75', col = 'darkgrey', fill = 'white') +
#   geom_label(x = .5, y = .5, label = '0.5', col = 'darkgrey', fill = 'white') +
#   geom_label(x = .5, y = .25, label = '0.25', col = 'darkgrey', fill = 'white') +
#   geom_label(x = .5, y = 0, label = '0', col = 'darkgrey', fill = 'white') +
#   
#   
#   geom_violin(size = 2, position = position_dodge(0.5)) +
#   
#  
# 
#   
#   geom_hline(yintercept = ci.norm(tree), color = 'black', lwd = .7, lty = 2) +
#   scale_fill_manual(values = colours) +
#   scale_color_manual(values = colours)+
#   theme_classic() +
#   ylab("") +
#   xlab("Normalised Colless Index") +
#   ylim(-0.5,1.1) +
#   theme(
#     legend.position = "bottom",
#     axis.text = element_blank(),
#     axis.title = element_blank(),
#     panel.grid = element_blank(),
#     plot.margin = unit(rep(-1,4), "cm") 
#   ) +
#   coord_polar(start = 0)
# p1
#ggsave(paste0(study,"_violin_exceed_threshold.pdf"), device = 'pdf', scale = .7)

summary(df_CI$model)



p2 <- ggplot(df_CI, aes(x = lambda_init-mu_init, y = index, col = log(Tree_length)))+
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~ model)
p2


counts <- df_CI %>%
  group_by(model, preservation) %>%
      summarise(length(which(index >= ci.norm(tree))))
colnames(counts) <- c('model', 'preservation', 'count')
counts$percentage <- counts$count/1000
#counts <- counts[1:14,]
counts$count
counts
tail(counts)
# per_past_threshold <- function(data, threshold, log = T){
#   #  models <- as.factor(data[,models])
#   levels <- levels(as.factor(data$model))
#   pres <- levels(data$preservation)
#   excess <- data.frame()
#   for(i in 1:length(levels)){
#     a <- which(data$model == levels[i])
#     s <- vector()
#     for(j in 1:length(pres)){
#       s <- c(s, (length(which(data$index[a] >= threshold & data$preservation[a] == pres[j]))))#/length(data$preservation[a] == pres[j]))
#     }
#     d <- data_frame(count = s, percentage = s/length(data$preservation[a] == pres[1]), model = levels[i], preservation = pres)
#     if(log == T){
#       d$percentage <- log((d$percentage)*100)
#     }
#     excess <- rbind(excess, d)
#   }
#   return(excess)
# }
# counts <- per_past_threshold(df_CI, ci.norm(tree), log = F)

#counts$model <- ordered(counts$model, levels = levels(df_CI$model))
counts$preservation <- ordered(counts$preservation, levels = levels(df_CI$preservation))
#counts$percentage[which(is.infinite(counts$percentage) == T)] <- 0
counts

save(counts, file = 'Clad_threshold_counst_CI.Rdata')
counts$count
empty_bar <- 3
to_add <- data.frame( matrix(NA, empty_bar*nlevels(counts$model), ncol(counts)) )
colnames(to_add) <- colnames(counts)
to_add$model <- rep(levels(counts$model), each=empty_bar)
counts <- rbind(counts, to_add)

counts <- counts %>% arrange(model)
extra <- counts[5:7,]
counts <- rbind(extra, counts)

counts$id <- seq(1, nrow(counts))

# Get the name and the y position of each label
label_data <- counts
number_of_bar <- nrow(label_data)
angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust <- ifelse( angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# prepare a data frame for base lines
base_data <- counts %>% 
  group_by(model) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

base_data[1,2] <- 4

#prepare a data frame for grid (scales)
grid_data <- base_data
grid_data$end <- grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start <- grid_data$start - 1
grid_data <- grid_data[-1,]

# Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.

counts$count <- counts$count/10


p3 <- ggplot(counts, aes(x=as.factor(id), y=log(count+1), fill=preservation)) +
  geom_bar(aes(x=as.factor(id), y=log(count+1), fill=preservation), stat="identity", alpha=1) +
  
  # geom_segment(aes(x = 2, y = log(2), xend = max(counts$id)-1, yend = log(2)), color ='grey', linetype = 'dashed')+
  geom_segment(aes(x = 2, y = log(6), xend = max(counts$id)-1, yend = log(6)), color ='lightblue')+
  geom_segment(aes(x = 2, y = log(26), xend = max(counts$id)-1, yend = log(26)), color ='grey', linetype = 'dashed')+
  geom_segment(aes(x = 2, y = log(51), xend = max(counts$id)-1, yend = log(51)), color ='grey', linetype = 'dashed')+
  geom_segment(aes(x = 2, y = log(76), xend = max(counts$id)-1, yend = log(76)), color ='grey', linetype = 'dashed')+
  geom_hline(yintercept = log(100), color = 'grey')  +
  
  
  geom_bar(aes(x=as.factor(id), y=log(count+1), fill=preservation), stat="identity", alpha=1) +
  
  
  #geom_segment(data=grid_data, aes(x = end, y = log(5), xend = start, yend = log(5)), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = log(15), xend = start, yend = log(15)), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  #geom_segment(data=grid_data, aes(x = end, y = log(25), xend = start, yend = log(25)), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
  # Add text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(counts$id),4), y = c(log(6), log(26), log(51), log(76)), label = c("ln(50)", "ln(250)", "ln(500)", "ln(750)") , color="grey", size=3 , angle=0, fontface="bold", hjust=0) +
  
  
  scale_fill_manual(values = colours) +
  ylim(-3,log(101)) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  ylab("") +
  coord_polar(start = 0) + 
  # geom_text(data=label_data, aes(x=id, y=percentage, label=preservation, hjust=hjust), 
  #           color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) +
  # 
  # # Add base line information
 geom_segment(data=base_data, aes(x = start, y = -0.01, xend = end, yend = -0.01), 
              colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  geom_text(data=base_data, aes(x = title, y = -.6, label=model),  
             colour = "black", alpha=0.8, size=4, hjust = c(-.5,.3,.5,.5,0.5,0.5,0.3,0.4,0.4,0.5), fontface="bold", inherit.aes = FALSE) 
  # 

  #xlab("% of simulated trees with imbalance >= empirical tree")
p3

counts$percentage
ggsave(paste0(study,"_ClaDS_exceed_threshold.pdf"), device = 'pdf', scale = .7)

require(ggridges)
qant <- function(x) quantile(x, probs = c(.05, 0.95))

p4 <- ggplot(df_CI, aes(x = index, y = model, col = preservation))+
  geom_violin() +
  scale_color_manual(values = colours) +
  scale_shape(solid = F) +
  stat_summary(fun = 'qant',
               geom = 'crossbar',
               width = 0.6,
               color = 'darkgrey') +
  stat_summary(fun = 'median',
               geom = 'crossbar',
               width = 0.8,
               color = 'black') +
  geom_vline(xintercept = ci.norm(tree), col = 'lightblue')+
  theme_classic()
p4

median_95 <- df_CI %>%
  group_by(model, preservation) %>%
    summarise(five = quantile(x = index, probs = 0.05),
              fifty = quantile(x = index, probs = 0.5),
              ninetyfive = quantile(x = index, probs = 0.95)
    )
median_95                               


p5 <- ggplot(median_95, aes(x = fifty, y = model, color = preservation))+
               geom_pointrange(aes(xmin = five, xmax = ninetyfive), position = position_dodge2(width = .4)) +
  scale_color_manual(values = colours) +
  geom_vline(xintercept = ci.norm(tree), col = 'lightblue') +
  theme_classic()
  
p5

