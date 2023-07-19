# Ok so lets get to data wrangling the luciferace data and setting up some metadata in a long format to start

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(ggpubr)
library(ggbreak)
library(scales)
library(cowplot)
library(multcompView)

luc_0710_df <- read.csv(file = "Luciferase_4h_WG (Modified)_20230710_162848.csv") # reads in data

colnames(luc_0710_df)[4:99] # well information columns 

luc_0710_long_df <- pivot_longer(luc_0710_df, cols = colnames(luc_0710_df)[4:99], # makes csv data into long format (easier for me to work with in ggplot)
                                     names_to = c("well"), values_to = c("luminescence"))

luc_0710_long_annotated_df <- luc_0710_long_df %>%  # annotates well information 
  mutate(promoter = case_when(str_detect(well, "3|10") == T ~ "G19",
                              str_detect(well, "4") == T ~ "EV",
                              str_detect(well, "5") == T ~ "35_S",
                              str_detect(well, "6|11") == T ~ "253",
                              str_detect(well, "7|12") == T ~ "255",
                              str_detect(well, "(?<=[:alpha:])2|9") == T ~ "256",
                              str_detect(well, "(1(?!.)|8)") == T ~ "249",)) %>% 
  mutate(biological_replicate = case_when(str_detect(well, "8|9|10|11|12|E4|F4|G4|H4|E5|F5|G5|H5") == T ~ "Replicate 1",
                                           str_detect(well, "(?<=[:alpha:])1|(?<=[:alpha:])2|3|4|5|6|7") == T ~ "Replicate 2")) %>% 
  filter(is.na(luminescence) == F)

# there are some wells that are pretty dead and likely had some mechanical error or did not get infiltrated so I am filtering them out 
# specifically wells G1, D7, B5 and maybe well D2
luc_0710_long_annotated_filtered_df <- luc_0710_long_annotated_df %>% 
  filter(well %in% c("H4","A4","A1","A8","H5") == F) %>% 
  group_by(well) %>% 
  reframe(luminescence = log(sum(luminescence)), promoter,biological_replicate) %>% 
  distinct(well, luminescence, promoter, biological_replicate)

# Lets add some factors
luc_0710_long_annotated_filtered_df$promoter <- as.factor(luc_0710_long_annotated_filtered_df$promoter)

# now order the factors so col is first
luc_0710_long_annotated_filtered_df$promoter = ordered(x=luc_0710_long_annotated_filtered_df$promoter, 
                                      c('35_S','G19','249','256','255','253','EV'))

# use to check levels of your factors
levels(luc_0710_long_annotated_filtered_df$promoter)


######################
# lets do a quick two way anova
luc_0710_long_annotated_filtered_summarized_df<-luc_0710_long_annotated_filtered_df %>% 
  group_by(well) %>% 
  filter(!promoter %in% c("35_S","G19"))

# https://statdoe.com/one-way-anova-and-box-plot-in-r/
anova <- aov(luminescence ~ promoter*biological_replicate, luc_0710_long_annotated_filtered_summarized_df)
summary(anova)
turkey <- TukeyHSD(anova)
cld<- multcompLetters4(anova, turkey)

table_of_anova <- luc_0710_long_annotated_filtered_summarized_df %>% 
  dplyr::group_by(promoter,biological_replicate) %>% 
  summarize(mean = mean(luminescence, na.rm = T), quant = quantile(luminescence, probs = 0.75)) %>% 
  arrange(desc(mean))

cld <- as.data.frame.list(cld$`promoter:biological_replicate`)
table_of_anova$cld <- cld$Letters

ggplot(data = luc_0710_long_annotated_filtered_summarized_df,
       mapping = aes(x=promoter, y = luminescence, fill = biological_replicate))+
  geom_boxplot()+
  geom_text(data = table_of_anova,
            aes(x = promoter, y = quant, group = biological_replicate,
                label = cld ),
            position = position_dodge(width = 1),
            size = 3, vjust=-1)+
  scale_y_continuous(labels = label_comma())+
  theme_linedraw()
  
###### lets do two one-way anovas ############## 

luc_plant_1 <- luc_0710_long_annotated_filtered_df %>% 
  filter(biological_replicate == "Replicate 1") %>% 
  filter(!promoter %in% c("G19","35_S","EV")) # removes non promoter groups


anova <- aov(luminescence ~ promoter, luc_plant_1)
summary(anova)
turkey <- TukeyHSD(anova)
cld<- multcompLetters4(anova, turkey)

table_of_anova <- luc_plant_1 %>% 
  dplyr::group_by(promoter) %>% 
  summarize(mean = mean(luminescence, na.rm = T), quant = quantile(luminescence, probs = 0.75)) %>% 
  arrange(desc(mean))

cld <- as.data.frame.list(cld$promoter)
table_of_anova$cld <- cld$Letters

plant_1_plot<-ggplot(data = luc_plant_1,
       mapping = aes(x=promoter, y = luminescence, fill = promoter))+
  geom_boxplot()+
  geom_jitter(alpha = 0.2)+
  geom_point()+
  geom_text(aes(label = well), hjust = 2)+
  geom_text(data = table_of_anova, aes(x = promoter, y = quant, label = cld ), size = 3, hjust=-1,vjust=-1)+
  scale_y_continuous(labels = label_comma())+
  labs(title = "Plant 1",
       caption = "
       
       
       ")+
  theme_linedraw()+
  theme(legend.position = "none")

## 

luc_plant_2 <- luc_0710_long_annotated_filtered_df %>% 
  filter(biological_replicate == "Replicate 2") %>% 
  filter(!promoter %in% c("G19","35_S","EV")) # removes non promoter groups


anova <- aov(luminescence ~ promoter, luc_plant_2)
summary(anova)
turkey <- TukeyHSD(anova)
cld<- multcompLetters4(anova, turkey)

table_of_anova <- luc_plant_2 %>% 
  dplyr::group_by(promoter) %>% 
  summarize(mean = mean(luminescence, na.rm = T), quant = quantile(luminescence, probs = 0.75)) %>% 
  arrange(desc(mean))

cld <- as.data.frame.list(cld$promoter)
table_of_anova$cld <- cld$Letters

plant_2_plot<-ggplot(data = luc_plant_2,
       mapping = aes(x=promoter, y = luminescence, fill = promoter))+
  geom_boxplot()+
  geom_point()+
  geom_text(aes(label = well), hjust = 2)+
  geom_jitter(alpha = 0.2)+
  geom_text(data = table_of_anova, aes(x = promoter, y = quant, label = cld ), size = 3, 
            vjust=-1, hjust=-1)+
  scale_y_continuous(labels = label_comma())+
  labs(title = "Plant 2",
       caption = "G19 confirmed to be contaminated by 35-S promoter in full plasmid sequencing and was ommited from analysis
       Data filtered for log(sum(luminescence)) for each well
       Anova+TukeyHSD test done for Luminescence ~ Promoter to generate letter labels
       Plants from week 4, week 2 transplants, Ag infiltrated Friday, experiment Monday
       Date: 07/10/23")+
  theme_linedraw()+
  theme(legend.position = "none")

grid_plot <-plot_grid(plant_1_plot,plant_2_plot, align = "h")

ggsave(filename = "promoter_Luc_experiment_log_sum_well_071023.png", grid_plot)





