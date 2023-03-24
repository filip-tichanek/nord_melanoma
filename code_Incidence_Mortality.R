rm(list=ls())
library(stringr)
library(tidyverse)
library(ggplot2)
library(cowplot)

# colors ---------------------------------------------------------

cola<-c(
  rgb(1,0.1,0.1),
  rgb(0.1,0.1,1),
  rgb(0, 0.6, 0.3),
  rgb(0.7,0.7,0.1))

colb<-c(
  rgb(1,0.1,0.1,alpha=0.2),
  rgb(0.1,0.1,1,alpha=0.2),
  rgb(0, 0.6, 0.3,alpha=0.2),
  rgb(0.7,0.7,0.1,alpha=0.2))

colc<-c(
  rgb(1,0.1,0.1,alpha=0.8),
  rgb(0.1,0.1,1,alpha=0.8),
  rgb(0, 0.6, 0.3,alpha=0.8),
  rgb(0.7,0.7,0.1,alpha=0.8))


# Data upload---------------------------------------------------------

urlfile <- "https://raw.githubusercontent.com/filip-tichanek/nord_melanoma/main/source_code/inc_mor_melanoma.csv"

melanoma_inc_mor <- read.csv(url(urlfile), sep=",")
colnam<-melanoma_inc_mor[,1]
melanoma_inc_mor<-data.frame(t(melanoma_inc_mor))[-1,]
colnames(melanoma_inc_mor) <- colnam
melanoma_inc_mor$year <- 1943:2020
### Sub-setting the years 1961-2020
melanoma_inc_mor<-melanoma_inc_mor[melanoma_inc_mor$year>1960,]
### Removing space character and converting characters to numbers
x=1;repeat{
  melanoma_inc_mor[,x]<-str_trim(melanoma_inc_mor[,x])
  melanoma_inc_mor[,x]<-as.numeric(melanoma_inc_mor[,x])
  x=x+1;if(x>17){break}}
summary(melanoma_inc_mor)


tidy_df <- melanoma_inc_mor %>% pivot_longer(-year, names_to = "outcome", values_to = "value")

# Plot the data

male_inc <- tidy_df %>% filter(
  grepl("incidence",outcome),
  grepl(", males",outcome)
)

male_inc <- male_inc %>% mutate(
  outcome = as.factor(outcome)
)


female_inc <- tidy_df %>% filter(
  grepl("incidence",outcome),
  grepl(", females",outcome)
)

female_inc <- female_inc %>% mutate(
  outcome = as.factor(outcome)
)



male_mor <- tidy_df %>% filter(
  grepl("mortality",outcome),
  grepl(", males",outcome)
)

male_mor <- male_mor %>% mutate(
  outcome = as.factor(outcome)
)


female_mor <- tidy_df %>% filter(
  grepl("mortality",outcome),
  grepl(", females",outcome)
)

female_mor <- female_mor %>% mutate(
  outcome = as.factor(outcome)
)




## new
p1 <- ggplot(male_inc, aes(x = year, y = value, color = outcome)) +
  geom_smooth(method = "loess",se=FALSE,lwd=0.8,span=0.4) +
  geom_smooth(data = subset(male_inc, outcome == "Finland, incidence, males"),
                method = "loess", level = 0.95, alpha = 0.5,lwd=0.8,span=0.4) +
  labs(x = "Year", y = "Incidence per 100,000 (ASR - World)") +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = c(0.3, 0.7),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.key.size = unit(1.2, "lines"),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, vjust = 0.2),
        axis.text.y = element_text(size = 12, hjust = 0.2),
        axis.title.x = element_text(size = 12, margin = margin(t = 1)),
        axis.title.y = element_text(size = 12, margin = margin(r = 1))
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 0.2)),
         fill = guide_legend(override.aes = list(alpha = 0.4)),
         shape = FALSE) +
  annotate("text",1990,30,label = "Males",size =5)+
  scale_color_manual(values = c(cola), labels = c( "Denmark","Finland","Norway", "Sweden")) +
  scale_fill_manual(values = c("Finland, incidence, males" = cola[2]))+
  scale_y_continuous(breaks = seq(0, 30, 2))+
  scale_x_continuous(breaks = seq(1970, 2020, 10))+
  ylim(0,30)+ 
  xlim(1961,2020)




p2 <- ggplot(female_inc, aes(x = year, y = value, color = outcome)) +
  geom_smooth(method = "loess",se=FALSE, lwd=0.8,span=0.4) +
  geom_smooth(data = subset(female_inc, outcome == "Finland, incidence, females"),
              method = "loess", level = 0.95, alpha = 0.5,lwd=0.8,span=0.4) +
  labs(x = "Year", y = NULL) +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = 'none',
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, vjust = 0.2),
        axis.text.y = element_text(size = 12, hjust = 0.2),
        axis.title.x = element_text(size = 12, margin = margin(t = 1)),
        axis.title.y = element_text(size = 12, margin = margin(r = 1))
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 0.2)),
         fill = guide_legend(override.aes = list(alpha = 0.4)),
         shape = FALSE) +
  annotate("text",1990,30,label = "Females",size =5)+
  scale_color_manual(values = c(cola), labels = c( "Denmark","Finland","Norway", "Sweden")) +
  scale_fill_manual(values = c("Finland, incidence, males" = cola[2]))+
  scale_y_continuous(breaks = seq(0, 30, 2))+
  scale_x_continuous(breaks = seq(1970, 2020, 10))+
  ylim(0,30)+ 
  xlim(1961,2020)



p3 <- ggplot(male_mor, aes(x = year, y = value, color = outcome)) +
  geom_smooth(method = "loess",se=FALSE,lwd=0.8,span=0.4) +
  geom_smooth(data = subset(male_mor, outcome == "Finland, mortality, males"),
              method = "loess", level = 0.95, alpha = 0.5,lwd=0.8,span=0.4) +
  labs(x = "Year", y = "Mortality per 100,000 (ASR - World)") +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = 'none',
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, vjust = 0.2),
        axis.text.y = element_text(size = 12, hjust = 0.2),
        axis.title.x = element_text(size = 12, margin = margin(t = 1)),
        axis.title.y = element_text(size = 12, margin = margin(r = 1))
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 0.2)),
         fill = guide_legend(override.aes = list(alpha = 0.4)),
         shape = FALSE) +
  scale_color_manual(values = c(cola), labels = c( "Denmark","Finland","Norway", "Sweden")) +
  scale_fill_manual(values = c("Finland, mortality, males" = cola[2]))+
  scale_y_continuous(breaks = seq(0, 30, 2))+
  scale_x_continuous(breaks = seq(1970, 2020, 10))+
  ylim(0,5)+ 
  xlim(1961,2020)


p4 <- ggplot(female_mor, aes(x = year, y = value, color = outcome)) +
  geom_smooth(method = "loess",se=FALSE,lwd=0.8,span=0.4) +
  geom_smooth(data = subset(female_mor, outcome == "Finland, mortality, females"),
              method = "loess", level = 0.95, alpha = 0.5,lwd=0.8,span=0.4) +
  labs(x = "Year", y = NULL) +
  theme(panel.background = element_rect(fill = "grey90"),
        legend.position = 'none',
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(size = 12, vjust = 0.2),
        axis.text.y = element_text(size = 12, hjust = 0.2),
        axis.title.x = element_text(size = 12, margin = margin(t = 1)),
        axis.title.y = element_text(size = 12, margin = margin(r = 1))
  ) +
  guides(color = guide_legend(override.aes = list(alpha = 0.2)),
         fill = guide_legend(override.aes = list(alpha = 0.2)),
         shape = FALSE) +
  scale_color_manual(values = c(cola), labels = c( "Denmark","Finland","Norway", "Sweden")) +
  scale_fill_manual(values = c("Finland, mortality, males" = cola[2]))+
  scale_y_continuous(NULL)+
  scale_x_continuous(breaks = seq(1970, 2020, 10))+
  ylim(0,5)+ 
  xlim(1961,2020)


plot_grid(p1, p2, p3,p4, labels=c("a", "b","c","d"),label_x = 0.16,
          label_y = 0.99, ncol = 2, axis='tb',
          label_size = 20,align = "v", nrow = 2, rel_heights = c(1, 1))

