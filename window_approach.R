### entropy 
library(tidyverse)

#### data in ####
setwd("/Users/andregomer/tempory_data/E1E2_revision/E1E2data_out")

entropy_df <- read_csv("entropy_df.csv")


setwd("/Users/andregomer/tempory_data/E1E2_revision/plots/")


##### colors ####
eqHV_colors <- c("horseA"='#081d58',"horseB"='#1d91c0', "horseC"='#7f0000',"horseD"='#b30000',"horseE"='#ef6548')
TV_colors <- c("horseF"="#3690c0", "horseG"="#74a9cf", "horseH"="#a6bddb","horseI"="#d0d1e6", 
               "inoculumII"="#212950", "inoculumI"="#56b3e6")
HCV_colors <- c("patientJ"="#6a51a3", "patientK"="#993404", "patientL"="#cc4c02")




#### function ####

slideFunct <- function(data, window, step){
  total <- length(data)
  spots <- seq(from=1, to=(total-window), by=step)
  result <- vector(length = length(spots))
  pos <- c()
  for(i in 1:length(spots)){
    result[i] <- mean(data[spots[i]:(spots[i]+window)])
    pos <- c(pos, i)
  }
  results <- data.frame(pos=pos, values=result)
  return(results)
}



# 
slide_output <- data.frame()
win <- 20
step <- 6


slide_df <- entropy_df %>% filter(annotation_track=="E1")
for (i in levels(factor(slide_df$sample_id))) {
  sample_loop <- slide_df %>% filter(sample_id==i)
  
  slide_data <- slideFunct(sample_loop$iSNV_freq, win, step)
  slide_data_consensus <- slideFunct(1- sample_loop$iSNV_freq, win, step)
  slide_data["sample"] <- rep(i, nrow(slide_data))
  slide_data["TVid"] <- rep(sample_loop$TVid[1], nrow(slide_data))
  slide_data["species"] <- rep(sample_loop$species[1], nrow(slide_data))
  slide_data["status"] <- rep(sample_loop$status[1], nrow(slide_data))
  slide_data["consensus_freq"] <- slide_data_consensus$values
  
  
  slide_output <- rbind(slide_output,slide_data)
}
slide_outputE1 <- slide_output
slide_outputE1$prot <- "E1"


# 
slide_output <- data.frame()
slide_df <- entropy_df %>% filter(annotation_track=="E2")
for (i in levels(factor(slide_df$sample_id))) {
  sample_loop <- slide_df %>% filter(sample_id==i)
  
  slide_data <- slideFunct(sample_loop$iSNV_freq, win, step)
  slide_data_consensus <- slideFunct(1- sample_loop$iSNV_freq, win, step)
  slide_data["sample"] <- rep(i, nrow(slide_data))
  slide_data["TVid"] <- rep(sample_loop$TVid[1], nrow(slide_data))
  slide_data["species"] <- rep(sample_loop$species[1], nrow(slide_data))
  slide_data["status"] <- rep(sample_loop$status[1], nrow(slide_data))
  slide_data["consensus_freq"] <- slide_data_consensus$values
  
  
  slide_output <- rbind(slide_output,slide_data)
}
slide_outputE2 <- slide_output
slide_outputE2$prot <- "E2"

slide_output <- rbind(slide_outputE1, slide_outputE2)



slide_output_stats <- slide_output %>%
  group_by(pos, TVid, species, status, prot) %>%
  summarise(mean_freq=mean(values), mean_consensus=mean(consensus_freq))






#### run function twice - E1 and E2 sep
## make metrics acute chronic HCV etc
slide_output_stats$class <- NA

slide_output_stats <-  slide_output_stats %>% mutate(class=replace(class, species=="Nat" & status == "chronic", "Natchronic"))
slide_output_stats <- slide_output_stats %>% mutate(class=replace(class, species=="Nat" & status == "acute", "Natacute"))
slide_output_stats <-  slide_output_stats  %>% mutate(class=replace(class, species=="HCV", "HCV"))
slide_output_stats <-  slide_output_stats  %>% mutate(class=replace(class, species=="TV", "TV"))
slide_output_stats <-  slide_output_stats  %>% mutate(class=replace(class, TVid=="TVD" | TVid=="TVP", "Inoculum"))


# set order
slide_output_stats$class <- factor(slide_output_stats$class, levels=c("TV", "Natacute","Natchronic", "HCV"))


#### plot ####


slide_output_stats %>%
  ggplot(aes(x=pos, y=mean_consensus, color=TVid))+
  geom_line()+
  theme_classic()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))+
  labs(x="Window", y="Mean consensus frequency")+
  scale_color_manual(values = c(HCV_colors, TV_colors, eqHV_colors), name="")+
  scale_x_continuous(breaks = seq(0, 550, 10))+
  facet_grid(class~prot)
ggsave("windowPlot.pdf", units="cm", height=8, width=21)


slide_output_stats %>%
  ggplot(aes(x=pos, y=mean_consensus, color=TVid))+
  geom_line()+
  theme_classic()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10),
        legend.position = "None")+
  labs(x="Window", y="Mean consensus frequency")+
  scale_color_manual(values = c(HCV_colors, TV_colors, eqHV_colors), name="")+
  scale_x_continuous(breaks = seq(0, 550, 10))+
  facet_grid(class~prot)
ggsave("windowPlot_fig.pdf", units="cm", height=5, width=21)





