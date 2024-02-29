### Mutation-tracing plot ####


###### library #####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)



##### in data ####
dnds_df <- read_csv("/Users/andregomer/tempory_data/E1E2_revision/data/20211201dnds_df.csv")
entropy_df <- read_csv("/Users/andregomer/tempory_data/E1E2_revision/data/20211203_entropy_data.csv")
AA_table <- read_csv('/Users/andregomer/tempory_data/E1E2_revision/E1E2data_out/AAcodonTable.csv')


##### set plotting variables ####
plotdir <- "/Users/andregomer/tempory_data/E1E2_revision/plots2"
setwd(plotdir)


##### ORF E1E2 #####
HVRs_eq <- 188
HVRe_eq <- 215
HVRs_HCV <- 193
HVRe_HCV <- 219

##### colors #####
cbPalette <-  c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#56B4E9")
eqHV_colors <- c('#081d58','#1d91c0', '#7f0000','#b30000','#ef6548')
TV_colors <- c("#3690c0", "#74a9cf", "#a6bddb","#d0d1e6", "#081d58", "#081d58")
HCV_colors <- c("#6a51a3", "#993404", "#cc4c02")


##### Set Factor levels in entropy data ######
entropy_df$TVid <- factor(entropy_df$TVid, levels = c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                        "horseF", "horseG", "horseH", "horseI",
                                                        "patientJ", "patientK", "patientL",
                                                        "InoculumI", "InoculumII"))

AA_table$TVid <-  factor(AA_table$TVid, levels = c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                      "horseF", "horseG", "horseH", "horseI",
                                                      "patientJ", "patientK", "patientL",
                                                      "inoculumI", "inoculumII"))

##### Set Factor levels in dNdS data ######
dnds_df$TVid <- factor(dnds_df$TVid, levels = c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                  "horseF", "horseG", "horseH", "horseI",
                                                  "patientJ", "patientK", "patientL",
                                                  "InoculumI", "InoculumII"))

entropy_df$sample_date <- factor(entropy_df$sample_date, levels = c("week01", "week02", "week03", "week04", "week05", "week06", 
                                                                      "week07", "week08", "week09", "week10", "week11", "week12", "week29",
                                                                      "week45", "week105", "week116", "week186", "week225", "week242", "week295",
                                                                      "week319"))

dnds_df$sample_date <- factor(dnds_df$sample_date, levels = c("week01", "week02", "week03", "week04", "week05", "week06", 
                                                                "week07", "week08", "week09", "week10", "week11", "week12", "week29",
                                                                "week45", "week105", "week116", "week186", "week225", "week242", "week295",
                                                                "week319"))

AA_table$sample_date <- factor(AA_table$sample_date, levels = c("week01", "week02", "week03", "week04", "week05", "week06", 
                                                                "week07", "week08", "week09", "week10", "week11", "week12", "week29",
                                                                "week45", "week105", "week116", "week186", "week225", "week242", "week295",
                                                                "week319"))




#### variant tracing AA-level EqHV ####
AA_table_EqHV_sel <- AA_table %>% filter(species == "Nat" & variant_freq >= .1)

AA_table_EqHV_pos <- unique(AA_table_EqHV_sel$Pos)
AA_table_EqHV_pos

AA_table_EqHV <- AA_table %>% filter(species == "Nat" & Pos %in% AA_table_EqHV_pos)


## bin low, high, consensus labels
iSNV_thres <- .03287
iSNV_low <- iSNV_thres
iSNV_medium <- .1
iSNV_high <- .5
AA_table_EqHV["mut_category"] <- "Low frequency"

# add mutation label 
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq >= iSNV_high, ("Consensus")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_high & variant_freq >= iSNV_medium, ("High frequency")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_medium & variant_freq >= iSNV_low, ("Medium frequency")))



#



AA_table_EqHV_spacer <- AA_table_EqHV[1:4,]
AA_table_EqHV_spacer$sample_id <- c("SPACER1", "SPACER2","SPACER3","SPACER4")
AA_table_EqHV_spacer$mut_category <- NA

AA_table_EqHV2 <- rbind(AA_table_EqHV_spacer, AA_table_EqHV)

AA_table_EqHV2$sample_id <- factor(AA_table_EqHV2$sample_id, levels=c(
  "horseA week01", "horseA week04", "SPACER1",
  "horseB week01", "horseB week07", "horseB week08", "SPACER2",
  "horseC week01", "horseC week04", "horseC week28", "horseC week54","horseC week82","SPACER3",
  "horseD week01", "horseD week08", "horseD week28", "horseD week54", "horseD week82", "SPACER4",
  "horseE week01", "horseE week08", "horseE week19", "horseE week37", "horseE week94"
  ))

EqHV_variants <- AA_table_EqHV2 %>%
  ggplot(aes(x=as.factor(Pos), y=(sample_id), fill=mut_category))+
  geom_point(size=6, alpha=.7, shape=21, color="black")+
  #geom_point(aes(shape=dnds_category,size=dnds_category), fill="transparent", color="black")+
  scale_shape_manual(values=c(95, 20, 3))+
  scale_size_manual(values=c(5,2,3))+
  scale_fill_manual(values=c("#67000d","#fdae61","transparent","#4393c3"))+
  geom_hline(yintercept = c(3, 7, 13, 19))+
  geom_vline(xintercept = c(16.5))+
  labs(x="Position", y="")+
  theme_classic()+
  annotate("rect", xmin=0, xmax=16.5, ymin=25.5, ymax=26.5, color="black", fill="#aaacb2")+
  annotate("text", x=8.25,y=26, color="white", label="E1", size=5)+
  annotate("rect", xmin=16.5, xmax=42, ymin=25.5, ymax=26.5, color="black", fill="#7e8084")+
  annotate("text", x=29.75,y=26, color="white", label="E2", size=5)

EqHV_variants
#facet_wrap(~TVid, ncol=1, scales="free_y", strip.position = "left", shrink=FALSE)
ggsave("Variantplot_EqHV.pdf", height = 6.5, width = 13.5) 
#








#### variant tracing AA-level HCV ####
AA_table_EqHV_sel <- AA_table %>% filter(species == "HCV" & variant_freq >= .1)

AA_table_EqHV_pos <- unique(AA_table_EqHV_sel$Pos)
AA_table_EqHV_pos

AA_table_EqHV <- AA_table %>% filter(species == "HCV" & Pos %in% AA_table_EqHV_pos)


## bin low, high, consensus labels
iSNV_low <- iSNV_thres
iSNV_medium <- .1
iSNV_high <- .5
AA_table_EqHV["mut_category"] <- "Low frequency"

# add mutation label 
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq >= iSNV_high, ("Consensus")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_high & variant_freq >= iSNV_medium, ("High frequency")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_medium & variant_freq >= iSNV_low, ("Medium frequency")))




AA_table_EqHV_spacer <- AA_table_EqHV[1:2,]
AA_table_EqHV_spacer$sample_id <- c("SPACER1", "SPACER2")
AA_table_EqHV_spacer$mut_category <- NA

AA_table_EqHV2 <- rbind(AA_table_EqHV_spacer, AA_table_EqHV)

AA_table_EqHV2$sample_id <- factor(AA_table_EqHV2$sample_id, levels=c("patientJ week01",  "patientJ week29",  "patientJ week45", "SPACER2",
                                                                      "patientK week01",  "patientK week105", "patientK week242","SPACER1",
                                                                      "patientL week01","patientL week04",  "patientL week116", "patientL week186",
                                                                      "patientL week225", "patientL week295", "patientL week319"))








barheight <- 16
barxposE1 <- 15.5
barxposE2 <- barxposE1+47.5

HCV_variants <- AA_table_EqHV2 %>%
  ggplot(aes(x=as.factor(Pos), y=(sample_id), fill=mut_category))+
  geom_point(size=6, alpha=.7, shape=21, color="black")+
  #geom_point(aes(shape=dnds_category,size=dnds_category), fill="transparent", color="black")+
  scale_shape_manual(values=c(95, 20, 3))+
  scale_size_manual(values=c(5,2,3))+
  scale_fill_manual(values=c("#67000d","#fdae61", "transparent", "#4393c3"))+
  geom_hline(yintercept = c(4, 8))+
  geom_vline(xintercept = c(15.5))+
  labs(x="Position", y="")+
  theme_classic()+
  annotate("rect", xmin=0, xmax=barxposE1, ymin=barheight, ymax=barheight+1, color="black", fill="#aaacb2")+
  annotate("text", x=barxposE1/2,y=barheight+.5, color="white", label="E1", size=5)+
  annotate("rect", xmin=barxposE1, xmax=barxposE2, ymin=barheight, ymax=barheight+1, color="black", fill="#7e8084")+
  annotate("text", x=(barxposE1+barxposE2)/2,y=barheight+.5, color="white", label="E2", size=5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#facet_wrap(~TVid, ncol=1, scales="free_y", strip.position = "left", shrink=FALSE)
HCV_variants
ggsave("Variantplot_HCV.pdf", height = 4.5, width = 17.5)
#







#### variant tracing AA-level TV ####
AA_table_EqHV_sel <- AA_table %>% filter(species == "TV" & variant_freq >= .1)

AA_table_EqHV_pos <- unique(AA_table_EqHV_sel$Pos)
AA_table_EqHV_pos

AA_table_EqHV <- AA_table %>% filter(species == "TV" & Pos %in% AA_table_EqHV_pos)



## bin low, high, consensus labels
iSNV_low <- iSNV_thres
iSNV_medium <- .1
iSNV_high <- .5
AA_table_EqHV["mut_category"] <- "Low frequency"

# add mutation label 
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq >= iSNV_high, ("Consensus")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_high & variant_freq >= iSNV_medium, ("High frequency")))
AA_table_EqHV <- AA_table_EqHV %>% mutate(mut_category=replace(mut_category, variant_freq <= iSNV_medium & variant_freq >= iSNV_low, ("Medium frequency")))




AA_table_EqHV_spacer <- AA_table_EqHV[1:4,]
AA_table_EqHV_spacer$sample_id <- c("SPACER1", "SPACER2", "SPACER3", "SPACER4")
AA_table_EqHV_spacer$mut_category <- NA

AA_table_EqHV2 <- rbind(AA_table_EqHV_spacer, AA_table_EqHV)

AA_table_EqHV2$sample_id <- factor(AA_table_EqHV2$sample_id, levels=c(
  "horseF week02", "horseF week04", "horseF week06", "horseF week08", "horseF week10", "horseF week11", "horseF week12", "SPACER1",
  "horseG week01", "horseG week03", "horseG week05", "horseG week07", "horseG week08", "SPACER2",
  "horseH week01", "horseH week03", "horseH week05", "horseH week06", "horseH week07", "SPACER3",
  "horseI week01","horseI week02", "SPACER4",
  "inoculumI I", "inoculumII II"))


### add horseF week01.....


barheight <- 26
barxposE1 <- 3.5
barxposE2 <- barxposE1+5



TV_variants <- AA_table_EqHV2 %>%
  ggplot(aes(x=as.factor(Pos), y=(sample_id), fill=mut_category))+
  geom_point(size=6, alpha=.7, shape=21, color="black")+
  #geom_point(aes(shape=dnds_category,size=dnds_category), fill="transparent", color="black")+
  scale_shape_manual(values=c(95, 20, 3))+
  scale_size_manual(values=c(5,2,3))+
  scale_fill_manual(values=c("#67000d","#fdae61", "transparent", "#4393c3"))+
  geom_hline(yintercept = c(8, 14, 20, 23))+
  geom_vline(xintercept = c(3.5))+
  labs(x="Position", y="")+
  theme_classic()+
  annotate("rect", xmin=0, xmax=barxposE1, ymin=barheight, ymax=barheight+1, color="black", fill="#aaacb2")+
  annotate("text", x=barxposE1/2,y=barheight+.5, color="white", label="E1", size=5)+
  annotate("rect", xmin=barxposE1, xmax=barxposE2, ymin=barheight, ymax=barheight+1, color="black", fill="#7e8084")+
  annotate("text", x=(barxposE1+barxposE2)/2,y=barheight+.5, color="white", label="E2", size=5)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#facet_wrap(~TVid, ncol=1, scales="free_y", strip.position = "left", shrink=FALSE)
TV_variants
ggsave("Variantplot_TV.pdf", height = 6.75, width = 4.5) 
#



















