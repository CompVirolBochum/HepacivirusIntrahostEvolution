###### Natural cohort variant plots #####
#


###### library #####
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)
library(gridExtra)
library(cowplot)



##### in data ####
setwd("/Users/andregomer/tempory_data/E1E2_revision/E1E2data_out")
AA_data <- read_csv("AA_df.csv")

dnds_df <- read_csv("dnds_df.csv")#

entropy_df <- read_csv("entropy_df.csv")

##### set plotting variables ####
plotdir <- "/Users/andregomer/tempory_data/E1E2_revision/plots/Figure2"
setwd(plotdir)



#### set thresholds & colors ####
iSNV_thres <- .03287
lSNV <- .01
hSNV <- .5
mSNV <- .1

color_dN <- "darkred"
color_dS <- "darkblue"


N_color <- "#01665e"
S_color <- "#8c510a"
total_color <-   "black"



##### colors #####
cbPalette <-  c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00", "#CC79A7", "#56b3e6",  "#212950")
eqHV_colors <- c('#081d58','#1d91c0', '#7f0000','#b30000','#ef6548')
TV_colors <- c("horseF"="#3690c0", "horseG"="#74a9cf", "horseH"="#a6bddb","horseI"="#d0d1e6", "inoculumII"="#212950", "inoculumI"="#56b3e6")
HCV_colors <- c("#6a51a3", "#993404", "#cc4c02")



##### Set Factor levels in entropy data ######
entropy_df$TVid <- factor(entropy_df$TVid, levels = c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                      "horseF", "horseG", "horseH", "horseI",
                                                      "patientJ", "patientK", "patientL", 
                                                      "inoculumI", "inoculumII"))


##### Set Factor levels in dNdS data ######
dnds_df$TVid <- factor(dnds_df$TVid, levels = c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                "horseF", "horseG", "horseH", "horseI",
                                                "patientJ", "patientK", "patientL",
                                                "inoculumI", "inoculumII"))

entropy_df$sample_date <- factor(entropy_df$sample_date, levels = c("week01", "week02", "week03", "week04", "week05", "week06",
                                                                    "week07", "week08", "week09", "week10", "week11", "week12", "week19", "week27","week28","week29", "week37",
                                                                    "week45", "week54", "week82", "week94","week105", "week116", "week186", "week225", "week242", "week295",
                                                                    "week319", "I", "II"))

dnds_df$sample_date <- factor(dnds_df$sample_date, levels =c("week01", "week02", "week03", "week04", "week05", "week06",
                                                             "week07", "week08", "week09", "week10", "week11", "week12", "week19", "week27","week28","week29", "week37",
                                                             "week45", "week54", "week82", "week94","week105", "week116", "week186", "week225", "week242", "week295",
                                                             "week319", "I", "II"))






##### Figure: Variant plot EqHV ####


# filter for naturally infection data only and set levels as facet labels
entropy_df_levels_eq <- entropy_df %>%
  filter(species=="Nat")
entropy_df_levels_eq <- levels(factor(entropy_df$TVid))
facet_labs_entropy_eq <- c("Horse A", "Horse B", "Horse C", "Horse D", "Horse E",
                           "Horse F", "Horse G", "Horse H", "Horse I",
                           "Patient J", "Patient K", "Patient L", "Inoculum I", "Inoculum II")
names(facet_labs_entropy_eq) <- entropy_df_levels_eq

### plots variant cloud
iSNV_eq <- entropy_df %>%
  filter(species=="Nat") %>%
  ggplot(aes(x=Position, y=iSNV_freq, color=as.factor(time_point)))+
  geom_point(size=.75)+
  geom_hline(yintercept = c(iSNV_thres), linetype="dotted")+
  geom_hline(yintercept = c(.5), linetype="dashed")+

  theme_classic()+
  xlab("Position")+
  ylab("iSNV frequency")+
  scale_color_manual(values = cbPalette, name="Timepoint")+
  scale_y_continuous(trans = "log2", breaks = c( .01, .1, .5, 1), limits = c(0.005, 1))+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)) )+
  facet_wrap(~TVid, nrow = 1, labeller = labeller(TVid = facet_labs_entropy_eq))
iSNV_eq

# isolate legend and remove from plot
iSNV_eq_legend <- get_legend(iSNV_eq)
iSNV_eq <- iSNV_eq + theme(legend.position = "None")



## color strips for acute and chronic infection

## make wrap label (top) in color
p1 = iSNV_eq+
  theme(strip.background = element_rect(fill="transparent"),
        legend.position = "None")+
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(strip.text = element_text(colour = 'white'))

dummy <- iSNV_eq+
  geom_rect(aes(fill=TVid), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  scale_fill_manual(values = eqHV_colors)+
  theme(legend.position = "None")



g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...)
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
strips <- grepl(pattern="strip_t", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips)
grid.newpage()
grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}
## ideally you'd remove the old strips, for now they're just covered
iSNV_eq_strip <- gtable_stack(new_strips,g1)
grid.newpage()
iSNV_eq2 <- grid.draw(iSNV_eq_strip)
## final plot
iSNV_eq2 <- as_ggplot(iSNV_eq_strip)




##### Figure: Entropy #####


## calc mean entropy per individuum ##
entropy_df_indiv <- entropy_df %>%
  group_by(Position, TVid, species) %>%
  summarise(Entropy=mean(as.numeric(`entropy(base e)`)))



## set levels for facet labels - Horse
entropy_df_levels_eq <- entropy_df_indiv %>%
  filter(species=="Nat")
entropy_df_levels_eq <- levels(factor(entropy_df_indiv$TVid))
facet_labs_entropy_eqMean <- c("Horse A", "Horse B", "Horse C", "Horse D", "Horse E",
                               "Horse F", "Horse G", "Horse H", "Horse I",
                               "Patient J", "Patient K", "Patient L", "Inoculum I", "Inoculum II")
names(facet_labs_entropy_eqMean) <- entropy_df_levels_eq


# plot entropy
entropy_df_indiv_eqHV <- entropy_df_indiv %>%
  filter(species=="Nat") %>%
  ggplot(aes(x=Position, y=Entropy, color=TVid))+
  geom_line()+
  scale_color_manual(values = eqHV_colors, label=c("Horse A", "Horse B", "Horse C", "Horse D", "Horse E"), name="Animals")+
  scale_y_continuous(breaks = seq(0, 1, .25), labels=seq(0, 1, .25), limits = c(0,1))+
  facet_wrap(~TVid, nrow=1, labeller = labeller(TVid = facet_labs_entropy_eqMean))+
  ylab("Entropy")+
  xlab("Position")+
  theme_classic()+
  theme(axis.text=element_text(size=10), axis.title=element_text(size=12),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)) )+
theme(legend.position="top")
entropy_df_indiv_eqHV

## legend
iSNV_rich_eqHV_legend <- get_legend(entropy_df_indiv_eqHV)

entropy_df_indiv_eqHV <- entropy_df_indiv_eqHV +
  theme(legend.position = "None")


# remove the strips for combined figure
entropy_df_indiv_eqHV2_nostrips <- entropy_df_indiv_eqHV+
  theme(strip.background = element_blank())+
  theme(plot.background = element_blank())+
  theme(strip.text = element_blank())









#### Figure: Divergence #####

#### threshold 3% line ####

dnds_df <- dnds_df %>% unite(TVid, time_point, Pos,col="ident" ,sep = "_",remove=FALSE)
entropy_df <- entropy_df %>% unite(TVid, time_point, Position,col="ident" ,sep = "_",remove=FALSE)
# join
dnds_entropy_df <- dnds_df %>% full_join(entropy_df, by=c("ident"="ident"),
                                            suffix = c("", ".y"))
# add levels
dnds_entropy_df$TVid <- factor(dnds_entropy_df$TVid, levels=c("horseA", "horseB", "horseC", "horseD", "horseE",
                                                                "horseF", "horseG", "horseH", "horseI",
                                                                "patientJ", "patientK", "patientL",
                                                                "InoculumI", "InoculumII"))



# add NS freq
dnds_entropy_df$obsNf <- dnds_entropy_df$obsN / dnds_entropy_df$Cov
dnds_entropy_df$obsSf <- dnds_entropy_df$obsS / dnds_entropy_df$Cov



## summerize per TVid
dnds_df_stats_richnessSNV <- dnds_entropy_df %>%
  filter(TVid %in% c("horseA", "horseB", "horseC", "horseD", "horseE")) %>%
  group_by(time_point, TVid, sample_date) %>%
  summarise(`Mean frequency`=mean(which(iSNV_freq>=iSNV_thres), na.rm = TRUE),
            richness=length(which(iSNV_freq>=iSNV_thres)), ratio=length(which(iSNV_freq>=iSNV_thres))/length(which(iSNV_freq>=iSNV_thres)),
            sum=sum(which(iSNV_freq>=iSNV_thres)),
            richnessLow=length(which(iSNV_freq>=lSNV)), sumLow=sum(which(iSNV_freq>=lSNV)),
            richnessConsensus=length(which(iSNV_freq>=hSNV)), sumConsensus=sum(which(iSNV_freq>=hSNV)),
            ratioConsensus=length(which(iSNV_freq>=hSNV))/length(which(iSNV_freq>=lSNV)),
            richnessN=length(which(obsNf>=iSNV_thres)),
            richnessS=length(which(obsSf>=iSNV_thres))
  )
dnds_df_stats_richnessSNV





### Divergence plot
richness_three <- dnds_df_stats_richnessSNV %>%
  ggplot(aes(group=TVid))+
  geom_line(aes(y=richness, x=time_point))+
  geom_point(aes(y=richness, x=time_point, size=ratio), alpha=.3)+
  geom_point(aes(y=richness, x=time_point, size=ratioConsensus), alpha=.6)+
  # add N S lines
  geom_line(aes(x=time_point, y= richnessN), color=N_color)+
  geom_point(aes(x=time_point, y= richnessN), color=N_color)+

  geom_line(aes(x=time_point, y= richnessS), color=S_color)+
  geom_point(aes(x=time_point, y= richnessS), color=S_color)+

  theme_classic()+
  scale_size_continuous(range = c(2, 8), name="Ratio")+
  scale_x_continuous(breaks = seq(1,6,1))+
  theme_classic()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))+
  theme(plot.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        strip.background = element_blank(), strip.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  labs(y = "Number iSNVs", x="Time point")+
  facet_grid(~TVid)
richness_three
#

ggsave("RichnessThree.pdf", width=24, height=6, units = "cm")

richness_three+
  theme(legend.position = "none" )




#### suppl figure: Divergence plot E1E2
## all - aka richness for 3% - SNV_thres
dnds_df_stats_richnessSNVE1E2 <- dnds_entropy_df %>%
  filter(TVid %in% c("horseA", "horseB", "horseC", "horseD", "horseE")) %>%
  group_by(time_point, TVid, sample_date, annotation_track) %>%
  summarise(`Mean frequency`=mean(which(iSNV_freq>=iSNV_thres), na.rm = TRUE),
            richness=length(which(iSNV_freq>=iSNV_thres)), ratio=length(which(iSNV_freq>=iSNV_thres))/length(which(iSNV_freq>=iSNV_thres)),
            sum=sum(which(iSNV_freq>=iSNV_thres)),
            richnessLow=length(which(iSNV_freq>=lSNV)), sumLow=sum(which(iSNV_freq>=lSNV)),
            richnessConsensus=length(which(iSNV_freq>=hSNV)), sumConsensus=sum(which(iSNV_freq>=hSNV)),
            ratioConsensus=length(which(iSNV_freq>=hSNV))/length(which(iSNV_freq>=lSNV)),
            richnessN=length(which(obsNf>=iSNV_thres)),
            richnessS=length(which(obsSf>=iSNV_thres))
  )
dnds_df_stats_richnessSNVE1E2




### three percent plot
richness_threeE1E2 <- dnds_df_stats_richnessSNVE1E2 %>%
  ggplot(aes(group=TVid))+
  geom_line(aes(y=richness, x=time_point))+
  geom_point(aes(y=richness, x=time_point, size=ratio), alpha=.3)+
  geom_point(aes(y=richness, x=time_point, size=ratioConsensus), alpha=.6)+
  # add N S lines
  geom_line(aes(x=time_point, y= richnessN), color=N_color)+
  geom_point(aes(x=time_point, y= richnessN), color=N_color)+

  geom_line(aes(x=time_point, y= richnessS), color=S_color)+
  geom_point(aes(x=time_point, y= richnessS), color=S_color)+

  theme_classic()+
  scale_size_continuous(range = c(2, 8), name="Ratio")+
  scale_x_continuous(breaks=seq(1,6,1))+
  theme_classic()+
  theme(axis.text=element_text(size=8), axis.title=element_text(size=10))+
  theme(plot.background = element_blank(),
        legend.background = element_rect(fill = "transparent"),
        panel.background = element_blank(),
        strip.background = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  labs(y = "Number SNVs", x="Time point")+
  facet_grid(annotation_track~TVid)
richness_threeE1E2
#

ggsave("RichnessThreeE1E2.pdf", width=24, height=10, units = "cm")

EqHVrichnessThreeE1E2 <- richness_threeE1E2+
  theme(legend.position = "none" )







##### Figure 1: Combine ####

## set axis parameter
grob_iSNV <- ggplotGrob(iSNV_eq2) # convert to gtable
grob_entropy <- ggplotGrob(entropy_df_indiv_eqHV2_nostrips) # convert to gtable

iSNV.widths <- grob_iSNV$widths[1:3] # extract the first three widths,
# corresponding to left margin, y lab, and y axis
entropy.widths <- grob_entropy$widths[1:3] # same for mpg plot


max.widths <- unit.pmax(iSNV.widths, entropy.widths) # calculate maximum widths

grob_iSNV$widths[1:3] <- max.widths # assign max. widths to gtable
grob_entropy$widths[1:3] <- max.widths # assign max widths to gtable





#### assemble figure 1.2
richness_three <- richness_three+
  theme(legend.position = "none")

entropy_richness_grid <- plot_grid(grob_entropy,richness_three,
                                   ncol=1, align = "v", axis="l", rel_heights = c(.65, 1))

axis <- align_plots(iSNV_eq2, entropy_richness_grid, align = "v", axis="l" )


EqHV_grid <- plot_grid(iSNV_rich_eqHV_legend, axis[[1]], axis[[2]], ncol=1, rel_heights = c(.25,.6,1))
EqHV_grid
ggsave("Figure1_Variant_plot_EqHV2_pp1.pdf", width=16, height=8)
ggsave("Figure1_Variant_plot_EqHV2_pp2.pdf", width=18, height=9)
ggsave("Figure1_Variant_plot_EqHV2_doc.pdf", width=16, height=10)
