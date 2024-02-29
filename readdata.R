##### read tables #####


#### library #####
library(tidyverse)



##### in ####
data_in <- "/Users/andregomer/tempory_data/E1E2_revision/data/E1E2data"





##### out ####
out <- "/Users/andregomer/tempory_data/E1E2_revision/E1E2data_out"
setwd(out)




#### read data ####
df <- data.frame()
entropy_df <- data.frame()
dnds_df <- data.frame()
AA_df <- data.frame()
coord_s <- 100


for (i in list.files(data_in)) {
  for (j in list.files(paste(data_in, i, sep = "/"))) {
    data_path <- paste(data_in, i, j, sep = "/")

    df_entropy <- read_delim(paste(paste(data_path, j, sep="/"),".toref_dedup_E1E2_entropy.txt", sep=""),
                     delim = "\t")

    df_dnds <- read_delim(paste(paste(data_path, j, sep="/"),".mapped.sorted_dedup_DNDS.txt", sep=""),
                     delim = "\t")

    df_AA <- read_delim(paste(paste(data_path, j, sep="/"),".toref_dedup_E1E2_AA.txt", sep=""),
                          delim = "\t")

    ### entropy
    df_entropy <-  df_entropy %>% filter(Position >= coord_s) # sel pos >99
    df_entropy$sample_id <- paste(unlist(str_split(j, "_"))[1], unlist(str_split(j, "_"))[2], sep=" ")
    df_entropy$sample_date <- unlist(str_split(j, "_"))[2]
    df_entropy$TVid <- unlist(str_split(j, "_"))[1]
    df_entropy$species <- i
    entropy_df <- rbind(entropy_df, df_entropy)


    ### dnds
    df_dnds
    df_dnds$sample_id <- paste(unlist(str_split(j, "_"))[1], unlist(str_split(j, "_"))[2], sep=" ")
    df_dnds$sample_date <- unlist(str_split(j, "_"))[2]
    df_dnds$TVid <- unlist(str_split(j, "_"))[1]
    df_dnds$species <- i
    dnds_df <- rbind(dnds_df, df_dnds)



    dnds_df <- rbind(dnds_df, df)


    ### AA
    df_AA$sample_id <- paste(unlist(str_split(j, "_"))[1], unlist(str_split(j, "_"))[2], sep=" ")
    df_AA$sample_date <- unlist(str_split(j, "_"))[2]
    df_AA$TVid <- unlist(str_split(j, "_"))[1]
    df_AA$species <- i
    AA_df <- rbind(AA_df, df_AA)

  }
}





###### add entropy data ######
entropy_df$iSNV_freq <- as.numeric(entropy_df$NonRefCnt) / entropy_df$Coverage


#### nat

# positions
entropy_df <- entropy_df %>% filter(((species == "Nat" | species == "TV") & Position <= 1668) |
                                      (TVid=="patientJ" & Position <= 1752) |
                                      (TVid=="patientK" & Position <= 1761) |
                                      TVid=="patientL" & Position <= 1749)


# annotation_track
entropy_df$annotation_track <- NA
entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, (species == "Nat" | species == "TV") &
                                                               Position >= 100 & Position <= 660, "E1"))
entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, (species == "Nat" | species == "TV") &
                                                               Position >= 661 & Position <= 1668, "E2"))

entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientJ" &
                                                               Position >= 100 & Position <= 675, "E1"))
entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientJ" &
                                                               Position >= 676 & Position <= 1761, "E2"))

entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientK" &
                                                               Position >= 100 & Position <= 675, "E1"))
entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientK" &
                                                               Position >= 676 & Position <= 1752, "E2"))

entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientL" &
                                                               Position >= 100 & Position <= 675, "E1"))
entropy_df <- entropy_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientL" &
                                                               Position >= 676 & Position <= 1749, "E2"))
# ORF_size
entropy_df$orf_size <- NA
entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                                               Position >= 100 & Position <= 660, 561))
entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                                               Position >= 661 & Position <= 1668, 1008))

entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & annotation_track=="E1", 576))
entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & annotation_track=="E2", 1086))

entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & annotation_track=="E1", 576))
entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & annotation_track=="E2", 1077))

entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & annotation_track=="E1", 576))
entropy_df <- entropy_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & annotation_track=="E2", 1074))



# time point
entropy_df$time_point <- NA

## add tp
# horseA
entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                         sample_id=="horseA week01" |
                                                         sample_id=="horseB week01" |
                                                         sample_id=="horseC week01" |
                                                         sample_id=="horseD week01" |
                                                         sample_id=="horseE week01" |
                                                         sample_id=="horseF week02" |
                                                         sample_id=="horseG week01" |
                                                         sample_id=="horseH week01" |
                                                         sample_id=="horseI week01" |
                                                         sample_id=="patientJ week01" |
                                                         sample_id=="patientK week01" |
                                                         sample_id=="patientL week01", 1))

entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseA week04" |
                                                         sample_id=="horseB week07" |
                                                         sample_id=="horseC week03" |
                                                         sample_id=="horseD week08" |
                                                         sample_id=="horseE week08" |
                                                         sample_id=="horseF week04" |
                                                         sample_id=="horseG week03" |
                                                         sample_id=="horseH week03" |
                                                         sample_id=="horseI week02" |
                                                         sample_id=="patientJ week29" |
                                                         sample_id=="patientK week105" |
                                                         sample_id=="patientL week04", 2))


entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseB week08" |
                                                         sample_id=="horseC week04" |
                                                         sample_id=="horseD week28" |
                                                         sample_id=="horseE week19" |
                                                         sample_id=="horseF week06" |
                                                         sample_id=="horseG week05" |
                                                         sample_id=="horseH week05" |
                                                         sample_id=="patientJ week45" |
                                                         sample_id=="patientK week242" |
                                                         sample_id=="patientL week116", 3))


entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                         sample_id=="horseC week28" |
                                                         sample_id=="horseD week54" |
                                                         sample_id=="horseE week37" |
                                                         sample_id=="horseF week08" |
                                                         sample_id=="horseG week07" |
                                                         sample_id=="horseH week06" |
                                                         sample_id=="patientL week186", 4))


entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                         sample_id=="horseC week54" |
                                                         sample_id=="horseD week82" |
                                                         sample_id=="horseE week94" |
                                                         sample_id=="horseF week10" |
                                                         sample_id=="horseG week08" |
                                                         sample_id=="horseH week07" |
                                                         sample_id=="patientL week225", 5))

entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseC week82" |
                                                         sample_id=="horseF week11" |
                                                         sample_id=="patientL week295", 6))



entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                         sample_id=="horseF week12" |
                                                         sample_id=="patientL week319", 7))

entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                             TVid=="inoculumI", 8))

entropy_df <- entropy_df %>% mutate(time_point=replace(time_point,
                                                       TVid=="inoculumII", 9))


# infection status
entropy_df$status <- "acute"
entropy_df <- entropy_df %>% mutate(status=replace(status, TVid %in% c("horseC", "horseD", "horseE", "inoculumI", "inoculumII") , "chronic"))
entropy_df <- entropy_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicNT"))
entropy_df <- entropy_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicT"))


write_csv(entropy_df, "entropy_df.csv")




###### dnds ######

# annotation_track
dnds_df$annotation_track <- NA
dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, (species == "Nat" | species == "TV") &
                                                               Pos >= 100 & Pos <= 660, "E1"))
dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, (species == "Nat" | species == "TV") &
                                                               Pos >= 661 & Pos <= 1668, "E2"))

dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientJ" &
                                                               Pos >= 100 & Pos <= 675, "E1"))
dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientJ" &
                                                               Pos >= 676 & Pos <= 1761, "E2"))

dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientK" &
                                                               Pos >= 100 & Pos <= 675, "E1"))
dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientK" &
                                                               Pos >= 676 & Pos <= 1752, "E2"))

dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientL" &
                                                               Pos >= 100 & Pos <= 675, "E1"))
dnds_df <- dnds_df %>% mutate(annotation_track=replace(annotation_track, TVid=="patientL" &
                                                               Pos >= 676 & Pos <= 1749, "E2"))
# ORF_size
dnds_df$orf_size <- NA
dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                                       Pos >= 100 & Pos <= 660, 561))
dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                                       Pos >= 661 & Pos <= 1668, 1008))

dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & annotation_track=="E1", 576))
dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & annotation_track=="E2", 1086))

dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & annotation_track=="E1", 576))
dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & annotation_track=="E2", 1077))

dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & annotation_track=="E1", 576))
dnds_df <- dnds_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & annotation_track=="E2", 1074))



# time point
dnds_df$time_point <- NA

## add tp
# horseA
dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseA week01" |
                                                         sample_id=="horseB week01" |
                                                         sample_id=="horseC week01" |
                                                         sample_id=="horseD week01" |
                                                         sample_id=="horseE week01" |
                                                         sample_id=="horseF week02" |
                                                         sample_id=="horseG week01" |
                                                         sample_id=="horseH week01" |
                                                         sample_id=="horseI week01" |
                                                         sample_id=="patientJ week01" |
                                                         sample_id=="patientK week01" |
                                                         sample_id=="patientL week01", 1))

dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseA week04" |
                                                         sample_id=="horseB week07" |
                                                         sample_id=="horseC week03" |
                                                         sample_id=="horseD week08" |
                                                         sample_id=="horseE week08" |
                                                         sample_id=="horseF week04" |
                                                         sample_id=="horseG week03" |
                                                         sample_id=="horseH week03" |
                                                         sample_id=="horseI week02" |
                                                         sample_id=="patientJ week29" |
                                                         sample_id=="patientK week105" |
                                                         sample_id=="patientL week04", 2))


dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseB week08" |
                                                       sample_id=="horseC week04" |
                                                         sample_id=="horseD week28" |
                                                         sample_id=="horseE week19" |
                                                         sample_id=="horseF week06" |
                                                         sample_id=="horseG week05" |
                                                         sample_id=="horseH week05" |
                                                         sample_id=="patientJ week45" |
                                                         sample_id=="patientK week242" |
                                                         sample_id=="patientL week116", 3))


dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseC week28" |
                                                         sample_id=="horseD week54" |
                                                         sample_id=="horseE week37" |
                                                         sample_id=="horseF week08" |
                                                         sample_id=="horseG week07" |
                                                         sample_id=="horseH week06" |
                                                         sample_id=="patientL week186", 4))


dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseC week54" |
                                                         sample_id=="horseD week82" |
                                                         sample_id=="horseE week94" |
                                                         sample_id=="horseF week10" |
                                                         sample_id=="horseG week08" |
                                                         sample_id=="horseH week07" |
                                                         sample_id=="patientL week225", 5))

dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseC week82" |
                                                       sample_id=="horseF week11" |
                                                         sample_id=="patientL week295", 6))



dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                       sample_id=="horseF week12" |
                                                         sample_id=="patientL week319", 7))

dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                 TVid=="inoculumI", 8))

dnds_df <- dnds_df %>% mutate(time_point=replace(time_point,
                                                 TVid=="inoculumII", 9))

# infection status
dnds_df$status <- "acute"
dnds_df <- dnds_df %>% mutate(status=replace(status, TVid %in% c("horseC", "horseD", "horseE") , "chronic"))
dnds_df <- dnds_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicNT"))
dnds_df <- dnds_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicT"))


write_csv(dnds_df, "dnds_df.csv")
















###### AA ######

# ORF_size
AA_df$orf_size <- NA
AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                             Protein=="E1", 561))
AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, (species == "Nat" | species == "TV") &
                                             Protein=="E1", 1008))

AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & Protein=="E1", 576))
AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientJ" & Protein=="E2", 1086))

AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & Protein=="E1", 576))
AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientK" & Protein=="E2", 1077))

AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & Protein=="E1", 576))
AA_df <- AA_df %>% mutate(orf_size=replace(orf_size, TVid=="patientL" & Protein=="E2", 1074))



# time point
AA_df$time_point <- NA

## add tp
# horseA
AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseA week01" |
                                                   sample_id=="horseB week01" |
                                                   sample_id=="horseC week01" |
                                                   sample_id=="horseD week01" |
                                                   sample_id=="horseE week01" |
                                                   sample_id=="horseF week02" |
                                                   sample_id=="horseG week01" |
                                                   sample_id=="horseH week01" |
                                                   sample_id=="horseI week01" |
                                                   sample_id=="patientJ week01" |
                                                   sample_id=="patientK week01" |
                                                   sample_id=="patientL week01", 1))

AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseA week04" |
                                                   sample_id=="horseB week07" |
                                                   sample_id=="horseC week03" |
                                                   sample_id=="horseD week08" |
                                                   sample_id=="horseE week08" |
                                                   sample_id=="horseF_week04" |
                                                   sample_id=="horseG week03" |
                                                   sample_id=="horseH week03" |
                                                   sample_id=="horseI week02" |
                                                   sample_id=="patientJ week29" |
                                                   sample_id=="patientK week105" |
                                                   sample_id=="patientL week04", 2))


AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                             sample_id=="horseB week08" |
                                                 sample_id=="horseC week04" |
                                                   sample_id=="horseD week28" |
                                                   sample_id=="horseE week19" |
                                                   sample_id=="horseF_week06" |
                                                   sample_id=="horseG week05" |
                                                   sample_id=="horseH week05" |
                                                   sample_id=="patientJ week45" |
                                                   sample_id=="patientK week242" |
                                                   sample_id=="patientL week116", 3))


AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseC week28" |
                                                   sample_id=="horseD week54" |
                                                   sample_id=="horseE week37" |
                                                   sample_id=="horseF_week08" |
                                                   sample_id=="horseG week07" |
                                                   sample_id=="horseH week06" |
                                                   sample_id=="patientL week186", 4))


AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseC week54" |
                                                   sample_id=="horseD week82" |
                                                   sample_id=="horseE week94" |
                                                   sample_id=="horseF_week10" |
                                                   sample_id=="horseG week08" |
                                                   sample_id=="horseH week07" |
                                                   sample_id=="patientL week225", 5))

AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                             sample_id=="horseC week82" |
                                                 sample_id=="horseF_week11" |
                                                   sample_id=="patientL week295", 6))



AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                                 sample_id=="horseF_week12" |
                                                   sample_id=="patientL week319", 7))

AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                             TVid=="inoculumI", 8))

AA_df <- AA_df %>% mutate(time_point=replace(time_point,
                                             TVid=="inoculumII", 9))




# infection status
AA_df$status <- "acute"
AA_df <- AA_df %>% mutate(status=replace(status, TVid %in% c("horseC", "horseD", "horseE") , "chronic"))
AA_df <- AA_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicNT"))
AA_df <- AA_df %>% mutate(status=replace(status, TVid %in% c("patientL") , "chronicT"))


write_csv(AA_df, "AA_df.csv")
















##### codon AA-table ####

## input dnds table
data <- dnds_df

data_codon_table_collapsed <- data.frame(
  Pos=data %>% dplyr::select(Pos),
  Ref=data %>% dplyr::select(Ref),
  RefCodon=data %>% dplyr::select(RefCodon),
  RefAA=data %>% dplyr::select(RefAA),
  dNdS=data %>% dplyr::select(dNdS),
  dNdSslide=data %>% dplyr::select(`Sliding(3)dNdS`),
  time_point=data %>% dplyr::select(time_point),
  sample_id=data %>% dplyr::select(sample_id),
  orf_size=data %>% dplyr::select(orf_size),
  annotation_track=data %>% dplyr::select(annotation_track),
  TVid=data %>% dplyr::select(TVid),
  status=data %>% dplyr::select(status),
  species=data %>% dplyr::select(species),
  Coverage=data %>% dplyr::select(CodonCov),
  sample_date=data %>% dplyr::select(sample_date),
  A=data %>%
    dplyr::select(dplyr::ends_with("A") & dplyr::contains("-")) %>%
    rowSums(),
  C=data %>%
    dplyr::select(dplyr::ends_with("C") & dplyr::contains("-")) %>%
    rowSums(),
  D=data %>%
    dplyr::select(dplyr::ends_with("D") & dplyr::contains("-")) %>%
    rowSums(),
  E=data %>%
    dplyr::select(dplyr::ends_with("E") & dplyr::contains("-")) %>%
    rowSums(),
  F=data %>%
    dplyr::select(dplyr::ends_with("F") & dplyr::contains("-")) %>%
    rowSums(),
  G=data %>%
    dplyr::select(dplyr::ends_with("G") & dplyr::contains("-")) %>%
    rowSums(),
  H=data %>%
    dplyr::select(dplyr::ends_with("H") & dplyr::contains("-")) %>%
    rowSums(),
  I=data %>%
    dplyr::select(dplyr::ends_with("I") & dplyr::contains("-")) %>%
    rowSums(),
  K=data %>%
    dplyr::select(dplyr::ends_with("K") & dplyr::contains("-")) %>%
    rowSums(),
  L=data %>%
    dplyr::select(dplyr::ends_with("L") & dplyr::contains("-")) %>%
    rowSums(),
  M=data %>%
    dplyr::select(dplyr::ends_with("M") & dplyr::contains("-")) %>%
    rowSums(),
  N=data %>%
    dplyr::select(dplyr::ends_with("N") & dplyr::contains("-")) %>%
    rowSums(),
  P=data %>%
    dplyr::select(dplyr::ends_with("P") & dplyr::contains("-")) %>%
    rowSums(),
  Q=data %>%
    dplyr::select(dplyr::ends_with("Q") & dplyr::contains("-")) %>%
    rowSums(),
  R=data %>%
    dplyr::select(dplyr::ends_with("R") & dplyr::contains("-")) %>%
    rowSums(),
  S=data %>%
    dplyr::select(dplyr::ends_with("S") & dplyr::contains("-")) %>%
    rowSums(),
  T=data %>%
    dplyr::select(dplyr::ends_with("T") & dplyr::contains("-")) %>%
    rowSums(),
  V=data %>%
    dplyr::select(dplyr::ends_with("V") & dplyr::contains("-")) %>%
    rowSums(),
  W=data %>%
    dplyr::select(dplyr::ends_with("W") & dplyr::contains("-")) %>%
    rowSums(),
  Y=data %>%
    dplyr::select(dplyr::ends_with("Y") & dplyr::contains("-")) %>%
    rowSums(),
  X=data %>%
    dplyr::select(dplyr::ends_with("X") & dplyr::contains("-")) %>%
    rowSums(),
  NX=data %>%
    dplyr::select(dplyr::ends_with("?") & dplyr::contains("-")) %>%
    rowSums())




data_codon_table_collapsed <- data_codon_table_collapsed[seq(3, nrow(data_codon_table_collapsed), 3),] # filter 3rd position

data_codon_table_collapsed$Pos <- (data_codon_table_collapsed$Pos-99) / 3

###### calc reference coverage ######
data_codon_table_collapsed["Ref_counts"] <- NA
AAcolnames <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y","X", "NX")
for (i in 1:nrow(data_codon_table_collapsed)) {
  for (j in AAcolnames) {
    if(j == data_codon_table_collapsed$RefAA[i]){
      data_codon_table_collapsed$Ref_counts[i] <- data_codon_table_collapsed[i,j]
      print(i)
      print(j)
    }
  }
}

# calc frequencies for consensus and variants
data_codon_table_collapsed["consensus_freq"] <- as.numeric(data_codon_table_collapsed$Ref_counts) / data_codon_table_collapsed$CodonCov
data_codon_table_collapsed["variant_freq"] <- 1 - data_codon_table_collapsed$consensus_freq


write_csv(data_codon_table_collapsed, "AAcodonTable.csv")
