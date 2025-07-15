#########
# SETUP #
#########

library(tidyverse)

#############################
# RECONSTRUCTING CPG STATES #
#############################

input_dir <- "./data/splices_bsbolt2bismark/"

# Extract sample file prefix (for later)
filename_SpliceA_OT <- list.files(path = input_dir, pattern = "*.SpliceA_OT.txt", recursive = TRUE)
sample_file_prefix <- sub("\\.SpliceA_OT\\.txt", "", filename_SpliceA_OT)

# Input methylation extract files
SpliceA_OT <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceA_OT.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceA_OB <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceA_OB.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceB_OT <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceB_OT.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceB_OB <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceB_OB.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceC_CTOT <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceC_CTOT.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceC_CTOB <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceC_CTOB.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceD_CTOT <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceD_CTOT.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()
SpliceD_CTOB <- paste0(input_dir, list.files(path = input_dir, pattern = "*.SpliceD_CTOB.txt", recursive = TRUE)) %>% read_tsv(col_names = FALSE) %>% as.data.frame()

# Add column names
MethColNames <- c("Cluster_ID", "Orientation", "Chr", "CpG_coordinate", "MethCall")
colnames(SpliceA_OT) <- MethColNames
colnames(SpliceA_OB) <- MethColNames
colnames(SpliceB_OT) <- MethColNames
colnames(SpliceB_OB) <- MethColNames
colnames(SpliceC_CTOT) <- MethColNames
colnames(SpliceC_CTOB) <- MethColNames
colnames(SpliceD_CTOT) <- MethColNames
colnames(SpliceD_CTOB) <- MethColNames

SpliceA_OT$Cluster_ID <- gsub("_.*","",SpliceA_OT$Cluster_ID)
SpliceA_OB$Cluster_ID <- gsub("_.*","",SpliceA_OB$Cluster_ID)
SpliceB_OT$Cluster_ID <- gsub("_.*","",SpliceB_OT$Cluster_ID)
SpliceB_OB$Cluster_ID <- gsub("_.*","",SpliceB_OB$Cluster_ID)
SpliceC_CTOT$Cluster_ID <- gsub("_.*","",SpliceC_CTOT$Cluster_ID)
SpliceC_CTOB$Cluster_ID <- gsub("_.*","",SpliceC_CTOB$Cluster_ID)
SpliceD_CTOT$Cluster_ID <- gsub("_.*","",SpliceD_CTOT$Cluster_ID)
SpliceD_CTOB$Cluster_ID <- gsub("_.*","",SpliceD_CTOB$Cluster_ID)

# Create columns for CpG_ID. Note that 4 nt were subtracted from coordinates due to inclusion of 'NNNN' at the start of the ref genome for H1/H2 (see above for explanation)
SpliceA_OT <- mutate(SpliceA_OT, CpG_ID = CpG_coordinate)
SpliceB_OB <- mutate(SpliceB_OB, CpG_ID = CpG_coordinate - 1)
SpliceC_CTOT <- mutate(SpliceC_CTOT, CpG_ID = CpG_coordinate)
SpliceD_CTOB <- mutate(SpliceD_CTOB, CpG_ID = CpG_coordinate - 1)

SpliceA_OB <- mutate(SpliceA_OB, CpG_ID = CpG_coordinate - 1)
SpliceB_OT <- mutate(SpliceB_OT, CpG_ID = CpG_coordinate)
SpliceC_CTOB <- mutate(SpliceC_CTOB, CpG_ID = CpG_coordinate - 1)
SpliceD_CTOT <- mutate(SpliceD_CTOT, CpG_ID = CpG_coordinate)

# Merge datasets corresponding to Splice A being OT
SpliceABCD <- list(SpliceA_OT,SpliceB_OB,SpliceC_CTOT,SpliceD_CTOB) %>% reduce(inner_join, by=c("Cluster_ID","Chr","CpG_ID"))

# Merge datasets corresponding to Splice A being OB. Note that the order in the list has changed -- essential for correct encoding of CpG states
SpliceABCD_v2 <- list(SpliceB_OT,SpliceA_OB,SpliceD_CTOT,SpliceC_CTOB) %>% reduce(inner_join, by=c("Cluster_ID","Chr","CpG_ID"))

# New column for CpG Site Meth call
SpliceABCD$CpG_MethCall <- paste(SpliceABCD$MethCall.x,
                              SpliceABCD$MethCall.y,
                              SpliceABCD$MethCall.x.x,
                              SpliceABCD$MethCall.y.y)

SpliceABCD_v2$CpG_MethCall <- paste(SpliceABCD_v2$MethCall.x,
                              SpliceABCD_v2$MethCall.y,
                              SpliceABCD_v2$MethCall.x.x,
                              SpliceABCD_v2$MethCall.y.y)

# Select only the relevant columns
key_info <- c("Cluster_ID", "Chr", "CpG_ID", "CpG_MethCall")
SpliceABCD_1 <- SpliceABCD[key_info]
SpliceABCD_2 <- SpliceABCD_v2[key_info]

# Merging both splice sets
SpliceABCD_merged <- rbind(SpliceABCD_1,SpliceABCD_2)

# Assign CpG-state and error labels
SpliceABCD_merged$CpG_State <-  gsub("Z Z Z Z","MM",
                                gsub("Z Z z z", "HH",
                                gsub("z z z z", "CC",
                                gsub("z Z z z", "CH",
                                gsub("Z z z z", "HC",
                                gsub("Z z z Z", "MC",
                                gsub("z Z Z z", "CM",
                                gsub("Z Z Z z", "HM",
                                gsub("Z Z z Z", "MH",
                                gsub("z z z Z", "Error",
                                gsub("z z Z Z", "Error",
                                gsub("z Z Z Z", "Error",
                                gsub("Z z Z z", "Error",
                                gsub("Z z Z Z", "Error",
                                gsub("z Z z Z", "Error",
                                gsub("z z Z z", "Error",
                                SpliceABCD_merged$CpG_MethCall))))))))))))))))


# Create output filename for merged splices
splices_output_filename <- paste0("./data/reconstruct_CpG_states/", sample_file_prefix, "_spliceABCD_merged.tsv")

# Write CpG states with read info to tsv file (comment out if not needed)
write_tsv(SpliceABCD_merged, file = splices_output_filename)

# Determine counts of each CpG state at each CpG site
CpG_state_quants <- SpliceABCD_merged %>%
  group_by(Chr, CpG_ID, CpG_State) %>%
  summarise(n = n()) %>%
  mutate(pct = n * 100 / sum(n))

# Create output filename for CpG state quants
CpG_state_quants_output_filename <- paste0("./data/reconstruct_CpG_states/", sample_file_prefix, "_CpG_state_quants.tsv")

# Write summary stats to tsv file
write_tsv(CpG_state_quants, file = CpG_state_quants_output_filename)

