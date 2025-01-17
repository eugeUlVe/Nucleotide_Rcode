# **************************************************************************
# *
# * Author:     
# *
# * Center for Biological Research Margarita Salas (CIB), CSIC 
# * Biocomputing Unit (BCU) of the National Center for Biotechnology (CNB), CSIC
# * Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
# *
# Short description: (FIRST DRAFT) This program makes boxplots for files generated by Rcommander_fileMarker.R
# 
# **************************************************************************

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 0){
  stop("USE: Rscript boxplots_NucleotideCombinations.R")
}

# Read files and read 
file_list <- list.files(path=getwd())

root_name <- "file_RCommanderNucleotideCombination_"

for(ele_list in 1:length(file_list)){
  if(grepl(root_name,file_list[ele_list], fixed = TRUE)){
    # print(file_list[ele_list])
    f_file <- scan(file_list[ele_list], what = character(), quiet = TRUE)
    
    #
    split_f_file <- strsplit(file_list[ele_list], "")
    # print(split_f_file)
    l_file_list <- length(split_f_file[[1]])
    
    # Due to the number of sequence might be different to the unit this has to be taken to account for the name length. NOTE: if file name root is changed fix the length 
    if(l_file_list > 29){
      dif_log <- l_file_list - 29
      dif_ini <- 13 + dif_log
      ini_pal_seq <- l_file_list - dif_ini
    }
    else
      ini_pal_seq <- l_file_list - 13
    # Paste from NucleotideCombination_[...]
    seq_name <- ""
    for(ele_file_list in ini_pal_seq:l_file_list){
      seq_name <- paste(seq_name, split_f_file[[1]][ele_file_list], sep = "")
    }
    # print(seq_name)
    #
    
    # f_file in array way
    vec_f_file <- c()
    for(ele_file in 1:length(f_file)){
      vec_f_file <- c(vec_f_file, f_file[ele_file])
    }
    # print(vec_f_file)
    
    # Put position for ID, All and group in vectors
    id_seq_pos <- grep(">", vec_f_file)
    
    vec_all <- c()
    vec_group <- c()
    for(ele_id in 1:length(id_seq_pos)){
      vec_all <- c(vec_all, as.numeric(f_file[id_seq_pos[ele_id] + 1]))
      vec_group <- c(vec_group, as.numeric(f_file[id_seq_pos[ele_id] + 2]))
    }
    # Make data frame to plot both boxplots
    df_set <- data.frame(All=vec_all, Group=vec_group)
    
    # BoxPlot 
    # 2 figures arranged in 2 rows and 1 column
    jpeg(paste('Rplots_BoxPlots_' , seq_name, '.jpg', sep = ""))
    par(mfrow=c(2,1))
    # All
    hist(df_set[,1], probability = TRUE, xlab = "All", ylab = "", col = "grey", axes = FALSE, main = "")
    lines(density(df_set[,1]), col = "red", lwd = 2)
    par(new = TRUE)
    boxplot(df_set[,1],  col = rgb(1, 0, 1, alpha = 0.15))
    points(mean(df_set[,1], na.rm=TRUE), col = "red", pch = 19)
    # Group
    hist(df_set[,2], probability = TRUE, xlab = "Group", ylab = "", col = "grey", axes = FALSE, main = "")
    lines(density(df_set[,1]), col = "blue", lwd = 2)
    par(new = TRUE)
    boxplot(df_set[,2],  col = rgb(0, 0, 1, alpha = 0.15))
    points(mean(df_set[,2], na.rm=TRUE), col = "blue", pch = 19)
    dev.off()
  }
}
print("Done")