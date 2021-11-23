# **************************************************************************
# *
# * Author:     
# *
# * Center for Biological Research Margarita Salas (CIB), CSIC 
# * Biocomputing Unit (BCU) of the National Center for Biotechnology (CNB), CSIC
# * Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
# *
# Short description: (FIRST DRAFT) This program makes a boxplot for Nucleotides Combinations that do not have the number of samples required for statistical analysis.
# 
# **************************************************************************

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
  stop("USE: Rscript noSampleSizeCondition_NucleotideCombinations.R <statistical_data.txt> <BioChemMolCli.tsv> <NumberOfSamples_CutOffPoint> <NameOf_ColumnSelectedTSVfile>")
}

# setwd("/home/eugenia/Documents/MIREIA_DATOS/RUN")
f_datos <- scan(args[1], what = character(), quiet = TRUE)
# f_datos <- scan("statistical_data.txt", what = character(), quiet = TRUE)
clinic_data <- read.csv(args[2], header = TRUE, sep = "\t", na.strings = NA)
# clinic_data <- read.csv("seqcovid_20210804.tsv", header = TRUE, sep = "\t", na.strings = NA)

# colum_data <- 11 # this is for Ct/Cq.PCR.gen.S column in the TSV file

# Position vector for grupos starting with >
inicio_grupo <- c()
for(separador_seq in 1:length(f_datos)){
  if(f_datos[separador_seq] == ">"){
    # print(f_datos[separador_seq])
    inicio_grupo <- c(inicio_grupo, separador_seq)
  }
}
# print(inicio_grupo)

# Selection for colum_data with values >= 10
cutoff_point <- as.numeric(args[3])
# cutoff_point <- 10
vec_times <- c()
nuevo_ini <- c()
nuevo_fin <- c()
seq_selec <- c()
for(ig in 1:(length(inicio_grupo)-1)){
  if(as.numeric(f_datos[inicio_grupo[ig]+4]) >= as.numeric(args[3])){
    # print(f_datos[inicio_grupo[ig]+1])
    vec_times <- c(vec_times, as.numeric(f_datos[inicio_grupo[ig]+4]))
    seq_selec <- c(seq_selec, f_datos[inicio_grupo[ig]+1])
    nuevo_ini <- c(nuevo_ini, inicio_grupo[ig])
    nuevo_fin <- c(nuevo_fin, inicio_grupo[ig+1])
  }
}
# print(seq_selec)

if(inicio_grupo[length(inicio_grupo)] == nuevo_ini[length(nuevo_ini)]){
  # For last position 
  last_posi <- inicio_grupo[length(inicio_grupo)]
  last_posi_data <- inicio_grupo[length(inicio_grupo)]+4
  if(f_datos[last_posi_data] >= as.numeric(args[3])){
    nuevo_ini <- c(nuevo_ini, inicio_grupo[length(inicio_grupo)])
    nuevo_fin <- c(nuevo_fin, length(f_datos)) 
    seq_selec <- c(seq_selec, f_datos[inicio_grupo[length(inicio_grupo)]+1])
  }
}
# print(seq_selec)

# ">" position in a vector
vec_solo_sm_pos <- c()
for(ele_f_datos in 1:length(f_datos)){
  if(f_datos[ele_f_datos] == ">"){
    vec_solo_sm_pos <- c(vec_solo_sm_pos, ele_f_datos)
  }
}
# print(vec_solo_sm_pos)

# Make 2 vectors for the start and end of the corresponding Sequence_
vec_start_Seq <- setdiff(vec_solo_sm_pos, vec_solo_sm_pos[length(vec_solo_sm_pos)])

vec_end_Seq <- c()
for(ele_vec_sm in 2:length(vec_solo_sm_pos)){
  vec_end_Seq <- c(vec_end_Seq, (vec_solo_sm_pos[ele_vec_sm] - 1))
}
# print(vec_end_Seq)

# Make a file 

encabezado <- paste("ID", paste("No_FulfillCondition_N.Combinations_", args[4], sep = ""), sep = "\t")
write(encabezado, file = "All_N.Combinations_NoSampleCondition.txt", append = TRUE)

# Number of intervals
l_interval <- length(vec_end_Seq)
ite_seq_selec <- 1
for(ele_seq in 1:l_interval){
  for(ele_interval in vec_start_Seq[ele_seq]:vec_end_Seq[ele_seq]){
    if(f_datos[vec_start_Seq[ele_seq] + 1] != seq_selec[ite_seq_selec]){
      detect <- 0
      for(ele_se_selec in 1:length(seq_selec)){
        if(f_datos[vec_start_Seq[ele_seq] + 1] == seq_selec[ele_se_selec]){
          detect <- detect + 1
        }
      }
      if(detect != 1){
        if(grepl(">COV",f_datos[ele_interval]) == TRUE){
          write(paste(f_datos[ele_interval],f_datos[ele_interval + 1], sep = "\t"), file = "All_N.Combinations_NoSampleCondition.txt", append = TRUE)
        }
        
        # write(f_datos[ele_interval], file = "All_Sequence_Times_SmallerThan.txt", append = TRUE)
      }
    }
  }
  # print(detect)
}

# For the last ">" 
dtect_2 <- 0
for(ele_last_sm in vec_solo_sm_pos[length(vec_solo_sm_pos)]:length(f_datos)){
  if(f_datos[vec_solo_sm_pos[length(vec_solo_sm_pos)]+1] != seq_selec[ite_seq_selec]){
    detect_2 <- 0
    for(ele_se_selec in 1:length(seq_selec)){
      if(f_datos[vec_start_Seq[ele_seq] + 1] == seq_selec[ele_se_selec]){
        detect_2 <- detect_2 + 1
      }
    }
    if(detect_2 != 1){
      if(grepl(">COV",f_datos[ele_interval]) == TRUE){
      write(f_datos[ele_interval], file = "All_N.Combinations_NoSampleCondition.txt", append = TRUE)
      }
    }
  }
}

# Make vector with Ct/Cq values
f_datos_no_condi <- scan("All_N.Combinations_NoSampleCondition.txt", what = character(), quiet = TRUE)

id_cov_pos <- grep(">COV", f_datos_no_condi)

ct_no_condi <- c()
for(ele_id_cov in 1:length(id_cov_pos)){
  ct_no_condi <- c(ct_no_condi, as.numeric(f_datos_no_condi[id_cov_pos[ele_id_cov] + 1]))
}
# print(ct_no_condi)

# Boxplot
jpeg(paste('Rplots_BoxPlots_No_Fulfill_Condition', '.jpg', sep = ""))
hist(ct_no_condi, probability = TRUE, xlab = paste("No_Fulfill_Condition_Group_", args[4], sep = ""), ylab = "", col = "grey", axes = FALSE, main = "")
lines(density(ct_no_condi), col = "red", lwd = 2)
par(new = TRUE)
boxplot(ct_no_condi,  col = rgb(0, 1, 0, alpha = 0.15))
points(mean(ct_no_condi, na.rm=TRUE), col = "red", pch = 19)
dev.off()


