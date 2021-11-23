# **************************************************************************
# *
# * Author:     
# *
# * Center for Biological Research Margarita Salas (CIB), CSIC 
# * Biocomputing Unit (BCU) of the National Center for Biotechnology (CNB), CSIC
# * Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
# *
# Short description: (FIRST DRAFT) This program creates files to for the variable selected from a database according to Nucleotide Combinations identified by NucleotideCombinations_detector.R.
# Wilcoxon rank sum test is applied.
# 
# **************************************************************************

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 5){
  stop("USE: Rscript Rcommander_fileMarker.R <statistical_data.txt> <BioChemMolCli.tsv> <ColumnSelectedTSVfile> <NumberOfSamples_CutOffPoint> <NameOf_ColumnSelectedTSVfile>")
}

# setwd("/home/eugenia/Documents/MIREIA_DATOS/RUN")
f_datos <- scan(args[1], what = character(), quiet = TRUE)
# f_datos <- scan("statistical_data.txt", what = character(), quiet = TRUE)
clinic_data <- read.csv(args[2], header = TRUE, sep = "\t", na.strings = NA)
# clinic_data <- read.csv("seqcovid_20210804.tsv", header = TRUE, sep = "\t", na.strings = NA)

colum_data <- as.numeric(args[3]) 
# colum_data <- 11 # this is for Ct/Cq.PCR.gen.S column in the TSV file

v_id <- c()
v_colum_data <- c()
for(value in 1:length(clinic_data[,as.numeric(args[3])])){
  if(is.na(clinic_data[value,as.numeric(args[3])]) == FALSE){
    # print(clinic_data[value,colum_data])
    if(as.numeric(clinic_data[value,as.numeric(args[3])]) != 0 && as.numeric(clinic_data[value,as.numeric(args[3])]) != 99){
      v_id <- c(v_id, clinic_data[value,1])
      v_colum_data <- c(v_colum_data, clinic_data[value,colum_data])
    }
  }
}

# Position vector for grupos starting with >
inicio_grupo <- c()
for(separador_seq in 1:length(f_datos)){
  if(f_datos[separador_seq] == ">"){
    # print(f_datos[separador_seq])
    inicio_grupo <- c(inicio_grupo, separador_seq)
  }
}
# print(inicio_grupo)

# Make table for R commander
table_st <- matrix(NA, nrow = length(v_id), ncol = 3)
# table_st <- matrix(NA, nrow = length(clinic_data[,1]), ncol = 3)

# Fill column 1 with ID
for(cd in 1:length(v_id)){
  sm_add <- paste(">", v_id[cd], sep = "")
  table_st[cd,1] <- sm_add 
  table_st[cd,2] <- v_colum_data[cd]
}

# Selection for colum_data with values >= 10
cutoff_point <- args[4]
# cutoff_point <- 10
vec_times <- c()
nuevo_ini <- c()
nuevo_fin <- c()
seq_selec <- c()
for(ig in 1:(length(inicio_grupo)-1)){
  if(as.numeric(f_datos[inicio_grupo[ig]+4]) >= as.numeric(args[4])){
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
  if(f_datos[last_posi_data] >= as.numeric(args[4])){
    nuevo_ini <- c(nuevo_ini, inicio_grupo[length(inicio_grupo)])
    nuevo_fin <- c(nuevo_fin, length(f_datos)) 
    seq_selec <- c(seq_selec, f_datos[inicio_grupo[length(inicio_grupo)]+1])
  }
}

v_pValue <- c()
# counter_files <- 1
for(nif in 1:length(nuevo_ini)){
  
  g_ini <- nuevo_ini[nif] + 5
  # g_ini <- nuevo_ini[nif] + 6  
  g_fin <- nuevo_fin[nif]
  
  for(cd in 1:length(v_id)){
    for(elemento in g_ini:g_fin){
      # if(elemento %% 2 == 0){
        if(table_st[cd,1] == f_datos[elemento]){
          table_st[cd,3] <- f_datos[elemento + 1]
          # print(f_datos[nuevo_ini[nif]+1])
          # print(paste(table_st[cd,1], table_st[cd,3], sep = " "))
        }
      # }
    }
  }
  colnames(table_st) <- c("ID", paste("All_", args[5], sep = ""), paste("Group_", args[5], sep = ""))
  # colnames(table_st) <- c("ID", "All Ct/Cq.PCR.gen.S", "Group Ct/Cq.PCR.gen.S")
  # View(table_st)

  # Make a file with table
  encabezado <- paste("ID", paste("All_", args[5], sep = ""), paste("Group_", args[5], sep = ""), sep = "\t")
  # encabezado <- paste("ID", "All Ct/Cq.PCR.gen.S", "Group Ct/Cq.PCR.gen.S", sep = "\t")
  write(encabezado, file = paste("file_RCommander", seq_selec[nif], ".txt", sep = ""), append = TRUE)
  for(cell in 1:nrow(table_st)){
    linea <- paste(table_st[cell,1], table_st[cell,2], table_st[cell,3], sep = "\t")
    write(linea, file = paste("file_RCommander", seq_selec[nif], ".txt", sep = ""), append = TRUE)
  }
  # counter_files <- counter_files + 1
  sink(paste("file_RCommander", seq_selec[nif], ".txt", sep = ""), append = TRUE)
  # Wilcox test
  print(wilcox.test(x = as.numeric(table_st[,2]), y = as.numeric(table_st[,3]), alternative = "greater", paired = FALSE, conf.int = 0.95))
  sink()
  #
  sink("w_test.txt")
  print(wilcox.test(x = as.numeric(table_st[,2]), y = as.numeric(table_st[,3]), alternative = "greater", paired = FALSE, conf.int = 0.95))
  sink()
  #
  # Read statistic test
  wtest <- scan("w_test.txt", what = character(), quiet = TRUE)
  
  for(wt in 1:length(wtest)){
    if(wtest[wt] == "p-value"){
      # print(wtest[wt+2])
      v_pValue <- c(v_pValue, format(as.numeric(wtest[wt+2]), scientific = FALSE))
    }
  }
  # Clean column 3 
  for(quitar in 1:length(v_id)){
    table_st[quitar,3] <- NA
  }
}
unlink("w_test.txt")
df_seq_pvalue <- data.frame("NucleotideCombination"=seq_selec, "n"=vec_times, "p-value"=v_pValue)
# print(df_seq_pvalue)

# Make a data.frame with column FDR 
df_seq_pvalue_amp <- cbind(df_seq_pvalue, "FDR")

# Order df_seq_pvalue_amp
order_df_seq_pvalue_amp <-df_seq_pvalue_amp[
  with(df_seq_pvalue_amp, order(df_seq_pvalue_amp[,3])),
]
# print(order_df_seq_pvalue_amp)

# Calculate FDR (this means, after order p values from smallest to greatest multiply the first by one, second by two, third by three and so on)
veces <- 0
for(dfa in 1:length(df_seq_pvalue_amp[,3])){
  veces <- veces + 1
  order_df_seq_pvalue_amp[dfa,4] <- format(as.numeric(order_df_seq_pvalue_amp[dfa,3]) * veces, scientific = FALSE)
}
print(order_df_seq_pvalue_amp)
