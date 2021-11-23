# **************************************************************************
# *
# * Authors:     
# *
# * Center for Biological Research Margarita Salas (CIB), CSIC 
# * Biocomputing Unit (BCU) of the National Center for Biotechnology (CNB), CSIC
# * Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
# *
# Short description: (FIRST DRAFT) This program makes fasta file for those Id's with query column in the dataset with values different to NA, 0 and 99.
# 
# **************************************************************************

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 4){
  stop("USE: Rscript seqcovid_sequenceFiltered_for_clinicalData.R <BioChemMolCli.tsv> <ColumnSelectedTSVfile> <sequences.fasta> <Output_DNA_sequences.fasta>")
}

data <- read.csv(args[1], header = TRUE, sep = "\t")
colum_data <- as.numeric(args[2])

# setwd("/home/eugenia/Documents/MIREIA_DATOS/seqcovid")
# data <- read.csv("seqcovid_20210804.tsv", header = TRUE, sep = "\t")
# colum_data <- 11 # this is for Ct/Cq.PCR.gen.S

v_id <- c()
for(id in 1:length(data[,1])){
  v_id <- c(v_id, data[id,1])
}
# print(v_id)

data_id_filter <- unique(data[,1])

sample_row <- c()
for(fila in 1:length(data$Ct.Cq.PCR.gen.S)){
  if(is.na(data[fila,as.numeric(args[2])]) == FALSE){
    if(data[fila,as.numeric(args[2])] != 0 && data[fila,as.numeric(args[2])] != 99){
    # print(fila)
    sample_row <- c(sample_row, fila)
    }
  }
}

samples_ID <- c()
for(s in 1:length(sample_row)){
  samples_ID <- c(samples_ID, data[sample_row[s],1])
}

# setwd("/home/eugenia/Documents/MIREIA_DATOS/seqcovid")
secuencias <- scan(args[3], what = character(), quiet = TRUE)
# secuencias <- scan("seqcovid_20210804.fasta", what = character(), quiet = TRUE)

for(cov in 1:length(samples_ID)){
  sm_sample_ID <- paste(">", samples_ID[cov], sep = "")
  for(line in 1:(length(secuencias)-1)){
    if(sm_sample_ID == secuencias[line]){
      write(secuencias[line], file = args[4], append = TRUE)
      write(secuencias[line + 1], file = args[4], append = TRUE)
      # write(secuencias[line], file = "available_sequence_geneS_filtrado.fasta", append = TRUE)
      # write(secuencias[line + 1], file = "available_sequence_geneS_filtrado.fasta", append = TRUE)
    }
  }
}
print("Done")





