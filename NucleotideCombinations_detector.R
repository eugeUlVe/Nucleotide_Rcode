# **************************************************************************
# *
# * Author:     
# *
# * Center for Biological Research Margarita Salas (CIB), CSIC 
# * Biocomputing Unit (BCU) of the National Center for Biotechnology (CNB), CSIC
# * Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
# *
# Short description: (FIRST DRAFT) This program selects nucleotides combinations from sequences aligned in a fasta file for statistical analysis using a dataset. 
# 
# **************************************************************************

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

if(length(args) < 9){
  stop("USE: Rscript NucleotideCombinations_detector.R <DNA_sequences.fasta> <porcentageMINaccepted> <porcentageMAXaccepted> <startPos_CDS> <endPos_CDS> <outputName.fasta> <BioChemMolCli.tsv> <ColumnSelectedTSVfile> <NameOf_ColumnSelectedTSVfile>")
}

f_read <- scan(args[1], what = character(), quiet = TRUE)

# Build the loci table
head_sm <- grep(">", f_read) 
num_rows_lt <- length(head_sm)

# Vector with sequence positions 
sec_pos <- c()
for(sp in 2:length(head_sm)){
  sec <- head_sm[sp] - 1
  sec_pos <- c(sec_pos, sec)
}
sec_pos <- c(sec_pos, length(f_read)) # add last sequence
# print(sec_pos)

# Number of columns for loci table
secu_long <- f_read[sec_pos[1]]
spl_secu_long <- strsplit(secu_long, "")
num_cols_lt <- length(spl_secu_long[[1]])

loci_table <- matrix(1, nrow = num_rows_lt, ncol = num_cols_lt)
# loci_table <- matrix(1, nrow = length(head_sm), ncol = as.numeric(args[2]))

# Make the loci table
for(ele in 1:length(sec_pos)){
  secu <- f_read[sec_pos[ele]]
  spl_secu <- strsplit(secu, "")
  for(n in 1:length(spl_secu[[1]])){
    # print(spl_secu[[1]][n])
    loci_table[ele, n] <- spl_secu[[1]][n]
  }
}
# View(loci_table)
print("Done loci table")

# Search for different nucleotides for a DNA locus
n_muestras <- num_rows_lt
n_nt <- num_cols_lt
por_cutoff_min <- args[2]
por_cutoff_max <- args[3]
inicio_cds <- args[4]
fin_cds <- args[5]

vector_posSeleccionadas <- c()
for(j in 1:ncol(loci_table)){
  num_a <- 0
  num_t <- 0
  num_g <- 0
  num_c <- 0
  num_n <- 0
  max_v_nts <- 0
  geno <- 0
  #
  num_r <- 0
  num_y <- 0
  num_s <- 0
  num_w <- 0
  num_k <- 0
  num_m <- 0
  num_b <- 0
  num_d <- 0
  num_h <- 0
  num_v <- 0
  num_gap <- 0
  #
  for(i in 1:nrow(loci_table)){
    geno <- geno + 1
    if(loci_table[i,j] == "a"){
      num_a <- num_a + 1
    }
    if(loci_table[i,j] == "t"){
      num_t <- num_t + 1
    }
    if(loci_table[i,j] == "g"){
      num_g <- num_g + 1
    }
    if(loci_table[i,j] == "c"){
      num_c <- num_c + 1
    }
    if(loci_table[i,j] == "n"){
      num_n <- num_n + 1
    }
    #
    if(loci_table[i,j] == "r"){
      num_r <- num_r + 1
    }
    if(loci_table[i,j] == "y"){
      num_y <- num_y + 1
    }
    if(loci_table[i,j] == "s"){
      num_s <- num_s + 1
    }
    if(loci_table[i,j] == "w"){
      num_w <- num_w + 1
    }
    if(loci_table[i,j] == "k"){
      num_k <- num_k + 1
    }
    if(loci_table[i,j] == "m"){
      num_m <- num_m + 1
    }
    if(loci_table[i,j] == "b"){
      num_b <- num_b + 1
    }
    if(loci_table[i,j] == "d"){
      num_d <- num_d + 1
    }
    if(loci_table[i,j] == "h"){
      num_h <- num_h + 1
    }
    if(loci_table[i,j] == "v"){
      num_v <- num_v + 1
    }
    if(loci_table[i,j] == "-"){
      num_gap <- num_gap + 1
    }
    
    #
  }
  v_nts <- c(num_a, num_t, num_g, num_c, num_r, num_y, num_s, num_w, num_k, num_m, num_b, num_d, num_h, num_v)

  d_cero <- 0
  quitar_cero <- c()
  for(n_nt in 1: length(v_nts)){
    if(v_nts[n_nt] != 0){
      d_cero <- d_cero + 1
      quitar_cero <- c(quitar_cero, v_nts[n_nt])
    }
  }

  if(length(quitar_cero) > 1){
  max_quitar_cero <- as.numeric(max(quitar_cero))
  por_max_quitar_cero <- as.numeric((max_quitar_cero/as.numeric(n_muestras))*100)
  min_quitar_cero <- as.numeric(min(quitar_cero))
  por_min_quitar_cero <- as.numeric((min_quitar_cero/as.numeric(n_muestras))*100)
  if(d_cero >= 2 && por_min_quitar_cero >= as.numeric(args[2]) && por_max_quitar_cero >= as.numeric(args[3])){
    if(j >= as.numeric(args[4]) && j <= as.numeric(args[5])){ 
    info = paste("pos: ",j , "; A: ", num_a, "; T: ", num_t,  "; G: ", num_g, "; C: ", num_c, "; N: ", num_n, "; R: ", num_r, "; Y: ", num_y, "; S: ", num_s, "; W: ", num_w, "; K: ", num_k, "; M: ", num_m, "; B: ", num_b, "; D: ", num_d, "; H: ", num_h, "; V: ", num_v, "; gaps: ", num_gap, sep = "")
    # print(info)
    print(paste(por_min_quitar_cero, por_max_quitar_cero, sep = ";"))
    vector_posSeleccionadas <- c(vector_posSeleccionadas, j)
     }
    }
  }
}
# print(vector_posSeleccionadas)

# Write a file with positions where nucleotides are exchanges
etiqueta <- paste("Position where there is a nucleotide exchanged: ", length(vector_posSeleccionadas), sep = "")
print(etiqueta)
#
counter_posSeleccionadas <- 0
for(vps in 1:length(vector_posSeleccionadas)){
  counter_posSeleccionadas <- counter_posSeleccionadas + 1
  linea_pos_sel <- paste("Position ", counter_posSeleccionadas, ": ", vector_posSeleccionadas[vps])
  write(linea_pos_sel, file = "SelectedPosition_NucleotideCombinations.txt", append = TRUE)
}
#

eje_y <- rep(0, length(vector_posSeleccionadas))
jpeg('Rplots_Position_nucleotidesExchanged.jpg')
plot(vector_posSeleccionadas, eje_y, "p", xlab = paste("Positions where nucleotides are excahnges in ", args[9], sep = ""))
dev.off()
print("Plot position detected done")
# summary(warnings())

# Make a vector with sequences created pasting position selected above
vec_haplo <- c()
for(aim_r in 1:num_rows_lt){
  haplo <- loci_table[aim_r, vector_posSeleccionadas[1]]
  for(aim_c in 2:length(vector_posSeleccionadas)){
    haplo <- paste(haplo, loci_table[aim_r, vector_posSeleccionadas[aim_c]], sep = "")
  }
  vec_haplo <- c(vec_haplo, haplo)
}
# print(vec_haplo)
# length(vec_haplo)

# Remove duplicates and count how many they are duplicated. Then sequences are printed in a fasta file.
vec_haplo_noRepes <- unique(vec_haplo)
# print(vec_haplo_noRepes)
# length(vec_haplo_noRepes)

count_sequences <- 0
for(vhn in 1:length(vec_haplo_noRepes)){
  count_repes <- 0
  count_sequences <- count_sequences + 1
  for(vh in 1:length(vec_haplo)){
    if(vec_haplo_noRepes[vhn] == vec_haplo[vh]){
      count_repes <- count_repes + 1
    }
  }
  # print(paste(count_repes, vec_haplo_noRepes[vhn], sep = " "))
  fragment <- paste(">NucleotideCombination_", count_sequences, "_", "times_", count_repes, sep = "")
  # fragment <- paste(">Sequence_", count_sequences, "_", "times_", count_repes, sep = "")
  write(fragment, file = args[6], append = TRUE)
  # write(fragment, file = "SequenceAlignment.fasta", append = TRUE)
  write(vec_haplo_noRepes[vhn], file = args[6], append = TRUE)
  # write(vec_haplo_noRepes[vhn], file = "SequenceAlignment.fasta", append = TRUE)
}
print("File counting combination with sequences is done")

# Identified each nucleotide combination to a ID sample
print("Plotting")
datos <- read.csv(args[7], header = TRUE, sep = "\t")
# datos <- read.csv("seqcovid_20210804.tsv", header = TRUE, sep = "\t")

colum_data <- args[8]
# colum_data <- 11 # this is for Ct/Cq.PCR.gen.S
n_filas_data <- length(datos[,1])
# print(length(datos[,1]))

# print(length(vector_posSeleccionadas))
# print(nrow(loci_table))
sm_sec <- grep(">", f_read)
# print(length(sm_sec))

counterSeq <- 0
eje_x_fn <- c()
eje_x_f <- c()
eje_y <- c()
count_grupos_f <- 0
count_grupos <- 0
eje_x <- c()
vec_combi_posi<- c()
for(combination in 1:length(vec_haplo_noRepes)){
  counterSeq <- counterSeq + 1
  count_grupos_f <- count_grupos_f + 1
  count_grupos <- count_grupos + 1
  grupo <- vec_haplo_noRepes[combination]
  s_combi <- strsplit(vec_haplo_noRepes[combination], "")
  # print(grupo)
  # print(s_combi)
  eje_y_grupo <- c()
  vec_grupo <- c()
  vec_grupo_f <- c()
  for(sample in 1:num_rows_lt){
    count <- 0
    for(nt in 1:length(s_combi[[1]])){
      if(s_combi[[1]][nt] == loci_table[sample, vector_posSeleccionadas[nt]]){
        count <- count + 1
      }
    }
    # print(count)
    if(count == length(vector_posSeleccionadas)){
      vec_grupo <- c(vec_grupo, f_read[sm_sec [sample]])
      # Select sequence ID for selected column for analysis which are not empty 
      for(id in 1:n_filas_data){
        sm_dato <- paste(">", datos[id,1], sep = "")
        if(f_read[sm_sec [sample]] == sm_dato){
          # Remove NA and O and 99 values
          if(is.na(datos[id,as.numeric(args[8])]) == FALSE){
            if(datos[id,as.numeric(args[8])] != 0 && datos[id,as.numeric(args[8])] != 99){
            vec_grupo_f <- c(vec_grupo_f, f_read[sm_sec [sample]])
            count_grupos_f <- count_grupos_f + 1
            eje_y <- c(eje_y, datos[id,as.numeric(colum_data)])
            eje_y_grupo <- c(eje_y_grupo, datos[id,as.numeric(colum_data)])
            }
          }
        }
      }
      # print(vec_grupo_f)
      # print(length(vec_grupo_f))
    }
  }
  # print(vec_grupo_f)
  # print(length(vec_grupo_f))
  # print(eje_y_grupo)
  # print(length(eje_y_grupo))
 
  # Write file for statistical analysis
  grupo_numComponents <- paste(">", paste("NucleotideCombination_", counterSeq, sep = ""), grupo, "times:", length(vec_grupo_f), sep = " ")
  write(grupo_numComponents, file = "statistical_data.txt", append = TRUE)
  
  for(pares in 1:length(vec_grupo_f)){
    dos_values <- paste(vec_grupo_f[pares], eje_y_grupo[pares])
    write(dos_values, file = "statistical_data.txt", append = TRUE)
   }
  
  grupos <- paste("grupo", count_grupos, sep = "")
  grupo_f <- paste("grupo", count_grupos_f, sep = "")
  eje_x_fn <- c(eje_x_fn, rep(grupo, length(vec_grupo_f))) # 
}

# Make a data frame for plotting
dt <- data.frame(Grupo=eje_x_fn, CargaViral=eje_y)

library(ggplot2)

ggplot(data = dt, aes(x = eje_x_fn, y = eje_y)) +  geom_point(size=1) + xlab("Nucleotide Combinations") + ylab(args[9]) 
# ggplot(data = dt, aes(x = eje_x_fn, y = eje_y)) +  geom_point(size=1) + xlab("Nucleotide Combinations") + ylab("Ct/Cq PCR gen S") 
ggsave("Rplots_NucleotideCombinations.pdf")
summary(warnings())
unlink("Rplots.pdf")
print("Plot is done")

