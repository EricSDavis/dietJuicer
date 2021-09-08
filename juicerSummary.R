dir <- getwd()
files <- list.files(dir, "_inter_30.txt", recursive=TRUE)
#Loop for extracting data from each of the inter.txt files: 
final <- list()
list <- list()
for (i in 1:length(files)){
  x <- read.delim(files[i], header=F)
  list[[1]] <- as.numeric(gsub(",", "", unlist(strsplit(as.character(x[1,]), " "))[5]))
  list[[2]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[2,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[3]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[3,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[4]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[4,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[5]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[5,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[6]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[6,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[7]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[7,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[8]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[8,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[9]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[9,]), ": "))[2], "[(]"))[2], "%"))[1])
  list[[10]] <- as.numeric(gsub(",", "", unlist(strsplit(unlist(strsplit(as.character(x[10,]), ": "))[2], " "))[1]))
  list[[11]] <- as.numeric(gsub(",", "", unlist(strsplit(unlist(strsplit(as.character(x[11,]), ": "))[2], " "))[1]))
  list[[12]] <- as.numeric(gsub(",", "", unlist(strsplit(unlist(strsplit(as.character(x[12,]), ": "))[2], " "))[1]))
  list[[13]] <- as.numeric(gsub(",", "", unlist(strsplit(unlist(strsplit(as.character(x[13,]), ": "))[2], " "))[1]))
  list[[14]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[14,]), "[(]"))[2], "%"))[1])
  list[[15]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[14,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[16]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[15,]), "[(]"))[2], "%"))[1])
  list[[17]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[15,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[18]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[16,]), "[(]"))[2], "%"))[1])
  list[[19]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[16,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[20]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[17,]), "[(]"))[2], "%"))[1])
  list[[21]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[17,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[22]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[20,]), "[(]"))[2], "%"))[1])
  list[[23]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[20,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[24]] <- as.numeric(unlist(strsplit(unlist(strsplit(as.character(x[21,]), "[(]"))[2], "%"))[1])
  list[[25]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[21,]), "[(]"))[2], "%"))[2], "/ "))[2])
  list[[26]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[22,]), "[(]"))[3], "/")), "%"))[1])
  list[[27]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[22,]), "[(]"))[3], "/")), "%"))[3])
  list[[28]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[23,]), "[(]"))[3], "/")), "%"))[1])
  list[[29]] <- as.numeric(unlist(strsplit(unlist(strsplit(unlist(strsplit(as.character(x[23,]), "[(]"))[3], "/")), "%"))[3])
  data <- do.call(rbind, list)
  final[[i]] <- data
}

data <- data.frame(do.call(cbind, final))
colnames(data) <- unlist(strsplit(files, "/"))[seq(3,(3*length(files)),by=3)]

names <- c("Sequenced Read Pairs", "Normal Paired", "Chimeric Paired", 
           "Collisions", "Low MAPQ Collissions", "Unmapped", "MAPQ",
           "Ligation Motif Present", "Alignable", "Unique Reads",
           "PCR Duplicates", "Optical Duplicates", "Library Complexity Estimate", 
           "Intra-Fragment Reads 1", "Intra-Fragment Reads 2", "Below MAPQ Threshold 1", 
           "Below MAPQ Threshold 2", "Hi-C Contacts 1", "Hi-C Contacts 2", 
           "Ligation Motif Present 1", "Ligation Motif Present 2",
           "Inter-chromosomal 1", "Inter-chromosomal 2", "Intra-chromosomal 1", 
           "Intra-chromosomal 2", "Short Range 1", "Short Range 2", "Long Range 1", "Long Range 2")
rownames(data) <- names
write.table(data, paste0(dir, "/juicerSummary.txt"), sep="\t")