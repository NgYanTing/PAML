# PAML automation script
# Automate, as far as possible, from having a gene list to getting CodeML output
# Use PIP, RPUSD3, PTPRC and OPN5 as test list of genes
library(plyr)
library(dplyr)
library(biomaRt)

# Load libraries
marthuman = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", GRCh=37)
martchimp = useEnsembl(biomart="ensembl", dataset="ptroglodytes_gene_ensembl", GRCh=37)
martgor = useEnsembl(biomart="ensembl", dataset="ggorilla_gene_ensembl", GRCh=37)
martorang = useEnsembl(biomart="ensembl", dataset="pabelii_gene_ensembl", GRCh=37)
martmac = useEnsembl(biomart="ensembl", dataset="mmulatta_gene_ensembl", GRCh=37)

setwd("C:/Users/NGYT/PAMLauto")
# Read in gene list
Genelist <- read.table("autotestgenelist.txt", header=T)
Humanlist <- subset(Genelist, Genelist$Species=="Human")
Chimplist <- subset(Genelist, Genelist$Species=="Chimp")
Gorlist <- subset(Genelist, Genelist$Species=="Gorilla")
Oranglist <- subset(Genelist, Genelist$Species=="Orangutan")
Maclist <- subset(Genelist, Genelist$Species=="Macaque")

# Get sequences
# Human
Hgene <- getSequence(id=Humanlist$Ensembl_gene_ID, type="ensembl_gene_id", seqType="coding", mart=marthuman)
# Keep the longest sequence if repeated
Hgene <- cbind(nchar(Hgene$coding), Hgene)
colnames(Hgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Hgene <- Hgene[order(Hgene$Seq_length, decreasing=TRUE),]
Hgene <- subset(Hgene, duplicated(Hgene$Ensembl_gene_ID)==FALSE)
rownames(Hgene) <- 1:nrow(Hgene)

# Chimp
Cgene <- getSequence(id=Chimplist$Ensembl_gene_ID, type="ensembl_gene_id", seqType="coding", mart=martchimp)
Cgene <- cbind(nchar(Cgene$coding), Cgene)
colnames(Cgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Cgene <- Cgene[order(Cgene$Seq_length, decreasing=TRUE),]
Cgene <- subset(Cgene, duplicated(Cgene$Ensembl_gene_ID)==FALSE)
rownames(Cgene) <- 1:nrow(Cgene)

# Gorilla
Ggene <- getSequence(id=Gorlist$Ensembl_gene_ID, type="ensembl_gene_id", seqType="coding", mart=martgor)
Ggene <- cbind(nchar(Ggene$coding), Ggene)
colnames(Ggene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Ggene <- Ggene[order(Ggene$Seq_length, decreasing=TRUE),]
Ggene <- subset(Ggene, duplicated(Ggene$Ensembl_gene_ID)==FALSE)
rownames(Ggene) <- 1:nrow(Ggene)

# Orangutan
Ogene <- getSequence(id=Oranglist$Ensembl_gene_ID, type="ensembl_gene_id", seqType="coding", mart=martorang)
Ogene <- cbind(nchar(Ogene$coding), Ogene)
colnames(Ogene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Ogene <- Ogene[order(Ogene$Seq_length, decreasing=TRUE),]
Ogene <- subset(Ogene, duplicated(Ogene$Ensembl_gene_ID)==FALSE)
rownames(Ogene) <- 1:nrow(Ogene)

# Macaca
Mgene <- getSequence(id=Maclist$Ensembl_gene_ID, type="ensembl_gene_id", seqType="coding", mart=martmac)
Mgene <- cbind(nchar(Mgene$coding), Mgene)
colnames(Mgene) <- c("Seq_length", "Coding_seq", "Ensembl_gene_ID")
Mgene <- Mgene[order(Mgene$Seq_length, decreasing=TRUE),]
Mgene <- subset(Mgene, duplicated(Mgene$Ensembl_gene_ID)==FALSE)
rownames(Mgene) <- 1:nrow(Mgene)

#All sequences together
Hgenex <- merge(Hgene, Humanlist, by="Ensembl_gene_ID")
Hgenex <- Hgenex[order(Hgenex$Gene_symbol),]
Cgenex <- merge(Cgene, Chimplist, by="Ensembl_gene_ID")
Cgenex <- Cgenex[order(Cgenex$Gene_symbol),]
Ggenex <- merge(Ggene, Gorlist, by="Ensembl_gene_ID")
Ggenex <- Ggenex[order(Ggenex$Gene_symbol),]
Ogenex <- merge(Ogene, Oranglist, by="Ensembl_gene_ID")
Ogenex <- Ogenex[order(Ogenex$Gene_symbol),]
Mgenex <- merge(Mgene, Maclist, by="Ensembl_gene_ID")
Mgenex <- Mgenex[order(Mgenex$Gene_symbol),]
onelist <- rbind(Hgenex, Cgenex, Ggenex, Ogenex, Mgenex)
rownames(onelist) <- 1:nrow(onelist)

# Assemble sequences into fasta format for alignment with msa package
setwd("C:/Users/NGYT/PAMLauto/fasta_holder")
for(i in 1:4){ # 4 genes
  fasfile <-  character(length = 10) # 10 lines per gene into separate files
  fasfile[c(T,F,F,F,F,F,F,F,F,F)] = sprintf("> %s_H", onelist$Gene_symbol[i])
  fasfile[c(F,T,F,F,F,F,F,F,F,F)] = onelist$Coding_seq[i]
  fasfile[c(F,F,T,F,F,F,F,F,F,F)] = sprintf("> %s_C", onelist$Gene_symbol[i+4])
  fasfile[c(F,F,F,T,F,F,F,F,F,F)] = onelist$Coding_seq[i+4]
  fasfile[c(F,F,F,F,T,F,F,F,F,F)] = sprintf("> %s_G", onelist$Gene_symbol[i+8])
  fasfile[c(F,F,F,F,F,T,F,F,F,F)] = onelist$Coding_seq[i+8]
  fasfile[c(F,F,F,F,F,F,T,F,F,F)] = sprintf("> %s_O", onelist$Gene_symbol[i+12])
  fasfile[c(F,F,F,F,F,F,F,T,F,F)] = onelist$Coding_seq[i+12]
  fasfile[c(F,F,F,F,F,F,F,F,T,F)] = sprintf("> %s_M", onelist$Gene_symbol[i+16])
  fasfile[c(F,F,F,F,F,F,F,F,F,T)] = onelist$Coding_seq[i+16]
  writeLines(fasfile, paste(onelist$Gene_symbol[i],".fasta",sep=""))
}

# Running msa package
library(msa)
files <- list.files(pattern="*.fasta", full.names=T, recursive=F)
# Assemble all phylip files into one for input into PhyML
for(i in files){
  SeqFile <- readDNAStringSet(i)
  msaout <- msa(SeqFile, method="Muscle", order="input")
  i <- sub("./", "", i)
  i <- sub(".fasta", "", i)
  write.phylip(msaout, file=paste("C:/Users/NGYT/PAMLauto/phylip_out/",i,".phylip",sep=""))
}

# Input for PhyML
setwd("C:/Users/NGYT/PAMLauto/phylip_out")
phyfiles <- list.files(pattern="*.phylip", full.names=T, recursive=F)
file.create("allalignedforphyml.phylip")
for(i in phyfiles){
  phy <- readLines(i)
  cat(phy, file="allalignedforphyml.phylip", append=TRUE, sep = "\n")
  cat("", file="allalignedforphyml.phylip", append=TRUE, sep = "\n")
}

# Seqfiles for CodeML
# Write out all sequences into individual phylip files for input into CodeML
# CodeML only accepts sequences whose total number of nucleotides is divisible by 3.
# Add "-" to any sequences whose total number of nucleotides is not divisible by 3 while
# maintaining interleaved format.
for(i in phyfiles){
  phy <- readLines(i)
  phy[1] <-  paste(phy[1], "I", sep=" ") # Adding this "I" tells CodeML the seq format is interleaved
  ntnum <- as.integer(unlist(strsplit(phy[1], " "))[3])
  if(ntnum %% 3 == 0){
    filename <- gsub("phylip","",i)
    filename <- gsub("[.///]","",filename)
    writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                          filename, ".phylip", sep=""))
  } else if ((ntnum+1)%%3 == 0){
    if(ntnum%%10 != 0){
      ntnum <- ntnum+1
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "-", sep="") # We have 5 species, so input "-" to last 5 lines
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "-", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "-", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "-", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "-", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if (ntnum%%10 == 0) {
      ntnum <- ntnum+1
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], " -", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], " -", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], " -", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], " -", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], " -", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    }
    
  } else if ((ntnum+2)%%3 == 0){
    if(ntnum%%10 == 0){
      ntnum <- ntnum+2
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], " --", sep="") # We have 5 species, so input "-" to last 5 lines
      phy[length(phy)-4] <- paste(phy[length(phy)-4], " --", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], " --", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], " --", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], " --", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if(ntnum%%10 != 0 & (((ntnum+2) %% 10) == 1)) {
      ntnum <- ntnum+2
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "- -", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "- -", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "- -", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "- -", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "- -", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if(ntnum%%10 != 0 & (((ntnum+2) %% 10) <=9)) {
      ntnum <- ntnum+2
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "--", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "--", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "--", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "--", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "--", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    }
  } else if ((ntnum+3)%%3 == 0){
    if(ntnum%%10 != 0 & (((ntnum+3) %% 10) <= 7)){
      ntnum <- ntnum+3
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "---", sep="") # We have 5 species, so input "-" to last 5 lines
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "---", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "---", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "---", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "---", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if(ntnum%%10 == 0 ) {
      ntnum <- ntnum+3
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], " ---", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], " ---", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], " ---", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], " ---", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], " ---", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if(ntnum%%10 != 0 & (((ntnum+3) %% 10) == 1)) {
      ntnum <- ntnum+3
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "-- -", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "-- -", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "-- -", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "-- -", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "-- -", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    } else if(ntnum%%10 != 0 & (((ntnum+3) %% 10) == 2)) {
      ntnum <- ntnum+3
      phy[1] <- sub((unlist(strsplit(phy[1], " "))[3]),as.character(ntnum),phy[1])
      phy[length(phy)] <- paste(phy[length(phy)], "- --", sep="")
      phy[length(phy)-4] <- paste(phy[length(phy)-4], "- --", sep="")
      phy[length(phy)-3] <- paste(phy[length(phy)-3], "- --", sep="")
      phy[length(phy)-2] <- paste(phy[length(phy)-2], "- --", sep="")
      phy[length(phy)-1] <- paste(phy[length(phy)-1], "- --", sep="")
      filename <- gsub("phylip","",i)
      filename <- gsub("[.///]","",filename)
      writeLines(phy, paste("C:/Users/NGYT/PAMLauto/seqfiles/",
                            filename, ".phylip", sep=""))
    }
  }
}

# Run script until here first! ############################################################################
# PhyML requires a lot of settings after the program is run.
# Therefore not run from within R.
# allalignedforphyml.phylip was put through PhyML (command line) to get a file that contains all the trees.
# Read tree file in
setwd("C:/Users/NGYT/PAMLauto/phylip_out")
trees <- readLines("allalignedforphyml.phylip_phyml_tree.txt")
# Delete all the numbers that specify branch lengths
# This allows CodeML to generate dN/dS values for each branch
setwd("C:/Users/NGYT/PAMLauto/trees")
trees <- gsub(":\\d.\\d+,", ":,", trees)
trees <- gsub(":\\d.\\d+)", ":)", trees)
trees <- gsub("\\d.\\d+:", ":", trees)
trees <- gsub(":", "", trees)
genenames <- vector(mode="character",length=length(trees))
for(i in 1:length(trees)){
  treesnames <- gsub("\\(", "", trees[i])
  treesnames <- unlist(strsplit(treesnames, split='_', fixed=TRUE))[1]
  writeLines(trees[i], paste(treesnames,".tree", sep=""))
  genenames[i] <- treesnames
}

# Run CodeML
# codeml.exe and the required codeml.ctl are together in a folder called CodeML
# Make a new ctl file for every run
# Pre-set the other settings in .ctl file first (such as seqtype, model, NSsites etc.)
# "rainbow" was used as a placeholder for name substitutions
for(i in 1:length(genenames)){
  ctlfile <- readLines("C:/Users/NGYT/PAMLauto/Codeml/rainbowcodeml.ctl")
  ctlfile[1] <- sub("rainbow", genenames[i],ctlfile[1])
  ctlfile[2] <- sub("rainbow", genenames[i],ctlfile[2])
  ctlfile[3] <- sub("rainbow", genenames[i],ctlfile[3])
  writeLines(ctlfile, "C:/Users/NGYT/PAMLauto/Codeml/codeml.ctl") # Make sure this is the same folder as codeml.exe
  shell("C:/Users/NGYT/PAMLauto/Run_CodeML.bat", intern=T)
}
# Outputs are the "mlc" files from CodeML that contain dN/dS values and other information

# Session Info
# R version 3.3.1 (2016-06-21)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 7 x64 (build 7601) Service Pack 1

# locale:
#   [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#   [4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
#   [1] msa_1.4.5           Biostrings_2.40.2   XVector_0.12.1      IRanges_2.6.1       S4Vectors_0.10.3    BiocGenerics_0.18.0
#   [7] biomaRt_2.28.0      dplyr_0.5.0         plyr_1.8.4         

# loaded via a namespace (and not attached):
#   [1] Rcpp_0.12.7          XML_3.98-1.4         assertthat_0.1       bitops_1.0-6         R6_2.2.0             DBI_0.5-1           
#   [7] magrittr_1.5         RSQLite_1.0.0        zlibbioc_1.18.0      tools_3.3.1          Biobase_2.32.0       RCurl_1.95-4.8      
#   [13] AnnotationDbi_1.34.4 tibble_1.2   