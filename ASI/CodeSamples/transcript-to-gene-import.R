### LOAD REQUIRED LIBRARIES but first install if they are not installed

library(biomaRt)
library(tximport)
library(rhdf5)
library(Rsamtools)
library(readr)
library(tximportData)


setwd("/Users/jonah/School/ASI/final/kallisto")

#biomart tx2 object code
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
k <- keys(txdb, keytype = "TXNAME")
View(k)
tx2gene <- select(txdb, k, columns=c("GENEID", "TXNAME","TXID"),"TXNAME" )
tx2gene <- as.data.frame(tx2gene)
nonas <- tx2gene[!is.na(tx2gene$GENEID),]

### USE TXIMPORT TO SUMMARIZE TRANSCRIPT COUNTS INTO GENE COUNTS
## For multiple samples, each named as a folder in the kallisto directory (can be abundance.h5 or abundance.tsv file)
accessions <- list.dirs(full.names=TRUE)[-c(1:2)]
kallisto.dir<-paste0(accessions)
kallisto.files<-file.path("abundance.tsv") #can also be abundance.tsv
names(kallisto.files)<- accessions
jonah.files<-c("pramehi-ctrl1/abundance.tsv","pramehi-ctrl2/abundance.tsv","pramehi-ctrl3/abundance.tsv","pramehi-sirna1/abundance.tsv","pramehi-sirna2/abundance.tsv","pramehi-sirna3/abundance.tsv","pramelo-ctrl2/abundance.tsv","pramelo-ctrl3/abundance.tsv","pramelo-oe2/abundance.tsv","pramelo-oe3/abundance.tsv")
jonah.files1 <- c("pramehi-ctrl1/quant.sf","pramehi-ctrl2/quant.sf","pramehi-ctrl3/quant.sf","pramehi-sirna1/quant.sf","pramehi-sirna2/quant.sf","pramehi-sirna3/quant.sf","pramelo-ctrl2/quant.sf","pramelo-ctrl3/quant.sf","pramelo-oe2/quant.sf","pramelo-oe3/quant.sf")

jonah.files_big <- c("pramelo-ctrl1/abundance.tsv","pramelo-ctrl2/abundance.tsv","pramelo-ctrl3/abundance.tsv","pramelo-oe1/abundance.tsv","pramelo-oe2/abundance.tsv","pramelo-oe3/abundance.tsv","pramehi-ctrl1/abundance.tsv","pramehi-ctrl2/abundance.tsv","pramehi-ctrl3/abundance.tsv","pramehi-sirna1/abundance.tsv","pramehi-sirna2/abundance.tsv","pramehi-sirna3/abundance.tsv")
jonah.files1_big <- c("pramehi-ctrl1/quant.sf","pramehi-ctrl2/quant.sf","pramehi-ctrl3/quant.sf","pramehi-sirna1/quant.sf","pramehi-sirna2/quant.sf","pramehi-sirna3/quant.sf","pramelo-ctrl1/quant.sf","pramelo-ctrl2/quant.sf","pramelo-ctrl3/quant.sf","pramelo-oe1/quant.sf","pramelo-oe2/quant.sf","pramelo-oe3/quant.sf")
tx.kallistoOrdered <- tximport(jonah.files, type = "kallisto", tx2gene = nonas, countsFromAbundance ="no", ignoreAfterBar = TRUE)


### GENERATE TWO COLUMN OUTPUT 
counts<-as.data.frame(tx.kallistoOrdered$counts)
len <- as.data.frame(tx.kallisto$len[row.names(tx.kallisto$len) %in% genes, ])
ids<-rownames(counts)

### ROUND VALUES (DESEQ2 DOES NOT LIKE FRACTIONS), AND WRITE TO OUTPUT FILE
getwd()
write.table(round(counts), "outputsOrdered_kallisto_small.txt")

