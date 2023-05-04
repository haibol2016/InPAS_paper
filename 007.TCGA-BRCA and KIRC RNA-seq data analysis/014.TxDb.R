#!/bin/env Rscript

args <- commandArgs(trailingOnly=TRUE)

library("GenomicFeatures")

# "~/work/mccb/genome/Homo_sapiens/human_gencode_v34/GRCh38.primary_assembly.genome.chr.length"
chrom_info <- read.delim(args[1],
                         header = FALSE, as.is = TRUE)
colnames(chrom_info) <- c("chrom", "length")
#chrom_info$isCircular  <- ifelse(chrom_info$seqnames == "chrM", TRUE, FALSE)

txdb <- makeTxDbFromGFF(file = args[2], 
                format="gtf", 
                dataSource= "Ensembl", 
                organism="Homo sapiens" , 
                taxonomyId=NA, 
                circ_seqs="chrM", 
                chrominfo=chrom_info, 
                miRBaseBuild=NA)

saveDb(txdb, file=args[3])

