library(tximport)

files <- dir(".", "quant.sf$", recursive = TRUE)
names(files) <- gsub("/quant.sf", "", files)

txi <- tximport(files, type = "salmon", txOut = TRUE)

counts <- txi$counts
len <- txi$length
count_per_base <- as.data.frame((counts *100)/ len)
count_per_base$total_counts <- rowSums(counts)
count_per_base$avg_length <- rowMeans(len)


biotype <- read.delim("/home/hl84w/arthur_mercurio/Haibo/InPAS/human_gencode_v34/Gencode.v34.transcript.biotype.MANE.txt",
                      header = FALSE, as.is = TRUE)
gene_tx_map <- read.delim("/home/hl84w/arthur_mercurio/Haibo/InPAS/human_gencode_v34/gencode.v34.primary_assembly.annotation.gene_tx.map.txt",
                          header = FALSE, as.is = TRUE)
colnames(gene_tx_map) <- c("gene_id", "tx_id", "symbol")
gene_tx_map_biotype <- merge(gene_tx_map,biotype, by.x = "tx_id",
                             by.y = "V1", all.y = TRUE)
colnames(gene_tx_map_biotype)[4:6] <- c("gene_type", "tx_type", "is_MANE")

dat <- merge(count_per_base, gene_tx_map_biotype, 
             by.x ="row.names",
             by.y = "tx_id", all.x = TRUE)
rownames(dat) <- dat[, 1]
dat <- dat[, -1]
write.table(dat, file = "Nine.cancer.averaged.base-wise.coverage.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)

## metadata
metadata <- read.delim("/home/hl84w/arthur_mercurio/Haibo/InPAS/docs_2/50-paires full SraRunInfo.txt",
                       as.is = TRUE, header = TRUE)
metadata <- metadata[, 1:2]
metadata$patient_id <- gsub("_[TN]$", "", metadata$SampleLabel, perl = TRUE)
dat <- dat[, c(metadata$Run, colnames(dat)[101:107])]

out_dir <- "Nine.cancer.abundant.long.tx"

if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

null <- lapply(seq(1,100, by =2), function(i){
  tmp <- dat[, c(i, i+1, 101:107)]
  keep <- rowSums(tmp[,1:2]) > 0
  tmp <- tmp[keep, ]
  tmp$total_basecov <-  log10(rowSums(tmp[,1:2]))
  
  ## per-base coverage < 10 and non-coding tx of protein-coding genes
  ## 128,534
  low <- tmp[10^tmp$total_basecov < 10 |
               (tmp$gene_type == "protein_coding" &
                  tmp$tx_type != "protein_coding"), ] 
  
  high <- tmp[!rownames(tmp) %in% rownames(low), ]
  high_long <- high[high$tx_type == "protein_coding" | high$avg_length > 200, ] 
  high_long <- high_long[!high_long$tx_type %in% c("snoRNA",
                                                   "miRNA", "Mt_rRNA",
                                                   "Mt_tRNA", "rRNA", 
                                                   "rRNA_pseudogene",
                                                   "scaRNA", "scRNA",
                                                   "snRNA","vaultRNA"), ]
  write.table(high_long, 
              file = file.path(out_dir, 
                               paste0(metadata$patient_id[i], ".averaged.base-wise.coverage.of.abundant.Tx.kept.txt")),
              sep = "\t", quote = FALSE, row.names = TRUE)
})








