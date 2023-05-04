## combine salmon output to get transcript-level expression count
library(tximport)
setwd("results/023.salmon.quant.out/002.HeLa.RNA-seq/")

files <- dir(".", "quant.sf$", recursive = TRUE)
names(files) <- gsub("/quant.sf", "", files)

txi <- tximport(files, type = "salmon", txOut = TRUE)

counts <- txi$counts
len <- txi$length
count_per_base <- as.data.frame((counts *100)/ len)
count_per_base$total_counts <- rowSums(counts)
count_per_base$avg_length <- rowMeans(len)
  
write.table(count_per_base, file = "HeLa.averaged.base-wise.coverage.txt",
            sep = "\t", quote = FALSE, row.names = TRUE)
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


keep <- rowSums(dat[,1:4]) > 0
dat <- dat[keep, ]
dat$total_basecov <-  log10(rowSums(dat[,1:4]))

## expression density per transcript biotype
library("ggplot2")
library("ggridges")
pdf("Fig 1.1.Expression.level.by.tx.biotypes.pdf", width = 10, height = 8)
ggplot(dat, aes(x = total_basecov, y = tx_type, fill = tx_type)) +
  geom_density_ridges() +
  theme_ridges() +
  xlab("log10(Base Cov)")+
  theme(legend.position = "none")

ggplot(dat, aes(y = total_basecov, x = tx_type, alpha = 0.5)) +
  geom_jitter(position = position_jitter(0.2), color = "black") + coord_flip()
dev.off()
  
ncov <- do.call("c", lapply(1:100, function(x){
  sum(10^dat$total_basecov <= x)
}))

cov <- data.frame(cov = 1:100, n= ncov)
pdf("Fig 2.coverage.above.thersholds.pdf", width = 8, height = 4)
ggplot(cov, aes(x = cov, y = ncov))+
  geom_point()+
  geom_line()
dev.off()

## per-base coverage < 10 and non-coding tx of protein-coding genes
## 88,208
low <- dat[10^dat$total_basecov < 10 |
             (dat$gene_type == "protein_coding" &
                dat$tx_type != "protein_coding") , ] 

high <- dat[!rownames(dat) %in% rownames(low), ]
high_long <- high[high$tx_type == "protein_coding" | high$avg_length > 200, ] 
high_long <- high_long[!high_long$tx_type %in% c("snoRNA",
                                                "miRNA", "Mt_rRNA",
                                                "Mt_tRNA", "rRNA", 
                                                "rRNA_pseudogene",
                                                "scaRNA", "scRNA",
                                                "snRNA","vaultRNA"), ]

pdf("Fig 1.1.Expression.level.of.kept.high.long.tx.by.tx.biotypes.pdf", width = 10, height = 6)
ggplot(high_long, aes(x = total_basecov, y = tx_type, fill = tx_type)) +
  geom_density_ridges() +
  theme_ridges() +
  xlab("log10(Base Cov)")+
  theme(legend.position = "none")

ggplot(high_long, aes(y = total_basecov, x = tx_type, alpha = 0.5)) +
  geom_jitter(position = position_jitter(0.2), color = "black") + coord_flip()
dev.off()

save(dat, high_long, file = "Filtering.tx.by.abundance.length.RData")
write.table(high_long, file = "00.final.high.long.Tx.txt",
            row.names = TRUE, quote = FALSE, sep = "\t")