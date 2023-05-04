if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("polyester")
library("polyester")
library("Biostrings")

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation)")

tx_label <-  read.delim("004.final.tx.gene.with.label.srt.txt", 
                        header = FALSE, as.is =TRUE)
colnames(tx_label) <- c("tx_id", "gene_id", "type")
tx_label$proportion <- 1

single <- tx_label[tx_label$type == "single", ]

long <- tx_label[tx_label$type == "long", ]
long <- long[order(long$gene_id), ]
short <- tx_label[tx_label$type == "short", ]
short <- short[order(short$gene_id), ]
stopifnot(nrow(long) == nrow(short))
stopifnot(all(long$gene_id ==short$gene_id))



tx <- seq_gtf(
  gtf = "003.tx4simulation_sorted.gtf",
  seqs = "GRCh38.genome",
  feature = "transcript",
  exononly = TRUE,
  idfield = "transcript_id",
  attrsep = "; "
)

writeXStringSet(tx, '004.tx4simulation.fa')

# ~20x coverage ----> reads per transcript = transcriptlength/readlength * 20
# here all transcripts will have ~equal FPKM


tx_len_count  <- data.frame(tx_id = names(tx),
                      length = width(tx), 
                      readspertx = round(20 * width(tx) / 100))
## condition 1
## assign expression proportion of isoform transcripts
set.seed(12233)
long$proportion <- runif(nrow(long), min = 0, max = 1)
short$proportion <- 1 - long$proportion
tx_label_p <- rbind(single, long, short)

tx_len_count_p <- merge(tx_len_count, tx_label_p, by = "tx_id",
                        all.x = TRUE)
rownames(tx_len_count_p) <- tx_len_count_p[, 1]
tx_len_count_p$reads_per_isoform <- round(tx_len_count_p$readspertx * tx_len_count_p$proportion)

tx_len_count_p_1 <- tx_len_count_p[tx_len_count$tx_id, ]

## condition 2
set.seed(34567)
long1 <- long
short1 <- short
long$proportion <- runif(nrow(long), min = 0, max = 1)
short$proportion <- 1 - long$proportion

tx_label_p2 <- rbind(single, long, short)

tx_len_count_p2 <- merge(tx_len_count, tx_label_p2, by = "tx_id",
                        all.x = TRUE)
rownames(tx_len_count_p2) <- tx_len_count_p2[, 1]
tx_len_count_p2$reads_per_isoform <- round(tx_len_count_p2$readspertx * tx_len_count_p2$proportion)

tx_len_count_p_2 <- tx_len_count_p2[tx_len_count$tx_id, ]

stopifnot(rownames(tx_len_count_p_2) == rownames(tx_len_count_p_1))
stopifnot(rownames(tx_len_count_p_2) == tx_len_count$tx_id)

readmat <- data.frame(cond_1 = tx_len_count_p_1$reads_per_isoform,
                      cond_2 = tx_len_count_p_2$reads_per_isoform)
rownames(readmat) <- rownames(tx_len_count_p_1)
readmat <- as.matrix(readmat)

# simulation call:
simulate_experiment_countmat(fasta = "004.tx4simulation.fa", 
                             gtf = NULL, 
                             seqpath = NULL,
                             readmat = readmat, 
                             outdir = "simulated_reads", 
                             fraglen = 250, 
                             fragsd = 25, 
                             readlen = 100, 
                             error_rate = 0.005, 
                             paired = TRUE,
                             seed = 458930)

for (cov in c(40, 60, 80, 100, 120, 140)) {
    readmat_t <- round(readmat * cov/20)
    simulate_experiment_countmat(fasta = "004.tx4simulation.fa", 
                                 gtf = NULL, 
                                 seqpath = NULL,
                                 readmat = readmat_t, 
                                 outdir =paste0(cov, "Xsimulated_reads") , 
                                 fraglen = 250, 
                                 fragsd = 25, 
                                 readlen = 100, 
                                 error_rate = 0.005, 
                                 paired = TRUE,
                                 seed = 458930)
}




