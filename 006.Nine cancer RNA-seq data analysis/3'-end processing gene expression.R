
library(tximport)


files <- dir(".", "quant.sf$", recursive = TRUE)
names(files) <- gsub("/quant.sf", "", files)

## gene-level count
tx2gene <- read.delim("~/marcus.ruscetti-umw/haibo/docs/gencode.v34.primary_assembly.annotation.gene_tx.map.txt",
                      header = FALSE)
tx2gene <- tx2gene[, c(2,1)]
txi <- tximport(files, type = "salmon", txOut = FALSE,
                tx2gene = tx2gene)
txi <- txi$counts
txi <- txi[rowSums(txi) > 0, ]
write.table(txi, file = "Kevin.raw.gene.count.table.txt",
            sep = "\t", quote = FALSE, row.names =TRUE)
			
txi <- tximport(files, type = "salmon", txOut = FALSE,
                tx2gene = tx2gene)			
txi <- do.call(cbind, txi[1:3])

colnames(txi) <- paste(colnames(txi), rep(c("TPM", "count", "length"), each =6), sep = "_")

write.table(txi, file = "Kevin.raw.TPM.count.length.txt",
            sep = "\t", quote = FALSE, row.names =TRUE)


## DEseq2 analysis of sample distances
library("ggplot2")
library("sva")
library("DESeq2")
library("edgeR")
library("limma")
library("pheatmap")
library("RColorBrewer")
library("dplyr")
library("ggfortify")
library("genefilter")
library("gplots")
library("doParallel")
library("WriteXLS")


intergenic_read <- c(1064257, 11332057, 13100186,80484,1650387,1493514)
intergenic_length <- 1270113221
dna_bp_cov <- intergenic_read/intergenic_length


cur_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\05.featurecount.table\final.count.table)"
setwd(cur_dir)

filename <- "Nine.cancer.type.gene.count.table.txt"

count_PE <- read.delim(filename, sep = "\t",
                       header = TRUE, as.is =TRUE, check.names = FALSE)
colnames(count_PE)
metadata <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\2. data for case studies\2. data for case studies\5.Nine cancer type-paired data/50-paires full SraRunInfo.txt)"
metadata <- read.delim(metadata)[, 1:6]
metadata$PHENOTYPE <- factor(metadata$PHENOTYPE)
metadata$tissue <- factor(metadata$tissue)
metadata$source_name <- gsub(" ", ".", metadata$source_name, perl = TRUE)
metadata$source_name <- factor(metadata$source_name)
count_PE <- count_PE[rowSums(count_PE) != 0, ]

## total reads assigned to gene features
total <- data.frame(total = colSums(count_PE))
total <- cbind(metadata, total)


pdf("Fig. 1. Total number of reads assigned to gene features.pdf",
    width = 16, height =5)
ggplot(total, aes(SampleLabel, total, fill = tissue))+
  geom_col() + 
  theme(axis.text.x = 
          element_text(angle=90,
                       hjust = 0, vjust = 0.5))
dev.off()


## remove lowly expressed genes
cpms = cpm(count_PE)
keep = rowSums(cpms > 1) >= 4
count_PE = count_PE[keep, ]
rm("cpms", "keep")

dim(count_PE)  ## 22525
count_PE <- count_PE[, metadata$Run]
all(colnames(count_PE) == metadata$Run)

dds <- DESeqDataSetFromMatrix(countData = round(as.matrix(count_PE)),
                              colData = metadata,
                              design = ~0 + PHENOTYPE + source_name)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType="parametric", 
                           maxit=1000)

### exploratory analysis
vsd <- vst(dds, blind = TRUE)
sampleDists <- dist(t(assay(vsd)))

#### Heatmap showing sample distances
#### samples SP4 and SP8 are very different from SP5-7 of the same treatment
distancePlot <- function(sampleDists, sampleNames)
{
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- sampleNames
  colnames(sampleDistMatrix) <- sampleNames
  colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  pheatmap(sampleDistMatrix,
           clustering_distance_rows = sampleDists,
           clustering_distance_cols = sampleDists,
           col = colors)
  sampleDistMatrix 
}

## one sample normal:SRR1313203.1 is an outlier, so it should be 
## removed and reanalyze it.
pdf("Fig 2. Heatmap showing sample distances.pdf", 
    width = 20, height = 20)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists,
                                 sampleNames = vsd$SampleLabel)
dev.off()


## PCA using adjusted value
pc <- prcomp(t(assay(vsd)), scale = T, 
             center = T, retx = TRUE)
pdf(file = "Fig 3. PCA.plot.pdf", 
    width = 22.0, height = 20)
autoplot(pc, data = metadata, 
         colour = 'tissue', size = 4)
dev.off()



library(corrplot)
metadata <- metadata[order(metadata$SampleLabel), ]
types_cancer <- gsub("_.+", "", metadata$SampleLabel, perl = TRUE)
cancer_type <- as.data.frame(table(types_cancer))

vsd_expr <- assay(vsd)[, metadata$Run]
all(colnames(vsd_expr) == metadata$Run)
colnames(vsd_expr) <- paste(metadata$SampleLabel, colnames(vsd_expr), sep = "_")
pdf(file = "5.sample.correlation.pdf", height =20, width = 20)
corrplot.mixed(cor(vsd_expr),
               tl.pos = 'lt',
               diag = "n",
               lower = "number", 
               upper = "circle",
               tl.col = "black")
dev.off()

library(corrplot)
library(RColorBrewer)
M <-cor(vsd_expr)
pdf(file = "6.sample.correlation.pdf", height =20, width = 20)
corrplot(M, method = "ellipse",
         type="upper", order="hclust",
         hclust.method = "ward.D2",
         tl.col = "black",
         col=brewer.pal(n=8, name="RdYlBu"))
dev.off()

## sva analysis
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1. Normalized.count.table.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

## using sva to reduce hidden variations
get.sva <- function (expr.data=NULL, meta.data=NULL)
{
  mod <- model.matrix(~ PHENOTYPE + source_name, data=meta.data)
  mod0 <- model.matrix(~ 1, data=meta.data)
  
  num.sva <-svaseq(as.matrix(expr.data), mod, mod0)$n.sv
  sv <- svaseq(as.matrix(expr.data), mod, mod0,n.sv=num.sva)$sv
  
  colnames(sv)<- paste0("sv",1:num.sva)
  
  meta.data.sva <-cbind(meta.data,sv ) 
  
  meta.data.sva
}

## adjust expression for hidden variations for EDA plots
get.adj <- function(expr.data=NULL, design=NULL,  meta.data=NULL)
{
  ## adjusted expression using voom()
  v <- voom(expr.data, design=design)
  fit <- lmFit(v, design=design)
  adj.exp <- v$E -fit$coefficients[, 11:29] %*% t(design[, 11:29])
  adj.exp
}

## get sva: 4 hidden variables
meta.sva <- get.sva(expr.data= normCounts, meta.data=metadata)
meta.sva$group <- factor(paste(meta.sva$PHENOTYPE, 
                               gsub("_.+$", "", meta.sva$SampleLabel, perl = TRUE),
                               sep="_"))
model <- model.matrix(~ 0 + group + sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 +sv8 +sv9+sv10+sv11 +sv12 +sv13 +sv14 +sv15 +sv16 +sv17 +sv18 +sv19,
                      data=meta.sva)


## adjust for sva
adj.exp <- get.adj(expr.data=normCounts, design=model,  meta.data=meta.sva)
write.table(adj.exp, "Table 2. sva-adjusted.gene.expression.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

#### Heatmap showing sample distances after adjusting for hidden variation
sampleDists <- dist(t(adj.exp))
pdf("Fig 2.2 Heatmap showing SVA-adjusted sample distances.pdf", 
    width = 20, height = 20)
sampleDistMatrix <- distancePlot(sampleDists = sampleDists, 
                                 sampleNames = paste0(vsd$SampleLabel))
dev.off()

## PCA using adjusted value
pc2 <- prcomp(t(adj.exp), scale = T, 
              center = T, retx = TRUE)
pdf(file = "Fig 3.2. PCA.plot.SVA.adjusted.samples.pdf", 
    width = 22.0, height = 20)
autoplot(pc2, data = metadata, 
         colour = c("source_name"), size = 4)
dev.off()

## SVA is much better
dds <- DESeqDataSetFromMatrix(countData = as.matrix(round(count_PE)),
                              colData = meta.sva,
                              design = ~0 + group + sv1 + sv2 + sv3 + sv4 + sv5 + sv6 + sv7 +sv8 +sv9+sv10+sv11 +sv12 +sv13 +sv14 +sv15 +sv16 +sv17 +sv18 +sv19)

dds <- estimateSizeFactors(dds)

dds <- estimateDispersions(dds, fitType="parametric", 
                           maxit=1000)
normCounts <- round(counts(dds, normalized = TRUE, replaced = FALSE))
write.table(normCounts, "Table 1.1 Normalized.count.table.SVA-adjusted.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

contrast.matrix <- matrix(c(1, rep(0, 8), -1, rep(0, 8), rep(0, 19),
                            rep(0, 1), 1, rep(0, 8), -1, rep(0,7), rep(0, 19),
                            rep(0, 2), 1, rep(0, 8), -1, rep(0,6), rep(0, 19),
                            rep(0, 3), 1, rep(0, 8), -1, rep(0,5), rep(0, 19),
                            rep(0, 4), 1, rep(0, 8), -1, rep(0,4), rep(0, 19),
                            rep(0, 5), 1, rep(0, 8), -1, rep(0,3), rep(0, 19),
                            rep(0, 6), 1, rep(0, 8), -1, rep(0,2), rep(0, 19),
                            rep(0, 7), 1, rep(0, 8), -1, rep(0,1), rep(0, 19),
                            rep(0, 8), 1, rep(0, 8), -1, rep(0, 19)
), nrow = 9, byrow = TRUE)
rownames(contrast.matrix) <- c("CSCC_TN", "ESCC_TN", "GAC_TN", 
                               "HCC_TN", "LUAD_TN", "LUSC_TN",
                               "PTC_TN", "SCLC_TN", "SRCC_TN")
memory.limit(size=56000)
dds <- nbinomWaldTest(dds, modelMatrix = NULL,
                      betaPrior=FALSE,
                      maxit = 50000, 
                      useOptim = TRUE, 
                      quiet = FALSE, 
                      useT = FALSE, useQR = TRUE)

gene_name <- read.delim("00.gencode.v34.geneid.name.mapping.txt", 
                        header = FALSE, as.is = TRUE)
colnames(gene_name) <- c("Gene_id", "Symbol")

output_DESeq <- function(i, dds, contrast.matrix, 
                         threshold, shrink)
{
  # i=1
  # threshold =1
  # shrink = FALSE
  res <- results(dds, alpha = 0.01,
                 contrast=contrast.matrix[i,], 
                 lfcThreshold=log2(threshold),
                 format = "DataFrame",
                 altHypothesis="greaterAbs") 
  if (shrink)
  {
    res <- lfcShrink(dds,
                     contrast = contrast.matrix[i,],
                     res = res,
                     format = "DataFrame",
                     lfcThreshold = log2(threshold),
                     type = "ashr")
  }
  res <- data.frame(res)
  res <- merge(res, gene_name, by= "Gene_id", all.x = TRUE)
  res 
}

## apply raw FLC
DESeq_out <- lapply(1:nrow(contrast.matrix), 
                    output_DESeq, dds = dds, 
                    contrast.matrix = contrast.matrix, 
                    threshold = 1, shrink = FALSE)

WriteXLS(x = DESeq_out, 
         ExcelFileName = "Table 3.2. DESeq.differential expressed genes.xlsx",
         row.names = FALSE, 
         SheetNames = rownames(contrast.matrix))


## 3' end processing genes
adj.exp <- read.delim("Table 2. sva-adjusted.gene.expression.txt")
rownames(adj.exp) <- gsub("\\.\\d+", "", rownames(adj.exp), perl = TRUE)
apa_genes <- read.delim(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS/GO0031124-mRNA 3' end processing genes.mart_export.txt)")[, 1:3]
apa_genes$Gene.stable.ID.version <- gsub("\\.\\d+", "", apa_genes$Gene.stable.ID.version, perl = TRUE)

apa_gene_expr <- adj.exp[rownames(adj.exp) %in% apa_genes$Gene.stable.ID.version, ]

apa_gene_expr <- apa_gene_expr[, sort(colnames(apa_gene_expr))]

expressed_apa_genes <- unique(apa_genes[apa_genes$Gene.stable.ID.version %in% rownames(apa_gene_expr), c(2,3) ])
exprsd_apa_genes <- expressed_apa_genes[, 2] 
names(exprsd_apa_genes) <- expressed_apa_genes[, 1]
rownames(apa_gene_expr) <- exprsd_apa_genes[rownames(apa_gene_expr)] 


metadata <- metadata[order(metadata$Run), ]
all(metadata$Run == colnames(apa_gene_expr))
colnames(apa_gene_expr) <- metadata$SampleLabel
metadata$source_name <- gsub("_.+$", "", metadata$SampleLabel,
                             perl = TRUE)
metadata <- metadata[order(metadata$source_name, metadata$SampleLabel), ]
apa_gene_expr <- apa_gene_expr[, metadata$SampleLabel]
all(metadata$SampleLabel == colnames(apa_gene_expr))

apa_gene_expr_fc <- apa_gene_expr[, seq(2, 100, by =2)] - apa_gene_expr[, seq(1, 100, by =2)]

apa_gene_expr_fc <- apa_gene_expr_fc[order(rownames(apa_gene_expr_fc)), ]

## order colnames to matching oncoprint orders
load("APA.events.across.9.cancer.types.RData")
colnames(apa_gene_expr_fc) <- gsub("_T$", "", 
                                   colnames(apa_gene_expr_fc), perl = TRUE)
colnames(dPDUI_cat)[34] <- "CSCC_T0010"
apa_gene_expr_fc <- apa_gene_expr_fc[, colnames(dPDUI_cat)]

annotation_col <- data.frame(Source = factor(gsub("_.+$", "", colnames(apa_gene_expr_fc),
                                                  perl = TRUE)))
rownames(annotation_col) <- colnames(apa_gene_expr_fc)

apa_gene_expr_fc_threshold <- as.matrix(apa_gene_expr_fc)
#apa_gene_expr_fc_threshold[abs(apa_gene_expr_fc_threshold)< log2(1.5)] <- 0
#apa_gene_expr_fc_threshold <- 
#      apa_gene_expr_fc_threshold[rowSums(apa_gene_expr_fc_threshold) >0, ]

apa_gene_expr_fc_threshold  <- data.frame(t(apa_gene_expr_fc_threshold))

apa_gene_expr_fc_threshold_split <- split(apa_gene_expr_fc_threshold, 
                                          f = annotation_col$Source)
apa_gene_expr_fc_threshold_split <- lapply(apa_gene_expr_fc_threshold_split,
                                           function(x) {
                                             colMeans(x)
                                           })
apa_gene_expr_fc_threshold_split <- do.call("cbind", apa_gene_expr_fc_threshold_split)

apa_gene_expr_fc_threshold_split <- apa_gene_expr_fc_threshold_split[, unique(gsub("_.+$", "", colnames(dPDUI_cat)))]
apa_gene_expr_fc_threshold_split <- apa_gene_expr_fc_threshold_split[sort(rownames(apa_gene_expr_fc_threshold_split)), ]
annotation_col_2 <- unique(annotation_col)
rownames(annotation_col_2) <-  gsub("_.+$", "", rownames(annotation_col_2),
                                    perl = TRUE)
cols <- rev(rainbow(9))
names(cols) <- annotation_col_2$Source
annotation_colors <- list(Source= cols)


genes <- rownames(apa_gene_expr_fc_threshold_split)
CPSF <- c(genes[10:13], "FIP1L1", "WDR33")
CSTF <- c(genes[16:19])
CFIm <- c(genes[14:15],"NUDT21")
CFIIm <- c("PCF11", "CLP1")
polyABP <- c("SYMPK", "PAPOLG", "PABPC1", 
             "PABPC4", "PABPN1", "RBBP6",
             "PPP1CA", "PPP1CB")
apa <-c(CPSF,CSTF,CFIm,CFIIm,polyABP)
other <- genes[!genes %in% apa]

apa_gene_expr_fc_threshold_split <- 
  apa_gene_expr_fc_threshold_split[c(apa, other), ]
library("pheatmap")
pdf("Fig 6.2 Heatmap showing expression of genes involved in mRNA 3' end processing.pdf",
    height = 5 , width = 2.8)
heat <- pheatmap(apa_gene_expr_fc_threshold_split,
                 scale = "none",
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 clustering_method ="ward.D2",
                 clustering_distance_rows = "correlation",
                 border_color= NA,
                 treeheight_row = 20,
                 treeheight_col = 10,
                 annotation_col = annotation_col_2,
                 annotation_colors =annotation_colors,
                 show_colnames= FALSE,
                 show_rownames= TRUE,
                 fontsize_row = 5)
print(heat)
dev.off()