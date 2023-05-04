
## APA genes compared to cancer genes downloaded from COSMIC
indir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types/)" 
cosmic <- read.delim(file.path(indir, "COSMIC cancer_gene_census.txt"), 
                     header = TRUE, as.is = TRUE)
cosmic$Ensembl_ID <- gsub(".+?(ENSG[^,]+).+", "\\1", cosmic$Synonyms, perl = TRUE)
cosmic$Ensembl_ID_s <- gsub("\\.\\d+", "", cosmic$Ensembl_ID, perl = TRUE)
apa <- read.delim(file.path(indir, "Significant.recurrent.APA.genes-annotated.txt"),
                  header = TRUE,  as.is = TRUE)
apa$Ensembl_ID_s <- gsub("\\.\\d+", "", apa$Gene_ID, perl = TRUE)

## merge
apa_cosmic <- merge(apa, cosmic, by = "Ensembl_ID_s", all.x = TRUE)

colnames(apa_cosmic)

apa_cosmic <- apa_cosmic[ , -c(3,4,7,29)]
colnames(apa_cosmic)[3:4] <- c("Gene_symbol", "Tx_biotype")

## add shortening/lengthening frequency
apa_freq <- read.delim(file.path(indir,"00.10percen.recurrent.APA.events.across.9.cancer.types.txt"),
                       header = TRUE, as.is = TRUE)

freq <- lapply(as.data.frame(t(apa_freq)), function(.x) {
  dat <- as.data.frame(table(.x))
  dat$.x <- as.character(dat$.x)
  dat$.x <- ifelse(dat$.x == "", "nsc", dat$.x)
  
  if (!"nsc" %in% dat$.x){dat <- rbind(dat, data.frame(.x = "nsc", Freq = 0 ))}
  else if (!"shortening" %in% dat$.x) {dat <- rbind(dat, data.frame(.x = "shortening", Freq = 0 ))}
  else if (!"lengthening" %in% dat$.x)  {dat <- rbind(dat, data.frame(.x = "lengthening", Freq = 0 ))}
  dat <- dat[order(dat$.x), ]
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1, drop = FALSE]
  dat <- t(dat)
  dat/50
  })
freq <- do.call(rbind, freq)
rownames(freq) <- rownames(apa_freq)
rownames(freq) <- gsub("ENSG\\d+\\.\\d+_", "", rownames(freq), perl = TRUE)

## merge with apa_cosmic

apa_cosmic <- merge(freq, apa_cosmic, by.x = "row.names", 
                    by.y = "Tx_ID", all.x = TRUE)

## pancancer.sig.mutated.genes

pancan <- read.delim(file.path(indir, "Pancancer.significantly.mutated.genes.txt"),
                     header = TRUE, as.is = TRUE)

pancan_genes <- unique(pancan[, 1])
apa_cosmic$Pan_cancer_hmg <- ifelse(apa_cosmic$Gene_symbol %in% pancan_genes, "Y", "N")

## 1032/1138 unique genes are in Pancancer HMG genes.
this <- apa_cosmic[, c("Ensembl_ID_s", "Pan_cancer_hmg")]
this <- unique(this)
dim(this)
table(this$Pan_cancer_hmg)


## overlapping with Xia et al. 2014 NC APA
xia_apa <- read.delim(file.path(indir, "DaPars.NC.dynamic.APA.txt"),
                      header = TRUE, as.is = TRUE)
xia_apa$Gene_symbol <- gsub("(NM_\\d+\\|)?(.+?)\\|.+", "\\2", xia_apa$Event_id, 
                            perl = TRUE)

apa_cosmic$Xia_APA_Genes <- ifelse(apa_cosmic$Gene_symbol %in% xia_apa$Gene_symbol, "Y", "N")
this <- apa_cosmic[, c("Ensembl_ID_s", "Xia_APA_Genes")]
this <- unique(this)
dim(this)
table(this$Xia_APA_Genes)

# N   Y 
# 946 192 


## lung cancer APA

lung_apa <- read.delim(file.path(indir, "Lung.cancer.APA.txt"),
                       header = FALSE, as.is = TRUE)

lung_apa$Gene_symbol <- gsub(" (repressed|enhanced)", "", lung_apa$V1, perl = TRUE)
lung_apa$Status <- gsub(".+ (repressed|enhanced)", "\\1_proximal", lung_apa$V1, perl = TRUE)
sum(apa_cosmic$Gene_symbol %in% lung_apa$Gene_symbol)

apa_cosmic <- merge(apa_cosmic, lung_apa, by= "Gene_symbol", all.x = TRUE )


library("WriteXLS")
WriteXLS(apa_cosmic, ExcelFileName = file.path(indir, "01.Annotate.apa.genes_3.xls"),
         row.names = FALSE,BoldHeaderRow = TRUE, FreezeRow = 1)

