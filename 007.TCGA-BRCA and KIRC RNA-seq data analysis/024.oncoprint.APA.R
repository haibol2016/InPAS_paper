
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\TCGA)")
files <- dir(".", "010.fisher.test.filtered.out.txt$", 
             recursive = TRUE, full.names = TRUE)
names(files) <- gsub(".+?/(.+?)/(.+)_010.+", "\\1_\\2", files, perl = TRUE)

dPDUI <- mapply(function(x, y){
    dat <- read.delim(x, header = TRUE, as.is = TRUE)
    dat <- dat[dat$PASS == "TRUE", ]
    dat <- dat[, c(6:7, 32)]
    dat$transcript <- paste(dat$gene, dat$transcript, sep = "_")
    dat <- dat[,-2]
    colnames(dat)[2] <- y
    dat
},files, names(files), SIMPLIFY = FALSE)

dPDUI_m <- Reduce(function(df1, df2){merge(df1, df2, 
                                           by = "transcript",
                                           all = TRUE)}, dPDUI)
rownames(dPDUI_m) <- dPDUI_m$transcript
dPDUI_m <- dPDUI_m[, -1]


## dPDUI = normal - tumor
dPDUI_cat <- lapply(dPDUI_m, function(x) { 
    x <- ifelse(!is.na(x), ifelse(x < 0, "lengthening", "shortening"), "") 
})

dPDUI_cat <- do.call("cbind", dPDUI_cat)
dPDUI_cat <- as.data.frame(dPDUI_cat)
rownames(dPDUI_cat) <- rownames(dPDUI_m)

## recurrent APA (10%)
dPDUI_cat_recurrent <- dPDUI_cat[rowSums(dPDUI_cat !="") > 18, ] 
dPDUI_m <- dPDUI_m[rowSums(dPDUI_cat !="") > 18, ]

dPDUI_cat <- dPDUI_cat_recurrent

freq <- do.call("rbind", lapply(dPDUI_cat, function(.x) {data.frame(table(.x))}))
freq <- freq[freq$.x !="", ]

freq$type <- gsub("_.+$", "", rownames(freq), perl = TRUE)
freq$sample<- gsub("\\.\\d+$", "", rownames(freq), perl = TRUE)

freq_short <- freq[freq$.x == "shortening", ]
freq_short <- freq_short[order(freq_short$type, freq_short$Freq), ]
freq_short_split <- split(freq_short, freq_short$type)
total_freq <- lapply(freq_short_split, function(.x) {sum(.x$Freq)})
sort(unlist(total_freq), decreasing =TRUE ) -> a


freq_short_split <- freq_short_split[names(a)]

freq_short_split <- lapply(freq_short_split, 
            function(x){x <- x[order(x$Freq, decreasing  = TRUE), ]})
freq_short <- do.call("rbind", freq_short_split)


dPDUI_cat <- dPDUI_cat[, freq_short$sample]


dPDUI_cat <- data.frame(t(dPDUI_cat))
short_freq <- lapply(dPDUI_cat, function(x) {
    dat <- as.data.frame(table(x))
    dat <- dat[-1, ,drop =FALSE]
    dat$x <- as.character(dat$x)
    if (!"shortening" %in% dat$x) { dat <- rbind(dat, data.frame(x = "shortening", Freq = 0))}
    if (!"lengthening" %in% dat$x) {dat <- rbind(dat, data.frame(x = "lengthening", Freq = 0))}
    dat
})
short_freq <- do.call("rbind", short_freq)
short_freq_short <- short_freq[short_freq$x == "shortening", ]
rownames(short_freq_short) <- gsub(".\\d+$", "", rownames(short_freq_short), perl = TRUE)
short_freq_long <- short_freq[short_freq$x == "lengthening", ]
rownames(short_freq_long) <- gsub(".\\d+$", "", rownames(short_freq_long), perl = TRUE)

colnames(short_freq_short)[2] <-  "shortening"
colnames(short_freq_long)[2] <-  "lengthening"
short_freq <- cbind(short_freq_short, short_freq_long)
short_freq <- with(short_freq, short_freq[order(-shortening, lengthening), ])

dPDUI_cat <- as.data.frame(t(dPDUI_cat))
dPDUI_cat <- dPDUI_cat[rownames(short_freq), ]
dPDUI_cat <- as.matrix(dPDUI_cat)

library("ComplexHeatmap")
pdf("Fig 1.Sorted.oncoprint.left.shortening.0.10.recurrent.rate.pdf", height = 5, width = 8)
col <- c(shortening = "red", lengthening = "blue")
cols <- rev(rainbow(2))
names(cols) <- unique(gsub("_.+$", "", colnames(dPDUI_cat),
                                                     perl = TRUE))
oncoPrint(dPDUI_cat,
          alter_fun = list(
              background = function(x, y, w, h)  NULL,
              shortening = function(x, y, w, h) {
                  grid.rect(x, y, w*0.5, h*0.2, 
                            gp = gpar(fill = col["shortening"], col = NA))
              },
              lengthening = function(x, y, w, h) {
                  grid.rect(x, y, w*0.5, h*0.2, 
                            gp = gpar(fill = col["lengthening"], col = NA)) 
              } 
          ), 
          top_annotation = HeatmapAnnotation(
              cbar =
                  anno_oncoprint_barplot(border = TRUE,  which = "column",
                                         bar_width = 0.5,
                                         height = unit(1.0, "cm")),
              # foo = anno_text(gsub("_.+$", "", colnames(dPDUI_cat), perl = TRUE),
              #                 location = 1, rot = 30,
              #                 just = "right", gp = gpar(fontsize = 8)),
              type = factor(gsub("_.+$", "", colnames(dPDUI_cat),
                                           perl = TRUE), 
                            levels = c("BRCA", "KIRC")),
              simple_anno_size = unit(0.3, "cm"),
              show_annotation_name = c(type = FALSE),
              col = list(type = cols)
          ),
          # right_annotation = rowAnnotation(
          #     rbar =
          #         anno_oncoprint_barplot(c("lengthening", "shortening"),
          #                                show_fraction = TRUE,
          #                                border = TRUE,
          #                                height = unit(1.0, "cm"),
          #                                axis_param = list(side = "top",
          #                                                  labels_rot = 90)),
          # ),
          left_annotation = rowAnnotation(
             lbar =
                  anno_oncoprint_barplot("shortening",
                                         show_fraction = TRUE,
                                         border = TRUE,
                                         height = unit(0.5, "cm"),
                                         axis_param = list(side = "top",
                                                           labels_rot = 90,
                                                           direction = "reverse"))
          ),
          #alter_fun_is_vectorized = FALSE,
          remove_empty_columns = TRUE, 
          #remove_empty_rows = TRUE,
          show_row_names = FALSE,
          column_order = colnames(dPDUI_cat),
          row_order = rownames(dPDUI_cat),
          show_pct = FALSE,
          use_raster = FALSE,
          show_column_names = FALSE,
          column_names_gp = gpar(fontsize = 6),
          heatmap_legend_param = list(title = r"(3' UTR)"),
          col = col)
dev.off()

write.table(dPDUI_cat, file = "00.0.10.recurrent.APA.events.across.2.cancer.types.txt",
            row.names =TRUE, sep = "\t", quote =FALSE)

save.image(file = "recurrent.APA.events.across.2.cancer.types.RData")
load("recurrent.APA.events.across.2.cancer.types.RData")

## KEGG pathway and REACTOME pathway analysis
apa_events <- dPDUI_cat

expr <- read.delim("00.all.expressed.genes.txt",
                   header = FALSE, as.is = TRUE)
rownames(expr) <- gsub("\\.\\d+", "", expr[,1], perl = TRUE)


library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("dplyr")
library("biomaRt")

ensembl = useEnsembl(biomart="ensembl",
                     dataset="hsapiens_gene_ensembl", 
                     version=107)
human_entrez_id <- getBM(attributes=c('ensembl_gene_id',
                                      'entrezgene_id'),
                         filters = 'ensembl_gene_id',
                         values =  rownames(expr), mart = ensembl)

apa_genes <- rownames(apa_events)
apa_genes <- gsub("\\.\\d+_.+", "", apa_genes,
                  perl = TRUE)

## 61
all_shorten_genes <- data.frame(gene = apa_genes[do.call("c", 
                                                         lapply(as.data.frame(t(apa_events)), 
                                                                function(.x) {
                                                                    all(.x =="" | .x =="shortening")
                                                                }))])


## 3
all_lengthen_genes <- data.frame(gene = apa_genes[do.call("c", 
                                                          lapply(as.data.frame(t(apa_events)), 
                                                                 function(.x) {
                                                                     all(.x =="" | .x =="lengthening")
                                                                 }))])

## 938
mix_genes <- data.frame(gene = apa_genes[!apa_genes %in% c(all_shorten_genes[, 1], 
                                                           all_lengthen_genes[, 1])])

## 938
any_shorten_gene <-  data.frame(gene = apa_genes[do.call("c", 
                                                    lapply(as.data.frame(t(apa_events)), 
                                                           function(.x) {
                                                               any(.x =="shortening")
                                                           }))])
## 880
any_lengthen_gene <-  data.frame(gene = apa_genes[do.call("c", 
                                                         lapply(as.data.frame(t(apa_events)), 
                                                                function(.x) {
                                                                    any(.x =="lengthening")
                                                                }))])
apa_gene_class <- list(short = all_shorten_genes,
                       long = all_lengthen_genes,
                       mix = mix_genes,
                       any_short = any_shorten_gene,
                       any_lengthen = any_lengthen_gene,
                       all = data.frame(gene = apa_genes))  ## 941
apa_gene_class <- lapply(apa_gene_class, function(.x){
    temp <- merge(.x, human_entrez_id, 
                  by.x = "gene", 
                  by.y = "ensembl_gene_id", all.x = TRUE)
    temp <- temp[!is.na(temp$entrezgene_id), ]
    temp
})


save(human_entrez_id, apa_gene_class, file = "Two.cancer.recurrent.sig.APA.RData")


library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("dplyr")
#library("biomaRt")
load("Two.cancer.recurrent.sig.APA.RData")
for (i in seq_along(apa_gene_class))
{
    kk <- enrichKEGG(
        gene = as.character(unique(apa_gene_class[[i]]$entrezgene_id)),
        organism = "hsa",
        keyType = "ncbi-geneid",
        pAdjustMethod = "BH",
        universe = as.character(human_entrez_id$entrezgene_id[!is.na(human_entrez_id$entrezgene_id)]),
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        use_internal_data = FALSE)
    
    df <- as.data.frame(kk)
    if (nrow(df) > 0) {
        write.table(kk, file = paste0(names(apa_gene_class)[i], "-KEGG pathway.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        pdf(paste0(names(apa_gene_class)[i],"-KEGG pathway.pdf"),
            width = 8, height = 6)
        print((dotplot(kk)))
        dev.off()
    }
    
    reactome <- enrichPathway(gene = as.character(unique(apa_gene_class[[i]]$entrezgene_id)),
                              organism = "human",
                              pAdjustMethod = "BH",
                              qvalueCutoff = 0.2,
                              universe = as.character(human_entrez_id$entrezgene_id[!is.na(human_entrez_id$entrezgene_id)]),
                              minGSSize = 10,
                              maxGSSize = 500,
                              readable = TRUE)
    df <- as.data.frame(reactome)
    if (nrow(df) > 0) {
        write.table(reactome, file = paste0(names(apa_gene_class)[i], "-Reactome pathway.txt"),
                    sep = "\t", quote = FALSE, row.names = FALSE)
        pdf(paste0(names(apa_gene_class)[i],"Reactome pathway.pdf"),
            width = 14, height = 8)
        print(dotplot(reactome))
        dev.off()
    }
}  


## 

apa_genes <- data.frame(Gene_ID = gsub("_.+", "", rownames(apa_events), perl = TRUE), 
                        Tx_ID = gsub(".+_", "", rownames(apa_events), perl = TRUE))
gencode_gene <- read.delim("../005.gencode.Tx.gene.name.type.txt", header = FALSE)
apa_genes <- merge(apa_genes, gencode_gene, by.x = "Tx_ID",
                   by.y = "V2", all.x = TRUE)

## APA genes compared to cancer genes downloaded from COSMIC
indir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types/)" 
cosmic <- read.delim(file.path(indir, "COSMIC cancer_gene_census.txt"), 
                     header = TRUE, as.is = TRUE)
cosmic$Ensembl_ID <- gsub(".+?(ENSG[^,]+).+", "\\1", cosmic$Synonyms, perl = TRUE)
cosmic$Ensembl_ID_s <- gsub("\\.\\d+", "", cosmic$Ensembl_ID, perl = TRUE)
apa <- apa_genes
apa$Ensembl_ID_s <- gsub("\\.\\d+", "", apa$Gene_ID, perl = TRUE)
apa <- apa[, -c(2:3)]
colnames(apa)[2:3] <- c("Gene_symbol", "Biotype")
## merge
apa_cosmic <- merge(apa, cosmic, by = "Ensembl_ID_s", all.x = TRUE)

colnames(apa_cosmic)

# apa_cosmic <- apa_cosmic[ , -c(3,4,7,29)]
# 
# 
# ## add shortening/lengthening frequency
# apa_freq <- read.delim(file.path(indir,"00.10percen.recurrent.APA.events.across.9.cancer.types.txt"),
#                        header = TRUE, as.is = TRUE)
# 
# freq <- lapply(as.data.frame(t(apa_freq)), function(.x) {
#     dat <- as.data.frame(table(.x))
#     dat$.x <- as.character(dat$.x)
#     dat$.x <- ifelse(dat$.x == "", "nsc", dat$.x)
#     
#     if (!"nsc" %in% dat$.x){dat <- rbind(dat, data.frame(.x = "nsc", Freq = 0 ))}
#     else if (!"shortening" %in% dat$.x) {dat <- rbind(dat, data.frame(.x = "shortening", Freq = 0 ))}
#     else if (!"lengthening" %in% dat$.x)  {dat <- rbind(dat, data.frame(.x = "lengthening", Freq = 0 ))}
#     dat <- dat[order(dat$.x), ]
#     rownames(dat) <- dat[, 1]
#     dat <- dat[, -1, drop = FALSE]
#     dat <- t(dat)
#     dat/50
# })
# freq <- do.call(rbind, freq)
# rownames(freq) <- rownames(apa_freq)
# rownames(freq) <- gsub("ENSG\\d+\\.\\d+_", "", rownames(freq), perl = TRUE)
# 
# ## merge with apa_cosmic
# 
# apa_cosmic <- merge(freq, apa_cosmic, by.x = "row.names", 
#                     by.y = "Tx_ID", all.x = TRUE)
# 
# ## pancancer.sig.mutated.genes

pancan <- read.delim(file.path(indir, "Pancancer.significantly.mutated.genes.txt"),
                     header = TRUE, as.is = TRUE)

pancan_genes <- unique(pancan[, 1])
apa_cosmic$Pan_cancer_hmg <- ifelse(apa_cosmic$Gene_symbol %in% pancan_genes, "Y", "N")

## 839/922 unique genes are in Pancancer HMG genes.
this <- apa_cosmic[, c("Ensembl_ID_s", "Pan_cancer_hmg")]
this <- unique(this)
dim(this)
table(this$Pan_cancer_hmg)

#N   Y 
#83 839 

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
# 506 416 


## lung cancer APA

lung_apa <- read.delim(file.path(indir, "Lung.cancer.APA.txt"),
                       header = FALSE, as.is = TRUE)

lung_apa$Gene_symbol <- gsub(" (repressed|enhanced)", "", lung_apa$V1, perl = TRUE)
lung_apa$Status <- gsub(".+ (repressed|enhanced)", "\\1_proximal", lung_apa$V1, perl = TRUE)
sum(apa_cosmic$Gene_symbol %in% lung_apa$Gene_symbol) # 422

apa_cosmic <- merge(apa_cosmic, lung_apa, by= "Gene_symbol", all.x = TRUE )


library("WriteXLS")
WriteXLS(apa_cosmic, ExcelFileName = file.path("./", "01.Annotate.apa.genes_3.xls"),
         row.names = FALSE,BoldHeaderRow = TRUE, FreezeRow = 1)



### IGV visualization

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\TCGA\KIRC)")

genes <- read.delim("TCGA-CZ-5988_010.fisher.test.filtered.out.txt", header = TRUE, as.is = TRUE)

genes <- genes[genes$PASS == "TRUE", 9]
sock <- SRAdb::IGVsocket()
for (i in genes) {
    SRAdb::IGVgoto(sock, region = i)
    takeSnapshot <- readline(prompt="Take a snapshot? Y or N \n")
    if (grepl("Y|y", takeSnapshot, perl =TRUE)) {
        SRAdb::IGVsnapshot(sock) 
    }
}


