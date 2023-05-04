
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types)")
files <- dir("./InPAS.APA.events", ".010.fisher.test.filtered.out.txt$", 
             recursive = TRUE, full.names = TRUE)
names(files) <- gsub(".+/(.+)\\.010.+", "\\1", files, perl = TRUE)

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
dPDUI_cat_recurrent <- dPDUI_cat[rowSums(dPDUI_cat !="") >= 5, ] 
dPDUI_cat <- dPDUI_cat_recurrent
dPDUI_m <- dPDUI_m[rowSums(dPDUI_cat !="") >= 5, ]

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
cols <- rev(rainbow(9))
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
                            levels = c("GAC", "HCC","LUAD", "LUSC",
                                       "CSCC", "ESCC", "SCLC", "SRCC", "PTC")),
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
          show_column_names = TRUE,
          column_names_gp = gpar(fontsize = 6),
          heatmap_legend_param = list(title = r"(3' UTR)"),
          col = col)
dev.off()

write.table(dPDUI_cat, file = "00.0.10.recurrent.APA.events.across.9.cancer.types.txt",
            row.names =TRUE, sep = "\t", quote =FALSE)

save.image(file = "recurrent.APA.events.across.9.cancer.types.RData")
load("recurrent.APA.events.across.9.cancer.types.RData")

## KEGG pathway and REACTOME pathway analysis
apa_events <- read.delim(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types/00.0.10.recurrent.APA.events.across.9.cancer.types.txt)")

expr <- read.delim(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types\9.cancer.final.count.table/Table 2. sva-adjusted.gene.expression.txt)",
                   header = TRUE, as.is = TRUE)
rownames(expr) <- gsub("\\.\\d+", "", rownames(expr), perl = TRUE)


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

## 101
all_shorten_genes <- data.frame(gene = apa_genes[do.call("c", 
                                                         lapply(as.data.frame(t(apa_events)), 
                                                                function(.x) {
                                                                    all(.x =="" | .x =="shortening")
                                                                }))])


## 13
all_lengthen_genes <- data.frame(gene = apa_genes[do.call("c", 
                                                          lapply(as.data.frame(t(apa_events)), 
                                                                 function(.x) {
                                                                     all(.x =="" | .x =="lengthening")
                                                                 }))])

## 1044
mix_genes <- data.frame(gene = apa_genes[!apa_genes %in% c(all_shorten_genes[, 1], 
                                                           all_lengthen_genes[, 1])])

## 1149
any_shorten_gene <-  data.frame(gene = apa_genes[do.call("c", 
                                                    lapply(as.data.frame(t(apa_events)), 
                                                           function(.x) {
                                                               any(.x =="shortening")
                                                           }))])
## 1061
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
                       all = data.frame(gene = apa_genes))  ## 1162
apa_gene_class <- lapply(apa_gene_class, function(.x){
    temp <- merge(.x, human_entrez_id, 
                  by.x = "gene", 
                  by.y = "ensembl_gene_id", all.x = TRUE)
    temp <- temp[!is.na(temp$entrezgene_id), ]
    temp
})


save(human_entrez_id, apa_gene_class, file = "nine.cancer.recurrent.sig.APA.RData")


library("ReactomePA")
library("clusterProfiler")
library("magrittr")
library("dplyr")
#library("biomaRt")
load("nine.cancer.recurrent.sig.APA.RData")
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
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types)")

apa_events <- read.delim(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Nine_cancer_types/00.0.10.recurrent.APA.events.across.9.cancer.types.txt)")

#lengthening  shortening 
#    3174        5156

apa_genes <- data.frame(Gene_ID = gsub("_.+", "", rownames(apa_events), perl = TRUE), 
                        Tx_ID = gsub(".+_", "", rownames(apa_events), perl = TRUE))
gencode_gene <- read.delim("../005.gencode.Tx.gene.name.type.txt", header = FALSE)
apa_genes <- merge(apa_genes, gencode_gene, by.x = "Tx_ID",
                   by.y = "V2", all.x = TRUE)



