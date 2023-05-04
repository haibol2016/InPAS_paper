library("EnhancedVolcano")
library("ggplot2")
library("ggrepel")
library("WriteXLS")
#options(ggrepel.max.overlaps = Inf)

## HeLa RNA-seq data analyzed using DaPars

in_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa)"
dapars <- read.delim(file.path(in_dir,
                               "01.HeLa.Github.DaPars.Final.sig.APA.txt"),
                     header = TRUE, as.is = TRUE
                     )

dapars_pos <- dapars[dapars$Strand == "+", ]  ## 2241
dapars_pos <- dapars_pos[order(dapars_pos$Gene_id, dapars_pos$start, -dapars_pos$end), ]
dapars_pos <- dapars_pos[!duplicated(dapars_pos[, c(2,10)]), ]


## uniq start
dapars_pos_uniq <- unique(dapars_pos[, c(2,10)])
nrow(dapars_pos_uniq) # 1545

dapars_neg <- dapars[dapars$Strand == "-", ]  ## 2004
dapars_neg <- dapars_neg[order(dapars_neg$Gene_id, dapars_neg$end, dapars_neg$start), ]
dapars_neg <- dapars_neg[!duplicated(dapars_neg[, c(2,11)]), ]

dapars_neg_uniq <- unique(dapars_neg[, c(2,11)])
nrow(dapars_neg_uniq) # 1403


## uniq entries
dapars <- rbind(dapars_pos, dapars_neg)  ## 2948
write.table(dapars, file = file.path(in_dir, 
                                     "01.HeLa.Github.DaPars.Final.sig.APA.deduplicates.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)


## InPAS fisher exact test
in_dir_inpas <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa\InPAS.last.exon)"
inpas <- read.delim(file.path(in_dir_inpas, "010.InPAS.sig.APA.last.exon.txt"),
                    header = TRUE, as.is = TRUE)
inpas <- inpas[!is.na(inpas$start), ]  ## 3105



library(VennDiagram)
library(RColorBrewer)
library(ggVennDiagram)
library(ggplot2)
x = list(InPAS = paste(inpas$gene, inpas$transcript, sep = "_"), 
         DaPars = paste(dapars$Gene_id, dapars$Transcript, sep = "_"))
venn.diagram(
  x,
  category.names = c("InPAS" , "DaPars"),
  #filename = file.path(in_dir, "HeLa.APA.InPAS.DaPars.svg"),
  #output = TRUE ,
  #imagetype="svg" ,
  lwd = 1,
  col=c("#440154ff", '#21908dff'),
  fill = c(alpha("#440154ff",0.3), alpha('#21908dff',0.3)),
  cex = 0.5
)
library(UpSetR)
upset(fromList(x), order.by = "freq")
write.table(paste(inpas$gene, inpas$transcript, sep = "_"), 
            file = file.path(in_dir, "InPAS.APA.ids.txt"),
            sep = "\t", quote =FALSE,
            row.names =FALSE, col.names = FALSE)
write.table(paste(dapars$Gene_id, dapars$Transcript, sep = "_"), 
            file = file.path(in_dir, "DaPars.APA.ids.txt"),
            sep = "\t", quote =FALSE,
            row.names =FALSE, col.names = FALSE)


## DeepPasta prediction versus polyAsite (https://polyasite.unibas.ch/atlas#2)

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa)")

plya <- read.delim("001.intesect.PAS.polyAsite2.0.txt",
                   header = TRUE, as.is = TRUE)

library(reshape2)

plya <- melt(plya, id.vars=c("cutoff"))

ggplot(plya, aes(x= cutoff, y = value, color = variable)) +
  geom_point() +
  geom_line()
  

## compare CP sites predicted by InPAS and DaPars relative to the PAC-seq polyA sites
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\HeLa)")

distance <- read.delim("001.InPAS.DaPars.pAPA.dist.diff.txt", header = TRUE)


