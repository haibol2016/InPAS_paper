
## download TCGA data
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\TCGA)")
manifest <- read.delim("gdc_manifest_20221003_135230.txt", header = TRUE)
metadata <- read.delim("gdc_sample_sheet.2022-10-03.tsv", header = TRUE)

## normal-tumor pairs
metadata <- metadata[metadata$Case.ID %in% metadata[metadata$Sample.Type == "Solid Tissue Normal", ]$Case.ID, ]

manifest <- manifest[manifest$filename %in% metadata$File.Name, ]

#write.table(manifest, file = "TCGA_KRIC_BRCA.manifest.2.download.txt",
#            sep = "\t", quote =FALSE, row.names = FALSE)
manifest <- manifest[order(manifest$id), ]
metadata <- metadata[order(metadata$File.ID), ]

manifest <- cbind(manifest, metadata)
write.table(manifest, file = "TCGA_KRIC_BRCA.manifest.metadata.txt",
            sep = "\t", quote =FALSE, row.names = FALSE)

BRCA <- manifest[manifest$Project.ID == "TCGA-BRCA", 1:5]
KIRC <- manifest[manifest$Project.ID == "TCGA-KIRC", 1:5]

write.table(BRCA, file = "TCGA_BRCA.manifest.metadata.txt",
            sep = "\t", quote =FALSE, row.names = FALSE)

write.table(KIRC, file = "TCGA_KIRC.manifest.metadata.txt",
            sep = "\t", quote =FALSE, row.names = FALSE)
