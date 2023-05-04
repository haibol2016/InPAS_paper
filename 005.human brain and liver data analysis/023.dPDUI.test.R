#!/usr/bin/env Rscript
library("BSgenome.Hsapiens.Gencode.hsa38")
library("EnsDb.Hsapiens.v86")
library("devtools")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

library(limma)
library("optparse")

option_list = list(
    make_option("--sqlite_db", type="character",  
                help="path to the sqlite database", metavar="character"),
    make_option("--eSet", type="character",
                help="path to eSet RDS file", metavar="character"),
    make_option("--TxDb_file", type="character",
                help="path to a TxDb sqlite database", metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

eSet <- readRDS(opt$eSet)
sqlite_db <- opt$sqlite_db
outdir <- opt$outdir
TxDb_file <- opt$TxDb_file

EnsDb <- EnsDb.Hsapiens.v86
## extract 3' UTR annotation from a TXDb
TxDb <- loadDb(TxDb_file)
genome <- BSgenome.Hsapiens.Gencode.hsa38

set_globals(genome = genome,
            TxDb = TxDb,
            EnsDb = EnsDb,
            outdir = outdir,
            chr2exclude = c("chrM", "MT",
                            "Pltd", "chrPltd"),
            lockfile = tempfile())

db_conn <- dbConnect(drv = RSQLite::SQLite(), dbname= sqlite_db)
metadata <- dbReadTable(db_conn, "metadata")
dbDisconnect(db_conn)
metadata$condition <- factor(metadata$condition)

groups <- split(metadata$tag, metadata$condition)
gp1 <- groups[[1]]
gp2 <- groups[[2]]

design <- model.matrix(~0 + condition, data = metadata)
colnames(design) <- gsub("condition", "", colnames(design))
contrast.matrix <- makeContrasts(contrasts = paste0(names(groups)[1], "-", names(groups)[2]), 
                                 levels = design)
save(design, contrast.matrix, file =file.path(outdir, "010.design.contrast.RData"))

fisher_test_out <- test_dPDUI(eset = eSet, 
                       method = "fisher.exact",
                       normalize = "none",
                       sqlite_db = sqlite_db)
saveRDS(fisher_test_out, file = file.path(outdir, "010.fisher.test.out.RDS"))

fisher_filter_out <- filter_testOut(res = fisher_test_out,
                             gp1 = gp1,
                             gp2 = gp2,
                             adj.P.Val_cutoff = 0.05,
                             dPDUI_cutoff = 0.2,
                             PDUI_logFC_cutoff = log2(1.5))
write.table(fisher_filter_out, file = file.path(outdir, "010.fisher.test.filtered.out.txt"),
           sep = "\t", quote = FALSE, row.names = FALSE)


