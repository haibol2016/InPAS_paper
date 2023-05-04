#!/usr/bin/env Rscript
library("BSgenome.Hsapiens.Gencode.hsa38")
library("EnsDb.Hsapiens.v86")

library("devtools")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

library("optparse")

option_list = list(
    make_option("--sqlite_db", type="character",  
                help="path to the sqlite database", metavar="character"),
    make_option("--TxDb_file", type="character",
                help="path to a TxDb sqlite database", metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

sqlite_db <- opt$sqlite_db
outdir <- opt$outdir
TxDb_file <- opt$TxDb_file
TxDb <- loadDb(TxDb_file)
EnsDb <- EnsDb.Hsapiens.v86
## extract 3' UTR annotation from a TXDb
genome <- BSgenome.Hsapiens.Gencode.hsa38

set_globals(genome = genome,
            TxDb = TxDb,
            EnsDb = EnsDb,
            outdir = outdir,
            chr2exclude = c("chrM", "MT",
                            "Pltd", "chrPltd"),
            lockfile = tempfile())

eSet <- get_UTR3eSet(sqlite_db,
                     normalize = "none",
                     singleSample = FALSE)
saveRDS(eSet, file= file.path(outdir, "009.eSet.RDS"))
