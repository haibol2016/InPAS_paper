#!/usr/bin/env Rscript

library("devtools")
library("BSgenome.Hsapiens.Gencode.hsa38")
library("EnsDb.Hsapiens.v86")

devtools::load_all("../InPAS-07062022")

library("optparse")

option_list = list(
    make_option("--utr3", type="character",
                help="path to a single .RDS file", metavar="character"),
    make_option("--sqlite_db", type="character",  
                help="path to the sqlite database", metavar="character"),
    make_option("--TxDb_file", type="character",
                help="path to a TxDb sqlite database", metavar="character"),
    make_option("--seqname", type="character",
                help="chromosome ID", metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"))

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser)

sqlite_db <- opt$sqlite_db
TxDb_file <- opt$TxDb_file
utr3 <- readRDS(opt$utr3)
seqname <- opt$seqname
chr.utr3 <- utr3[[seqname]]
outdir <- opt$outdir

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

null <- get_regionCov(chr.utr3 = chr.utr3,
              sqlite_db = sqlite_db,
              outdir = outdir,
              phmm = FALSE)
