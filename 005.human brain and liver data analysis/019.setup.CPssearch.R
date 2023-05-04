#!/usr/bin/env Rscript
library("devtools")
library("RSQLite")
library("EnsDb.Hsapiens.v86")
library("BSgenome.Hsapiens.Gencode.hsa38")
library("optparse")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

option_list = list(
    make_option("--utr3", type="character",  
                help="path to a utr3 RDS file", 
                metavar="character"),
    make_option("--sqlite_db", type="character",  
                help="path to a sqlite database file", 
                metavar="character"),
    make_option("--seqname", type="character",  
                help="chr name", 
                metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"),
    make_option("--TxDb", type="character", default=NULL, 
                help="TXDb sqlite file name", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
outdir <- opt$outdir
TxDb_file <- opt$TxDb
sqlite_db <- opt$sqlite_db
seqname <- opt$seqname
utr3 <- readRDS(opt$utr3)

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
null <- setup_CPsSearch (sqlite_db = sqlite_db,
                            genome = getInPASGenome(), 
                            chr.utr3 = utr3[[seqname]],
                            seqname = seqname,
                            background = "same_as_long_coverage_threshold",
                            TxDb = getInPASTxDb(),
                            hugeData = TRUE,
                            outdir = getInPASOutputDirectory(),
                            silence = FALSE,
                            minZ = 2,
                            cutStart = 0,
                            MINSIZE = 10,
                            coverage_threshold = 30)

