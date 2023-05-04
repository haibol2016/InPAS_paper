#!/usr/bin/env Rscript
library("devtools")
library("RSQLite")
library("EnsDb.Hsapiens.v86")
library("BSgenome.Hsapiens.Gencode.hsa38")
library("optparse")
library("future")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

option_list = list(
    make_option("--metadata", type="character",  
                help="path to a single metadata file", 
                metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"),
    make_option("--TxDb", type="character", default=NULL, 
                help="TXDb sqlite file name", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
metadata <- opt$metadata
outdir <- opt$outdir
TxDb_file <- opt$TxDb

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
sqlite_db <- setup_sqlitedb(metadata, outdir)

plan(multicore)
extract_UTR3Anno(sqlite_db,
                 TxDb = getInPASTxDb(),
                 edb = getInPASEnsDb(),
                 genome = getInPASGenome(),
                 outdir = getInPASOutputDirectory(),
                 chr2exclude = getChr2Exclude(),
                 MAX_EXONS_GAP = 10000L,
                 UTR3 = "last.exon",
                 future.chunk.size = 4)

