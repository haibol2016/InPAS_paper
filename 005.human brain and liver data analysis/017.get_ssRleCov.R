#!/usr/bin/env Rscript
library("devtools")
library("RSQLite")
library("EnsDb.Hsapiens.v86")
library("BSgenome.Hsapiens.Gencode.hsa38")
library("optparse")
library("future.apply")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

option_list = list(
    make_option("--bedgraph", type="character",  
                help="path to a bedgraph file", 
                metavar="character"),
    make_option("--tag", type="character",  
                help="a name tag of a sample", 
                metavar="character"),
    make_option("--sqlite_db", type="character",  
                help="path to a sqlite database", 
                metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"),
    make_option("--TxDb", type="character", default=NULL, 
                help="TXDb sqlite file name", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
outdir <- opt$outdir
TxDb_file <- opt$TxDb
tag <- opt$tag
bedgraph <- opt$bedgraph
sqlite_db <- opt$sqlite_db

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

plan(multicore)    
null <- get_ssRleCov(bedgraph = bedgraph,
             tag = tag,
             genome = getInPASGenome(),
             sqlite_db = sqlite_db,
             future.chunk.size = 10,
             outdir = getInPASOutputDirectory(),
             chr2exclude = getChr2Exclude())
