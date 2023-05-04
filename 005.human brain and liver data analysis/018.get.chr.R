#!/usr/bin/env Rscript
library("devtools")
library("optparse")

devtools::load_all("/home/hl84w/arthur_mercurio/Haibo/InPAS-07062022")

option_list = list(
    make_option("--utr3", type="character", 
                help="path to UTR3 RDS",
                metavar="character"),
    make_option("--sqlite_db", type="character",  
                help="path to a sqlite database", 
                metavar="character"),
    make_option("--outdir", type="character", 
                help="output directory", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)
outdir <- opt$outdir
sqlite_db <- opt$sqlite_db
utr3 <- readRDS(opt$utr3)

chr <- get_chromosomes(utr3, sqlite_db)
write.table(chr, file=file.path(outdir, "02.seqnames.txt"),
            row.names = FALSE, col.names = FALSE,
            quote = FALSE)
