#!/usr/bin/env Rscript

library("devtools")
library("BSgenome.Hsapiens.Gencode.hsa38")
library("EnsDb.Hsapiens.v86")
library("optparse")

devtools::load_all("../InPAS-07062022")

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

seqname <- opt$seqname
sqlite_db <- opt$sqlite_db
utr3 <- readRDS(opt$utr3)
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


## load the Naive Bayes classifier model from the cleanUpdTSeq package
library(cleanUpdTSeq)
data(classifier)
genome <- BSgenome.Hsapiens.Gencode.hsa38

null <- search_CPs(seqname = seqname,
                       sqlite_db = sqlite_db,
                       genome = getInPASGenome(), 
                       MINSIZE = 10, 
                       window_size = 100,
                       search_point_START = 1,
                       search_point_END = NA,
                       cutEnd = NA,
                       filter.last = TRUE,
                       adjust_distal_polyA_end = FALSE,
                       long_coverage_threshold = 5,
                       PolyA_PWM = NA, 
                       classifier = classifier,
                       classifier_cutoff = 0.8,
                       shift_range = 100,
                       step = 2,
                       outdir = getInPASOutputDirectory(),
                       silence = FALSE,
                       cluster_type = "lsf",
                       template_file = "./batchtools.lsf.tmpl",
                       mc.cores = 8,
                       future.chunk.size = 100,
                       resources = list(walltime = 3600*8, ncpus = 8,
                                        mpp = 1024*4, queue = "long",
                                        memory = 4 * 8 * 1024),
                      DIST2ANNOAPAP = 200,
                      DIST2END = 200)
