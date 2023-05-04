library(BSgenome)
library(devtools)
library(optparse)

## get per-chromosome fasta file and compressed in .gz format:
## chr1.fa.gz, chr2.fa.gz,...
generate_multifasta <- function(path_single_fasta_genome, out_dir)
{
  if (missing(path_single_fasta_genome) || missing(out_dir))
  {
    stop("path_single_fasta_genome and out_dir are required!")
  }
  if (!file.exists(path_single_fasta_genome)){
    stop("path_single_fasta_genome doesn't exists!")
  }
  if (!dir.exists(out_dir))
  {
    dir.create(out_dir, recursive = TRUE)
  }
  
  if (grepl(".(fa|fasta).gz$", path_single_fasta_genome))
  {
    in_fasta <- gzfile(path_single_fasta_genome, open ="rt")
  } else if (grepl(".(fa|fasta)$", path_single_fasta_genome)) {
    in_fasta <- file(path_single_fasta_genome, open = "r")
  } else {
    stop("It seems the genome sequence file is not a fasta file ",
         "which should with an extension .fa, .fasta, .fa.gz or .fasta.gz")
  }
  
  f <- ""
  while (length({line <- readLines(in_fasta, n = 1, warn = FALSE)}) > 0)
  {
    line <- trimws(line)
    if (grepl("^>", line))
    {
      f <- gzfile(file.path(out_dir, 
                            gsub("^>(chr)?([^\\s]+).*", "chr\\2.fa.gz",
                                 line, perl = TRUE)),
                  "w")
      writeLines(gsub("^>(chr)?([^\\s]+).*",">chr\\2", line, perl = TRUE), f)
    } else {
      writeLines(line, f)
    }
  }
  close(f)
  close(in_fasta)
  out_dir
}

## generate a seed file for building a BSgenome package
generate_seed_file <- function(path_to_multifasta, 
                               latin_name, 
                               common_name, 
                               genome_build, 
                               seed_file_name, 
                               fasta_url, 
                               source = "Ensembl", 
                               version = "1.0.0")
{
  if (missing(path_to_multifasta) || missing(latin_name) ||
      missing(common_name) || missing(genome_build) || 
      missing(seed_file_name) || missing(fasta_url))
  {
    stop("All arguments except source and version are required!")
  }
  if (!dir.exists(path_to_multifasta)) {
    stop("Path to multifasta ", path_to_multifasta, " doesn't exist!")
  }
  chr_fa_files <- dir(path_to_multifasta, ".fa.gz$")
  if (length(chr_fa_files) < 1)
  {
    stop("There is no multiple fasta files in the directory ", 
         path_to_multifasta, "!")
  }
  seed_dir <- dirname(path_to_multifasta)
  if (!dir.exists(seed_dir))
  {
    dir.create(seed_dir, recursive = TRUE)
  }
  sink(seed_file_name)
  BSgenomeObjname <- gsub("^(.).*\\s+(.+)", "\\1\\2", latin_name)
  package_name <- paste("BSgenome", BSgenomeObjname, source, genome_build, sep = ".")
  cat(paste0("Package: ", package_name, "\n"))
  cat(paste("Title: Full genome sequences for", latin_name,
            paste0("(", source), "version", paste0(genome_build,")\n")))
  cat(paste("Description: Full genome sequences for", latin_name,
            paste0("(", common_name, ")"), "as provided by",
            source, paste0("(", genome_build, ")"), 
            "and stored in Biostrings objects.\n"))
  cat(paste0("Version: ", version, "\n"))
  cat(paste0("organism: ", latin_name, "\n"))
  cat(paste0("common_name: ", common_name,"\n"))
  cat(paste0("provider: ", source, "\n"))
  cat("release_date: May, 2007\n")
  cat(paste0("genome: ", genome_build, "\n"))
  cat(paste0("source_url: ", fasta_url, "\n"))
  cat(paste0("BSgenomeObjname: ", BSgenomeObjname, "\n"))
  cat(paste0("organism_biocview: ", 
             gsub("\\s+", "_", latin_name, perl = TRUE), "\n"))
  chromosome_names <- gsub(".fa.gz$", "", dir(path_to_multifasta, "fa.gz$"))
  cat(paste0('seqnames: c("', paste(chromosome_names, collapse = '","'), '")\n'))
  
  circ_seqs <- c("chrM", "MT", "Pltd", "chrPltd")
  circ_seqs <- circ_seqs[circ_seqs %in% chromosome_names]
  cat(paste0('circ_seqs: c("', paste(circ_seqs, collapse = '","'), '")\n'))
  cat(paste0("seqs_srcdir: ", path_to_multifasta, "\n"))
  cat(paste0("seqfiles_suffix: .fa.gz\n"))
  sink()
  c(package_name, seed_file_name)
}

forge_install_BSgenome <- function(seed_file = NULL, 
                                   dest_dir = ".")
{
  if (is.null(seed_file))
  {
    stop("seed_file are required!")
  }
  if (!dir.exists(dest_dir))
  {
    dir.create(dest_dir, recursive = TRUE)
  }
  pkgname <- gsub("Package:\\s*([^\\s]+)", 
                  "\\1",
                  trimws(readLines(seed_file, n =1)))
  
  ## crate a BSgenome package from a seed file
  forgeBSgenomeDataPkg(seed_file, destdir = dest_dir)
  
  ## check, build and install package
  ## OR install using command line
  # R CMD build BSgenome.Hsapiens.Ensembl.GRCh38
  # R CMD check BSgenome.Hsapiens.Ensembl.GRCh38_1.0.0.tar.gz
  # R CMD INSTALL BSgenome.Hsapiens.Ensembl.GRCh38_1.0.0.tar.gz
  devtools::check(file.path(dest_dir, pkgname))
  devtools::build(file.path(dest_dir, pkgname))
  devtools::install(file.path(dest_dir, pkgname))
}


## entry point

option_list = list(
  make_option("--path_single_fasta_genome", type="character",  
              help="path to the reference genome fasta file", 
              metavar="character"),
  make_option("--out_dir", type="character", 
              help="output directory", metavar="character"),
  make_option("--latin_name", type="character",  
              help="the Latin name of the species", 
              metavar="character"),
  make_option("--common_name", type="character",  
              help="the commonly used name of the species", 
              metavar="character"),  
  make_option("--genome_build", type="character",  
              help="the genome build version", 
              metavar="character"),
  make_option("--seed_file_name", type="character", default=NULL, 
              help="the seed file name to generate", metavar="character"),
  make_option("--fasta_url", type="character", default=NULL, 
            help="URL for downloading the reference genome fasta file", metavar="character"),
  make_option("--source", type="character", default=NULL, 
            help="source for downloading the reference genome fasta file", metavar="character"),
  make_option("--version", type="character", default=NULL, 
            help="the version of the BSgenome package to be built", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

path_single_fasta_genome <- opt$path_single_fasta_genome
out_dir <- path_to_multifasta <- opt$out_dir
latin_name <- opt$latin_name
common_name <- opt$common_name
genome_build <- opt$genome_build
seed_file_name <- opt$seed_file_name
fasta_url <- opt$fasta_url
source <- opt$source
version <- opt$version


mulifasta_dir <-
    generate_multifasta(path_single_fasta_genome =path_single_fasta_genome,
                        out_dir = out_dir)
package_seed <- generate_seed_file(path_to_multifasta = path_to_multifasta,
                                   latin_name = latin_name,
                                   common_name = common_name,
                                   genome_build = genome_build,
                                   seed_file_name = seed_file_name,
                                   fasta_url = fasta_url,
                                   source = source,
                                   version = version)
forge_install_BSgenome(package_seed[2])