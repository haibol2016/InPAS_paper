library("ggplot2")

dir_1 <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\01.polyester.simulation.out)"
tx_label <- read.delim(file.path(dir_1, "004.final.tx.gene.with.label.txt"), header = FALSE)
colnames(tx_label) <- c("transcript", "gene", "type")

utrs <- read.delim(file.path(dir_1, "valid.utrs.txt"), header = FALSE)
colnames(utrs) <- c("chr", "start", "end", "width", "strand",
                    "exon_rank", "transcript", "utr3_type",
                    "gene", "symbol")
tx_label <- merge(tx_label, utrs[, -9], by = "transcript", all.x = TRUE)

is_long <- FALSE
if (is_long)
{
  cur_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\01.polyester.simulation.out\00.long.single.final.CPsites/nonpolished)"
  setwd(cur_dir)
  fs <- dir(".", ".long.single.non.polished.CPsites.txt$")
  names(fs) <- gsub(".long.single.non.polished.CPsites.txt", "",
                    fs, perl = TRUE)
  out <- "../00.annotated.long.single.CPsites.out"
} else {
  cur_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\01.polyester.simulation.out\00.short.single.final.CPsites/nonpolished)"
  setwd(cur_dir)
  fs <- dir(".", ".short.single.non.polished.CPsites.txt$")
  names(fs) <- gsub(".short.single.non.polished.CPsites.txt", "",
                    fs, perl = TRUE)
  out <- "../00.annotated.short.single.CPsites.out"
}
if (!dir.exists(out)){
  dir.create(out)
}

null <- mapply(function(f, nf){
  dat <- read.delim(f, header = TRUE)
  dat <- merge(dat, tx_label, by = "gene", all.x =TRUE)
  dat <- dat[, -c(7:13, 26, 29:32)]

  false_pos <- sum(dat$type =="single" & !is.na(dat$Predicted_Proximal_APA))
  true_pos <- sum(dat$type =="long" & !is.na(dat$Predicted_Proximal_APA))
  false_neg <- sum(dat$type =="long") -  true_pos
  true_neg <- sum(dat$type =="single") - false_pos

  confusion.matrix <- data.frame(predicted_TRUE = c(true_pos, false_pos),
             predicted_FALSE = c(false_neg, true_neg))
  rownames(confusion.matrix) <- c("actual_pos", "actual_neg")
  
  write.table(confusion.matrix, file = file.path(out, 
                                                 paste0(nf, ".counfusion.matrix.txt")), 
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  dat_distal_cp <- dat[dat$type =="single" | dat$type =="long", ]
  dat_distal_cp$diff <- ifelse(dat_distal_cp$strand.x == "+", 
                               dat_distal_cp$preliminary_distal_APA - dat_distal_cp$end.y,
                               dat_distal_cp$start.y - dat_distal_cp$preliminary_distal_APA)
  distal_cp_diff <- data.frame(diff = dat_distal_cp$diff, site = "distal")
  
  dat_2isoform_proximal_cp <- dat[dat$type !="single" & !is.na(dat$Predicted_Proximal_APA), ]
  dat_2isoform_proximal_cp <- split(dat_2isoform_proximal_cp, 
                                    f = dat_2isoform_proximal_cp$type)
  
  dat_2isoform_proximal_cp[[2]]$diff <- ifelse(dat_2isoform_proximal_cp[[2]]$strand.x == "+", 
                                               dat_2isoform_proximal_cp[[2]]$Predicted_Proximal_APA - dat_2isoform_proximal_cp[[2]]$end.y,
                                               dat_2isoform_proximal_cp[[2]]$start.y - dat_2isoform_proximal_cp[[2]]$Predicted_Proximal_APA)
  proximal_cp_diff <- data.frame(diff = dat_2isoform_proximal_cp[[2]]$diff, site = "proximal")
  plot_dat <- rbind(distal_cp_diff, proximal_cp_diff)
  plot_dat <- plot_dat[abs(plot_dat$diff) <= 2000, ]
  plot_dat$site <- factor(plot_dat$site, level = c("distal", "proximal"))
  
  pdf(file.path(out, paste0(nf, ".predicted.sites-known.sites.pdf")))
  p <- ggplot(plot_dat, aes(x = diff, fill = diff)) +
    geom_histogram(bins = 40) +
    xlim(-2000, 2000) +
    facet_wrap(.~site, scales = "free")
  print(p)
  dev.off()
  
  dat$diff <- NA
  for (i in 1:nrow(dat))
  {
    if (dat[i,]$type =="single" || dat[i, ]$type =="long"){
      dat[i,]$diff <- 
      ifelse(dat[i,]$strand.x == "+", 
             dat[i,]$preliminary_distal_APA - dat[i,]$end.y,
             dat[i,]$start.y - dat[i,]$preliminary_distal_APA)
    } else if (dat[i,]$type == "short") {
      ifelse(dat[i,]$strand.x == "+", 
             dat[i,]$Predicted_Proximal_APA - dat[i,]$end.y,
             dat[i,]$start.y - dat[i,]$Predicted_Proximal_APA)
    }
  }
  write.table(dat, file = file.path(out, f), sep = "\t",
              quote = FALSE, row.names = FALSE)
}, fs, names(fs), SIMPLIFY = FALSE)

