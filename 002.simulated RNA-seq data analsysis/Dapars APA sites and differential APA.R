library("ggplot2")
library("patchwork")
## long UTR as input
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\Dapars_long)")

f <- dir(".", "DaPars.out_All_Prediction_Results.txt.reformatted$")
names(f) <- gsub(".DaPars.out_All_Prediction_Results.txt.reformatted", "", f, perl = TRUE)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )

null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by.x = "Gene_ID", by.y = "gene", all.x = TRUE)
  cp$diff_dist_proximal[cp$utr_type == "short"] <- {ifelse(cp$strand.x[cp$utr_type == "short"] == "+",
                                                           cp$Predicted_Proximal_APA[cp$utr_type == "short"] - cp$end.y[cp$utr_type == "short"],
                                                           cp$end.y[cp$utr_type == "short"] - cp$Predicted_Proximal_APA[cp$utr_type == "short"])}
                                                            
  write.table(cp, file = paste0(.y, "diff.predicted.pAPA-real.pAPA.txt"),
              sep = "\t",
              quote = FALSE, row.names = FALSE)
  ggplot(cp, aes(x= diff_dist_proximal)) + geom_histogram(bins =100) + 
    xlab("Predicted - simulated pAPA (bp)") +
    ggtitle(gsub("01.polyester.(\\d+X).long.CPsites", "\\1", .y)) +
    theme(plot.title = element_text(hjust = 0.5))
  
}, f, names(f), SIMPLIFY = FALSE)

pdf(file = "Given.long.predict.proximalAPA.diff.pdf", height = 12, width = 12)
wrap_plots(null, nrow =4)
dev.off()


## short UTR or single UTR as input
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\DaPars_short.single)")

f <- dir(".", "DaPars.out_All_Prediction_Results.txt.reformatted$")
names(f) <- gsub(".DaPars.out_All_Prediction_Results.txt.reformatted", "", f, perl = TRUE)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )


null <- mapply(function(.x, .y){
  .x = f[1]
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by.x ="Gene_ID", by.y = "gene", all.x = TRUE)
  cp$diff_dist_distal[cp$utr_type == "long"] <- {ifelse(cp$strand.x[cp$utr_type == "long"] == "+",
                                                           cp$Predicted_Distal_APA[cp$utr_type == "long"] - cp$end.y[cp$utr_type == "long"],
                                                           cp$end.y[cp$utr_type == "long"] - cp$Predicted_Distal_APA[cp$utr_type == "long"])}
  
  write.table(cp, file = paste0(.y, "diff.predicted.dAPA-real.dAPA.txt"),
              sep = "\t",
              quote = FALSE, row.names = FALSE)
  ggplot(cp, aes(x= diff_dist_distal)) + geom_histogram(bins =100) + 
    xlab("Predicted - simulated dAPA (bp)") +
    ggtitle(gsub("01.polyester.(\\d+X).short.single.CPsites", "\\1", .y)) +
    theme(plot.title = element_text(hjust = 0.5))
}, f, names(f), SIMPLIFY = FALSE)

pdf(file = "Given.short.predict.distalAPA.diff.pdf", height = 12, width = 12)
wrap_plots(null, nrow =4)
dev.off()




##
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\Dapars_long)")

f <- dir(".", ".DaPars.out_All_Prediction_Results.txt.reformatted$")
names(f) <- gsub(".DaPars.out_All_Prediction_Results.txt.reformatted", "", f)

utr_ends <- read.delim("../005.bi-isoform.utrs.expression.PDUI.txt",
                       header = TRUE)
utr_ends <- utr_ends[, c(2,13,14,15)]


null <- mapply(function(.x, .y){
  #.x= f[1]
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by.x = "Gene_ID", by.y = "gene_id", all.x = TRUE)
  
  # write.table(cp, file = paste0(.y, ".plus.PDUI.txt"), 
  #             sep = "\t",
  #             quote = FALSE, row.names = FALSE)
  #table(cp$PASS.x, cp$PASS.y)
  # cor(x = cp$dPDUI.x, y = cp$dPDUI.y, use = "complete.obs")
  # ggplot(cp, aes(x = PDUI_Group_diff, y = dPDUI)) + geom_point(alpha = 0.3, shape = 1) +
  #          xlab("Predicted dPDUI") +
  #          ylab ("Simulated dPDUI") +
  #          ggtitle(gsub("01.polyester.(\\d+X).long.diff.APA",
  #                                            "\\1", .y)) +
  #          theme(plot.title = element_text(hjust = 0.5))
  table(cp$Pass_Filter, cp$PASS)
}, f, names(f), SIMPLIFY = FALSE)
pdf("Given.long.correlation between calculated and theoretic dPDUI.pdf", width = 10, height = 6)
wrap_plots(null, nrow =2)
dev.off()



## APA discovery rate
short <- c(250, 329, 499, 522, 524, 531, 526, 528)/566 *100
long <- c(252, 328, 494, 520, 523, 526,523,525)/566 *100
depth <- c(20, 30, 40, 60, 80, 100, 120, 140)

adj_discovery_rate <- data.frame(recovery=c(short, long), 
                             depth = c(depth, depth),
                             input = rep(c("short", "long"), each =8))

short <- c(247, 329, 499, 522, 524, 531, 526, 528)/566 *100
long <- c(252, 328, 498, 523, 525, 530,525,528)/566 *100
depth <- c(20, 30, 40, 60, 80, 100, 120, 140)

noadj_discovery_rate <- data.frame(recovery=c(short, long), 
                                 depth = c(depth, depth),
                                 input = rep(c("short", "long"), each =8))
discovery_rate <- rbind(adj_discovery_rate, noadj_discovery_rate)

discovery_rate$adjust <- rep(c("YES", "NO"), each = 16)
discovery_rate$method <- "InPAS"

short <- c(0, 0, 0, 0, 0, 0, 0, 0)
long <- c(250, 269, 395, 521, 524, 529, 526, 527)/566 *100
depth <- c(20, 30, 40, 60, 80, 100, 120, 140)
dapars <- data.frame(recovery = c(short, long),
                     depth = c(depth, depth),
                     input = rep(c("short", "long"), each =8),
                     adjust = "NO",
                     method = "DaPars")
discovery_rate <- rbind(discovery_rate, dapars)

pdf(file = "../APA.recovery.rate.versus.depth.pdf", height = 3, width = 4)
ggplot(discovery_rate, aes(x= depth, y = recovery, 
                           color = input, 
                           shape = adjust,
                           linetype = method)) +
  geom_point(size = 1.5) +
  geom_line() + 
  xlab("Mean depth") +
  ylab("Recovered APAs (%)") +
  labs(shape = "Proximal CPs adjustment", 
       color = "Provided 3' UTR",
       linetype = "Method") +theme_light()+
  theme(legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, hjust= 1, vjust = 0.5)) + 
  scale_x_continuous(breaks = c(20, 30, 40, 60, 80, 100, 120, 140))  +
  scale_y_continuous(breaks = c(0, 20, 30, 40, 50, 60, 70, 80, 90, 100))
  
dev.off()  



### proximal APA length comparison.
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long)")
inpas_f <- dir(".", ".predicted.pAPA-real.pAPA.txt$", full.names = TRUE)
names(inpas_f) <- gsub(".long.CPsitesdiff.predicted.pAPA-real.pAPA.txt", "", inpas_f)

in_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\Dapars_long)"
dapars_f <- dir(in_dir, "simulated_readsdiff.predicted.pAPA-real.pAPA.txt", full.names = TRUE)
names(dapars_f) <- gsub("^0", "", gsub(".+/(.+)simulated_readsdiff.predicted.pAPA-real.pAPA.txt", "\\1", dapars_f))

pdf(file = "Dapars.InPAS.pAPA-true pAPA.pdf", width =16,
    height = 16)
par(mfrow =c(4,4))
null <- lapply(seq_along(inpas_f), function(.x){
  #.x =4
  inpas <- read.delim(inpas_f[.x], header = TRUE, as.is = TRUE)[, c(25, 28, 29)]
  inpas <- inpas[inpas$utr_type == "short", ]
  colnames(inpas) <- c("transcript", "utr_type", "dist_diff_inpas")
  dapars <- read.delim(dapars_f[.x], header = TRUE, as.is = TRUE)[, c(22,25, 26)]
  dapars <- dapars[dapars$utr_type == "short", ]
  colnames(dapars) <- c("transcript", "utr_type", "dist_diff_dapars")
  
  dist_diff <- merge(inpas, dapars, by = "transcript")
  plot(dist_diff[, c(3,5)], main = names(dapars_f)[.x],
       xlab = "InPAS predicted - true APA (bp)",
       ylab = "DaPars predicted - true APA (bp)")
  plot(dist_diff[, c(3,5)], xlim = c(-500, 500), ylim = c(-500, 500), 
       main = names(dapars_f)[.x],
       xlab = "InPAS predicted - true APA (bp)",
       ylab = "DaPars predicted - true APA (bp)")
  
})
dev.off()



### proximal APA length comparison.
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")
inpas_f <- dir(".", ".predicted.pAPA-real.pAPA.txt$", full.names = TRUE)
names(inpas_f) <- gsub(".long.CPsitesdiff.predicted.pAPA-real.pAPA.txt", "", inpas_f)

in_dir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\Dapars_long)"
dapars_f <- dir(in_dir, "simulated_readsdiff.predicted.pAPA-real.pAPA.txt", full.names = TRUE)
names(dapars_f) <- gsub("^0", "", gsub(".+/(.+)simulated_readsdiff.predicted.pAPA-real.pAPA.txt", "\\1", dapars_f))

pdf(file = "Dapars.InPAS.adj.pAPA-true pAPA.pdf", width =16,
    height = 16)
par(mfrow =c(4,4))
null <- lapply(seq_along(inpas_f), function(.x){
  #.x =4
  inpas <- read.delim(inpas_f[.x], header = TRUE, as.is = TRUE)[, c(25, 28, 29)]
  inpas <- inpas[inpas$utr_type == "short", ]
  colnames(inpas) <- c("transcript", "utr_type", "dist_diff_inpas")
  dapars <- read.delim(dapars_f[.x], header = TRUE, as.is = TRUE)[, c(22,25, 26)]
  dapars <- dapars[dapars$utr_type == "short", ]
  colnames(dapars) <- c("transcript", "utr_type", "dist_diff_dapars")
  
  dist_diff <- merge(inpas, dapars, by = "transcript")
  plot(dist_diff[, c(3,5)], main = names(dapars_f)[.x],
       xlab = "InPAS NBC predicted - true APA (bp)",
       ylab = "DaPars predicted - true APA (bp)")
  plot(dist_diff[, c(3,5)], xlim = c(-500, 500), ylim = c(-500, 500), 
       main = names(dapars_f)[.x],
       xlab = "InPAS NBC predicted - true APA (bp)",
       ylab = "DaPars predicted - true APA (bp)")
  
})
dev.off()



### proximal APA length comparison.
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")
inpas_f_NBC <- dir(".", ".predicted.pAPA-real.pAPA.txt$", full.names = TRUE)
names(inpas_f_NBC) <- gsub(".long.CPsitesdiff.predicted.pAPA-real.pAPA.txt", "", inpas_f_NBC)


indir <- r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long)"
inpas_f <- dir(indir, ".predicted.pAPA-real.pAPA.txt$", full.names = TRUE)
names(inpas_f) <- gsub(".+?\\.([0-9]+X)", "\\1", 
                       gsub(".+/(.+?).long.CPsitesdiff.predicted.pAPA-real.pAPA.txt", "\\1", inpas_f))
names(inpas_f) <- gsub("^0", "", names(inpas_f))

pdf(file = "InPAS.NBC.versus.InPAS.pAPA-true pAPA.pdf", width =16,
    height = 16)
par(mfrow =c(4,4))
null <- lapply(seq_along(inpas_f), function(.x){
  inpas <- read.delim(inpas_f[.x], header = TRUE, as.is = TRUE)[, c(25, 28, 29)]
  inpas <- inpas[inpas$utr_type == "short", ]
  colnames(inpas) <- c("transcript", "utr_type", "dist_diff_inpas")
  inpas_nbc <- read.delim(inpas_f_NBC[.x], header = TRUE, as.is = TRUE)[, c(25,28,29)]
  inpas_nbc <- inpas_nbc[inpas_nbc$utr_type == "short", ]
  colnames(inpas_nbc) <- c("transcript", "utr_type", "dist_diff_inpas_nbc")
  
  dist_diff <- merge(inpas, inpas_nbc, by = "transcript")
  plot(dist_diff[, c(3,5)], main = names(inpas_f)[.x],
       xlab = "InPAS predicted - true APA (bp)",
       ylab = "InPAS NBC predicted - true APA (bp)")
  plot(dist_diff[, c(3,5)], xlim = c(-500, 500), ylim = c(-500, 500), 
       main = names(inpas_f)[.x],
       xlab = "InPAS predicted - true APA (bp)",
       ylab = "InPAS NBC predicted - true APA (bp)")
  abline(h = 0, col = "red")
  abline(v = 0, col = "blue")
})
dev.off()

