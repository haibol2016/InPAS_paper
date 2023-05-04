library("ggplot2")
library("patchwork")
## long UTR as input
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")

f <- dir(".", "long.CPsites.txt$")
names(f) <- gsub(".txt", "", f, perl = TRUE)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )

null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by = "gene", all.x = TRUE)
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


## long APA
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")

f <- dir(".", "long.diff.APA.txt$")
names(f) <- gsub("diff.APA.txt", "", f)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )

null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- cp[cp$PASS == "TRUE", ]
  cp <- merge(cp, utr_ends, by = "transcript", all.x = TRUE)
  
  write.table(cp, file = paste0(.y, "sig.diff.APA.txt"), 
              sep = "\t",
              quote = FALSE, row.names = FALSE) 
}, f, names(f), SIMPLIFY = FALSE)




## short UTR or single UTR as input
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\short.single_20_adj_pAPA)")

f <- dir(".", "short.single.CPsites.txt$")
names(f) <- gsub(".txt", "", f, perl = TRUE)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )


null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by = "gene", all.x = TRUE)
  cp$diff_dist_distal[cp$utr_type == "long"] <- {ifelse(cp$strand.x[cp$utr_type == "long"] == "+",
                                                           cp$Predicted_Distal_APA[cp$utr_type == "long"] - cp$end.y[cp$utr_type == "long"],
                                                           cp$end.y[cp$utr_type == "long"] - cp$Predicted_Distal_APA[cp$utr_type == "long"])}
  cp <- cp[cp$utr_type == "long", ]
  
  sum(!is.na(cp$diff_dist_distal) & abs(cp$diff_dist_distal) <= 50)/962
  # write.table(cp, file = paste0(.y, "diff.predicted.dAPA-real.dAPA.txt"),
  #             sep = "\t",
  #             quote = FALSE, row.names = FALSE)
  # ggplot(cp, aes(x= diff_dist_distal)) + geom_histogram(bins =100) + 
  #   xlab("Predicted - simulated dAPA (bp)") +
  #   ggtitle(gsub("01.polyester.(\\d+X).short.single.CPsites", "\\1", .y)) +
  #   theme(plot.title = element_text(hjust = 0.5))
}, f, names(f), SIMPLIFY = FALSE)

num_correct <- do.call("c", null)
num_correct <- data.frame(correct = num_correct,
                          depth = c(20, 30, 40, 60, 80, 100, 120,140))



pdf(file = "Given.short.predict.distalAPA.diff.pdf", height = 12, width = 12)
wrap_plots(null, nrow =4)
dev.off()


## single short APA
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\short.single_20_adj_pAPA)")

f <- dir(".", "short.single.diff.APA.txt$")
names(f) <- gsub("diff.APA.txt", "", f)

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )

null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- cp[cp$PASS == "TRUE", ]
  cp <- merge(cp, utr_ends, by = "transcript", all.x = TRUE)
  
  write.table(cp, file = paste0(.y, "sig.diff.APA.txt"), 
              sep = "\t",
              quote = FALSE, row.names = FALSE) 
}, f, names(f), SIMPLIFY = FALSE)


##
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")

f <- dir(".", ".long.diff.APA.txt$")
names(f) <- gsub(".txt", "", f)

utr_ends <- read.delim("../005.bi-isoform.utrs.expression.PDUI.txt",
                       header = TRUE)
utr_ends <- utr_ends[, c(2,13,14,15)]


null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by.x = "gene", by.y = "gene_id", all.x = TRUE)
  
  # write.table(cp, file = paste0(.y, ".plus.PDUI.txt"), 
  #             sep = "\t",
  #             quote = FALSE, row.names = FALSE)
  #table(cp$PASS.x, cp$PASS.y)
  # cor(x = cp$dPDUI.x, y = cp$dPDUI.y, use = "complete.obs")
  # ggplot(cp, aes(x = dPDUI.x, y = dPDUI.y)) + geom_point(alpha = 0.3, shape = 1) +
  #          xlab("Predicted dPDUI") +
  #          ylab ("Simulated dPDUI") +
  #          ggtitle(gsub("01.polyester.(\\d+X).long.diff.APA",
  #                                            "\\1", .y)) +
  #          theme(plot.title = element_text(hjust = 0.5))
  table(cp$PASS.x, cp$PASS.y)
}, f, names(f), SIMPLIFY = FALSE)
pdf("Given.long.correlation between calculated and theoretic dPDUI.pdf", width = 10, height = 6)
wrap_plots(null, nrow =2)
dev.off()


## short.single
setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\short.single_20_adj_pAPA)")

f <- dir(".", ".short.single.diff.APA.txt$")
names(f) <- gsub(".txt", "", f)

utr_ends <- read.delim("../005.bi-isoform.utrs.expression.PDUI.txt",
                       header = TRUE)
utr_ends <- utr_ends[, c(2,13,14,15)]


null <- mapply(function(.x, .y){
  cp <- read.delim(.x, header = TRUE, as.is = TRUE)
  cp <- merge(cp, utr_ends, by.x = "gene", by.y = "gene_id", all.x = TRUE)
  
  write.table(cp, file = paste0(.y, ".plus.PDUI.txt"),
              sep = "\t",
              quote = FALSE, row.names = FALSE)
  cor(x = cp$dPDUI.x, y = cp$dPDUI.y, use = "complete.obs")
  
  # 
  # ggplot(cp, aes(x = dPDUI.x, y = dPDUI.y)) + geom_point(alpha = 0.3, shape = 1) +
  #   xlab("Predicted dPDUI") +
  #   ylab ("Simulated dPDUI") +
  #   ggtitle(gsub("01.polyester.(\\d+X).short.single.diff.APA",
  #                "\\1", .y)) +
  #   theme(plot.title = element_text(hjust = 0.5))
  
  table(cp$PASS.x, cp$PASS.y)
}, f, names(f), SIMPLIFY = FALSE)
pdf("Given.short.correlation between calculated and theoretic dPDUI.pdf", width = 10, height = 6)
wrap_plots(null, nrow =2)
dev.off()







## proximal CPS accuracy give long

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long_20_adj_pAPA)")
f <- dir(".", "long.CPsitesdiff.predicted.pAPA-real.pAPA.txt$")

names(f) <- gsub("^0", "", gsub("01.polyester.+?(\\d+X).+", "\\1", f))

long_inpas_adj <- sapply(seq_along(f), function(.x){
 
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- dat[dat$utr_type == "short", ]
  sum(!is.na(dat$diff_dist_proximal) & abs(dat$diff_dist_proximal) <= 50)/ 962
})


setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\long)")
f <- dir(".", "long.CPsitesdiff.predicted.pAPA-real.pAPA.txt$")

names(f) <- gsub("^0", "", gsub("01.polyester.+?(\\d+X).+", "\\1", f))

long_inpas <- sapply(seq_along(f), function(.x){
  
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- dat[dat$utr_type == "short", ]
  sum(!is.na(dat$diff_dist_proximal) & abs(dat$diff_dist_proximal) <= 50)/ 962
})

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\Dapars_long)")
f <- dir(".", "diff.predicted.pAPA-real.pAPA.txt$")

names(f) <- gsub("^0", "", gsub("simulated_readsdiff.predicted.pAPA-real.pAPA.txt", "", f))

long_dapars <- sapply(seq_along(f), function(.x){
  
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- dat[dat$utr_type == "short", ]
  sum(!is.na(dat$diff_dist_proximal) & abs(dat$diff_dist_proximal) <= 50)/ 962
})

dat <- data.frame(accuracy = c(long_inpas, long_inpas_adj, long_dapars) *100,
                  depth = rep(c(20, 30, 40, 60, 80, 100, 120,140), 3),
                  method = rep(c("InPAS", "InPAS", "DaPars"), each = 8),
                  adjust = rep(c("NO", "YES", "NO"), each = 8))
#pdf(file = "../accuracy of proximal.CPs.versus.depth.pdf", height = 3, width = 4)
p1 <- ggplot(dat, aes(x = depth, y = accuracy, color = method, linetype = adjust))+
  geom_point(size = 1.5) +
  geom_line() + 
  xlab("Mean depth") +
  ylab("Correct proximal CPS (%)") +
  labs(linetype = "Proximal CPS adj.", color = "Method") +
  theme_bw() + theme(panel.grid.major =  element_blank(),
                     panel.grid.minor =  element_blank(),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6))+
  scale_x_continuous(breaks = c(20, 30, 40, 60, 80, 100, 120,140))
p1
#dev.off()


## false positive rate

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\short.single)")
f <- dir(".", "short.single.sig.diff.APA.txt$")

names(f) <- gsub("^0", "", gsub("01.polyester.(\\d+X).short.single.sig.diff.APA.txt", "\\1", f))

single_Inpas <- sapply(seq_along(f), function(.x){
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- dat[dat$utr_type == "single", ]
  nrow(dat)/5689 * 100
})

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\short.single_20_adj_pAPA)")
f <- dir(".", "short.single.sig.diff.APA.txt$")

names(f) <- gsub("^0", "", gsub("01.polyester.(\\d+X).short.single.sig.diff.APA.txt", "\\1", f))

single_Inpas_adj <- sapply(seq_along(f), function(.x){
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- dat[dat$utr_type == "single", ]
  nrow(dat)/5689 * 100
})

setwd(r"(C:\Users\liuh\OneDrive - University Of Massachusetts Medical School\Desktop\InPAS\08022022_final\Polyester.simulation\DaPars_short.single)")

f <- dir(".", "simulated_reads.DaPars.out_All_Prediction_Results.txt.reformatted")

utr_ends <- read.delim("../004.tx4simulation.short.long.single.utr.end.txt",
                       header = FALSE)
colnames(utr_ends) <- c("gene", "transcript", "end", "strand", "utr_type" )

single_dapars <- sapply(seq_along(f), function(.x){
  dat <- read.delim(f[.x], header = TRUE, as.is = TRUE)
  dat <- merge(dat, utr_ends, by.x = "Tx_ID" , by.y = "transcript", all.x = TRUE)
  
  dat <- dat[dat$utr_type == "single", ]
  nrow(dat)/5689 *100
})



dat <- data.frame(false_rate = c(single_Inpas, single_Inpas_adj, single_dapars),
                  depth = rep(c(20, 30, 40, 60, 80, 100, 120,140), 3),
                  method = rep(c("InPAS", "InPAS", "DaPars"), each = 8),
                  adjust = rep(c("NO", "YES", "NO"), each = 8))

#pdf(file = "../False positive rate.versus.depth.pdf", height = 3, width = 4)
p2 <- ggplot(dat, aes(x = depth, y = false_rate, color = method, linetype = adjust))+
  geom_point(size = 1.5) +
  geom_line() + 
  xlab("Mean depth") +
  ylab("APA FDR (%)") +
  labs(linetype = "pCPS adjustment", color = "Method") +
  theme_bw() + theme(panel.grid.major =  element_blank(),
                     panel.grid.minor =  element_blank(),
                     legend.title = element_text(size = 8),
                     legend.text = element_text(size = 6))+
  scale_x_continuous(breaks = c(20, 30, 40, 60, 80, 100, 120,140))
p2
#dev.off()


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

short <- c(0, 0, 0, 0, 0,0, 0, 0)/566 *100
long <- c(250, 269, 395, 521, 524, 529,526,527)/566 *100
depth <- c(20, 30, 40, 60, 80, 100, 120, 140)

dapars_discovery_rate <- data.frame(recovery=c(short, long), 
                                   depth = c(depth, depth),
                                   input = rep(c("short", "long"), each =8))
dapars_discovery_rate$adjust <- rep("NO", 16)

discovery_rate <- rbind(discovery_rate, dapars_discovery_rate)
discovery_rate$method <- rep(c("InPAS", "InPAS", "DaPars"), each = 16)

p3 <- ggplot(discovery_rate, aes(x= depth, y = recovery, 
                                 color = method, 
                                 linetype = adjust,
                                 shape = input)) +
  geom_point(size = 1.5) +
  geom_line() + 
  xlab("Mean depth") +
  ylab("Diff. APA accuracy (%)") +
  labs(linetype = "pCPs adjustment", color = "Method", shape = "Provided 3' UTRs") + theme_bw() +
  theme(panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  scale_x_continuous(breaks = c(20, 30, 40, 60, 80, 100, 120,140)) +
  scale_y_continuous(breaks = c(0, 40, 50, 60, 70, 80, 90, 100))
p3

pdf(file = "../pCP accuracy, false positive rate, dAPA recovery.rate.versus.depth.pdf", height = 3, width = 8)
p <- p1+ p2 + plot_layout(guides = "collect")
wrap_plots(p, p3, guides = "auto") & theme(axis.text.x= element_text(size = 6, 
                                                                          vjust = 0.5,
                                                                          hjust = 1,
                                                                          angle = 90),
                                               axis.text.y= element_text(size = 6),
                                               axis.title=element_text(size=8)) 

dev.off()


## given short to predict long
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
} 
gg_color_hue(4)[3]


num_correct <- rbind(num_correct, data.frame(correct = 0, 
                                             depth = c(20, 30, 40,
                                                       60, 80, 100,
                                                       120, 140)))
num_correct[,1] <- num_correct[,1]*100
num_correct$method <- rep(c("InPAS", "DaPars"), each = 8)

p4 <- ggplot(num_correct, aes(x = depth, y = correct, color = method))+
  geom_point(size = 1.5, shape= 2) +
  geom_line() + 
  xlab("Mean depth") +
  ylab("dCPS accuracy (%)") +
  labs(color = "Method") +
  theme_bw() +
  theme(panel.grid.major =  element_blank(),
        panel.grid.minor =  element_blank(),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6)) +
  scale_x_continuous(breaks = c(20, 30, 40, 60, 80, 100, 120, 140)) +
  scale_y_continuous(breaks = c(0, 20, 30, 50, 60, 70, 80, 90))

  
pdf(file = "../dCP and pCP accuracy, false positive rate, dAPA recovery.rate.versus.depth.pdf",
    height = 4.5, width = 7)

layout <- "AABB
           CCDD"
wrap_plots(p1, p4, p2, p3, guides = "auto", nrow = 2,
           design = layout) + plot_annotation(tag_level = "A") & theme(axis.text.x= element_text(size = 6, 
                                                                     vjust = 0.5,
                                                                     hjust = 1,
                                                                     angle = 90),
                                           axis.text.y= element_text(size = 6),
                                           axis.title=element_text(size=8),
                                           plot.tag = element_text(size = 8, face = "bold")) 

dev.off()



## HeLa APA
library("ggplot2")
dist <- read.delim("001.InPAS.DaPars.pAPA.dist.diff.txt", 
                   header = TRUE, as.is = TRUE)

pdf("002.InPAS.vs.DaPars.dAPA2nearest.PolyAsites.dist.pdf", height = 4, width = 5.5)
ggplot(dist, aes(x= InPAS, y = DaPars, alpha = 0.05)) +
  geom_point(size =0.01) + xlim(-1000, 1000) +ylim(-1000, 1000) +
  xlab("InPAS (bp)") + ylab("DaPars (bp)")+
  geom_abline(slope =1, intercept = 0, color = "red") +theme_bw()
dev.off()