library( "dplyr" )
library( "DESeq2" )
library( "biomaRt" )
library( "pheatmap" )
library( "clusterProfiler" )
library( "org.Dm.eg.db" )
library( "pathview" )
library( "readr")
library( "tibble" )
library( "forcats" )
library( "ggplot2" ) 

#customized results for domain type
FBgnLookup <- read_tsv("~/git/polycomb-development/metadata/polyATranscripts_MEtype.txt", col_names=TRUE)

#-------DEseq2 pipeline------
getResults <- function(sample, reference) {
  sampledds <- DESeqDataSetFromHTSeqCount( data.frame(sample), "~/forGEO_021220/", ~ condition )
  sampledds$condition <- relevel(sampledds$condition, ref = reference)
  sampledds <- DESeq(sampledds)
  sample_res <- results(sampledds)
  sample_res$FBgn <- rownames(sample_res)
  sample_res <- as.data.frame(sample_res)
  sample_res <- merge(sample_res, FBgnLookup, by = "FBgn" )
  return(sample_res)
}

#-------scatter plot used for figure 7A--------
plotScatter <- function(sample_res) {
  scatterPlot <- ggplot(sample_res, aes(y=log2FoldChange, x=log10(baseMean), color= ME_type))+
    geom_point(size=3)+
    theme_gray()+
    scale_color_manual(values=c("magenta", "black", "green"))+
    geom_hline(yintercept = 0, linetype="dashed")+
    ylim(-10,10)+
    xlim(0.5, 6)+
    theme(legend.position="none")
  return(scatterPlot)
}

#-------boxplots for figure 7A---------
plotBox <- function(sample_res, geneName) {
  coloredBoxPlot <- ggplot(sample_res, aes(x = ME_type, y=log2FoldChange, fill= ME_type))+
    geom_boxplot(notch=TRUE, outlier.shape = NA, range = NA)+
    geom_point(data = sample_res[which(sample_res$geneName == geneName),], shape = 19, size = 12)+
    ylim(-4, 3)+
    geom_hline(yintercept = 0, linetype="dashed")+
    scale_fill_manual(values=c("magenta", "black", "green"))+
    theme_gray()+
    theme(legend.position="none")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(coloredBoxPlot)
}
#-------function for formatting ribbon plot data for figure 7B-------
ribbonRow <- function(sample_res, ME_color, stageName, RNAi) {
  quantiles <- sapply(subset(sample_res, ME_type == ME_color, select = "log2FoldChange"), quantile, na.rm=TRUE)
  MedianFC <- quantiles["50%",]
  CI95 <- 1.58*(quantiles["75%",] - quantiles["25%",])/sqrt(nrow(subset(sample_res, ME_type == ME_color)))
  newRow <- data.frame(stageName, ME_color, RNAi, MedianFC, CI95)
  names(newRow) <- c("Stage", "ME_type", "RNAi", "MedianFC", "CI95")
  return(newRow)
}


#-------bar graphs for figure 7B--------
barGraph <- function(sample_res) {
  Activequantiles <- sapply(subset(sample_res, ME_type == "Active", select = "log2FoldChange"), quantile, na.rm=TRUE)
  ActiveMedianFC <- Activequantiles["50%",]
  ActiveCI95 <- 1.58*(Activequantiles["75%",] - Activequantiles["25%",])/sqrt(nrow(subset(sample_res, ME_type == "Active")))
  Inactivequantiles <- sapply(subset(sample_res, ME_type == "Inactive", select = "log2FoldChange"), quantile, na.rm=TRUE)
  InactiveMedianFC <- Inactivequantiles["50%",]
  InactiveCI95 <- 1.58*(Inactivequantiles["75%",] - Inactivequantiles["25%",])/sqrt(nrow(subset(sample_res, ME_type == "Inactive")))
  PcGquantiles <- sapply(subset(sample_res, ME_type == "PcG", select = "log2FoldChange"), quantile, na.rm=TRUE)
  PcGMedianFC <- PcGquantiles["50%",]
  PcGCI95 <- 1.58*(PcGquantiles["75%",] - PcGquantiles["25%",])/sqrt(nrow(subset(sample_res, ME_type == "PcG")))
  barData <- data.frame(c("Active", "Inactive", "PcG"), c(ActiveMedianFC, InactiveMedianFC, PcGMedianFC), c(ActiveCI95, InactiveCI95, PcGCI95))
  names(barData) <- c("ME_type", "MedianFC", "CI95")
  print(barData)
  BarPlot <- ggplot(barData) +
    geom_bar( aes(x=ME_type, y=MedianFC, fill=ME_type), stat="identity") +
    scale_fill_manual(values = c("magenta", "black", "green"))+
    geom_errorbar( aes(x=ME_type, ymin=MedianFC - CI95, ymax=MedianFC + CI95), width=0.4, colour="black", size=1.3)+
    ylim(-3.1,1)+
    theme(legend.position="none")+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(BarPlot)
}

#----load metadata------
setwd("~/git/polycomb-development/metadata")
Ezbam <- read_tsv( "bam_ez.tsv" )
Ez <- read_tsv( "ez.tsv" )
scm <- read_tsv( "scm.tsv" )
sce <- read_tsv( "sce.tsv" )
Pc <- read_tsv( "Pc.tsv" )
pcl <- read_tsv( "pcl.tsv" )
Jarid <- read_tsv( "Jarid.tsv" )
bam_EzYO <- read_tsv( "bam_ez_YO.tsv" )
bam_YO <- read_tsv( "bam_luc_YO.tsv")
bam_WO <- read_tsv( "bam_w_WO.tsv" )
bam_Pc_WO <- read_tsv( "bam_pc_WO.tsv" )
bam_pcl_WO <- read_tsv( "bam_pcl_WO.tsv" )
bam_sce_WO <- read_tsv( "bam_sce_WO.tsv" )
bam_scm_WO <- read_tsv( "bam_scm_WO.tsv" )
bam_jarid2_WO <- read_tsv( "bam_jarid2_WO.tsv" )

#--------generate DEseq results---------
Ezbam_res <- getResults(Ezbam, "none")
Ez_res <- getResults(Ez, "luc")
scm_res <- getResults(scm, "w")
sce_res <- getResults(sce, "w")
Pc_res <- getResults(Pc, "w")
pcl_res <- getResults(pcl, "w")
Jarid_res <- getResults(Jarid, "w")
bamYO_res <- getResults(bam_YO, "none")
bamEzYO_res <- getResults(bam_EzYO, "none")
bamWO_res <- getResults(bam_WO, "none")
bamPcWO_res <- getResults(bam_Pc_WO, "none")
bampclWO_res <- getResults(bam_pcl_WO, "none")
bamsceWO_res <- getResults(bam_sce_WO, "none")
bamscmWO_res <- getResults(bam_scm_WO, "none")
bamjarid2WO_res <- getResults(bam_jarid2_WO, "none")

#---make scatter plots for figure 7A---------
plotScatter(Ezbam_res)
plotScatter(Ez_res)
plotScatter(scm_res)
plotScatter(sce_res)
plotScatter(Pc_res)
plotScatter(pcl_res)
plotScatter(Jarid_res)


#---make box plots for figure 7A---------
plotBox(sample_res = Ezbam_res, geneName = "E(z)")
plotBox(sample_res = Ez_res, geneName = "E(z)")
plotBox(sample_res = scm_res, geneName = "Scm")
plotBox(sample_res = sce_res, geneName = "Sce")
plotBox(sample_res = Pc_res, geneName = "Pc")
plotBox(sample_res = pcl_res, geneName = "Pcl")
plotBox(sample_res = Jarid_res, geneName = "Jarid2")

#---reference stats for main text-------
summary(subset(Ezbam_res, ME_type == "Active"))
summary(subset(Ezbam_res, ME_type == "Inactive"))
summary(subset(Ezbam_res, ME_type == "PcG"))
summary(subset(Ez_res, ME_type == "Active"))
summary(subset(Ez_res, ME_type == "Inactive"))
summary(subset(Ez_res, ME_type == "PcG"))


#----------generate ribbon plots for firgure 7B---------------
ThreeColorBam <- data.frame(c(1,1,1), c("Active", "Inactive", "PcG"), c("none", "none", "none"), c(0,0,0), c(0,0,0))
names(ThreeColorBam) <- c("Stage", "ME_type", "RNAi", "MedianFC", "CI95")

ribbonPlot <- data.frame(ribbonRow(sample_res = bamWO_res, ME_color = "Active", stageName = 3 , RNAi = "none"))
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamWO_res, ME_color = "Inactive", stageName = 3 , RNAi = "none"))    
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamWO_res, ME_color = "PcG", stageName = 3 , RNAi = "none"))
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamYO_res, ME_color = "Active", stageName = 2 , RNAi = "none"))   
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamYO_res, ME_color = "Inactive", stageName = 2 , RNAi = "none"))    
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamYO_res, ME_color = "PcG", stageName = 2 , RNAi = "none"))
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamEzYO_res, ME_color = "Active", stageName = 2 , RNAi = "Ez"))   
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamEzYO_res, ME_color = "Inactive", stageName = 2 , RNAi = "Ez"))    
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = bamEzYO_res, ME_color = "PcG", stageName = 2 , RNAi = "Ez"))
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = Ezbam_res, ME_color = "Active", stageName = 1 , RNAi = "Ez"))   
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = Ezbam_res, ME_color = "Inactive", stageName = 1 , RNAi = "Ez"))    
ribbonPlot <- rbind(ribbonPlot, ribbonRow(sample_res = Ezbam_res, ME_color = "PcG", stageName = 1 , RNAi = "Ez"))
ribbonPlot <- rbind(ribbonPlot, ThreeColorBam)

ggplot(subset(ribbonPlot, RNAi == "none"), aes(x=Stage, y=MedianFC, color = ME_type))+
  geom_line()+
  scale_color_manual(values = c("magenta", "black", "green"))+
  geom_ribbon(aes(ymin=(MedianFC - CI95), ymax=(MedianFC + CI95), x=Stage, fill = ME_type), color=NA, alpha = 0.2)+
  scale_fill_manual(values = c("magenta", "black", "green"))+
  xlim(1,3)+
  ylim(-3.1,1)+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

ggplot(subset(ribbonPlot, RNAi == "Ez"), aes(x=Stage, y=MedianFC, color = ME_type))+
  geom_line()+
  scale_color_manual(values = c("magenta", "black", "green"))+
  geom_ribbon(aes(ymin=(MedianFC - CI95), ymax=(MedianFC + CI95), x=Stage, fill = ME_type), color=NA, alpha = 0.2)+
  scale_fill_manual(values = c("magenta", "black", "green"))+
  xlim(1,3)+
  ylim(-3.1,1)+
  theme(legend.position="none")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

#---generate bar graphs for figure 7B----
barGraph(Ezbam_res)
barGraph(bamYO_res)
barGraph(bamEzYO_res)
barGraph(bamWO_res)
barGraph(bampclWO_res)
barGraph(bamjarid2WO_res)
barGraph(bamscmWO_res)
barGraph(bamPcWO_res)
barGraph(bamsceWO_res)



