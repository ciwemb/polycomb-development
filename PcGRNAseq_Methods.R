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

# customized results for domain type
FBgnLookup <- read_tsv("~/genomes/dm6/polyATranscripts_MEtype.txt", col_names=TRUE)
customDomains <- function(sample_dds) {
  sample_res <- results(sample_dds)
  sample_res$FBgn <- rownames(sample_res)
  sample_res <- as.data.frame(sample_res)
  sample_res <- merge(sample_res, FBgnLookup, by = "FBgn" )
  return(sample_res)
}
# plots used for figure 2
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
plotBox <- function(sample_res) {
  coloredBoxPlot <- ggplot(sample_res, aes(x = ME_type, y=log2FoldChange, fill= ME_type))+
    geom_boxplot(notch=TRUE, outlier.shape = NA)+
    ylim(-3, 3)+
    geom_hline(yintercept = 0, linetype="dashed")+
    scale_fill_manual(values=c("magenta", "black", "green"))+
    theme_gray()+
    theme(legend.position="none")
  return(coloredBoxPlot)
}

#generate DEseq results
setwd("~/rnaseq/PcGRNAi/htseq/")
Ezbam <- read_tsv( "~/rnaseq/Fastq/Bam_ovary/htseq/bam_ez.tsv" )
Ez <- read_tsv( "~/rnaseq/Fastq/Young_ovary/htseq/ez.tsv" )
scm <- read_tsv( "scm.tsv" )
sce <- read_tsv( "sce.tsv" )
Pc <- read_tsv( "Pc.tsv" )
pcl <- read_tsv( "pcl.tsv" )
Jarid <- read_tsv( "Jarid.tsv" )

Ezbamdds <- DESeqDataSetFromHTSeqCount( data.frame(Ezbam), "~/rnaseq/Fastq/Bam_ovary/htseq", ~ condition )
Ezbamdds$condition <- relevel(Ezbamdds$condition, ref = "none")
Ezbamdds <- DESeq(Ezbamdds)
Ezdds <- DESeqDataSetFromHTSeqCount( data.frame(Ez), "~/rnaseq/Fastq/Young_ovary/htseq", ~ condition )
Ezdds$condition <- relevel(Ezdds$condition, ref = "luc")
Ezdds <- DESeq(Ezdds)
scmdds <- DESeqDataSetFromHTSeqCount( data.frame(scm), "~/rnaseq/PcGRNAi/htseq/", ~ condition )
scmdds$condition <- relevel(scmdds$condition, ref = "w")
scmdds <- DESeq(scmdds)
scedds <- DESeqDataSetFromHTSeqCount( data.frame(sce), "~/rnaseq/PcGRNAi/htseq/", ~ condition )
scedds$condition <- relevel(scedds$condition, ref = "w")
scedds <- DESeq(scedds)
Pcdds <- DESeqDataSetFromHTSeqCount( data.frame(Pc), "~/rnaseq/PcGRNAi/htseq/", ~ condition )
Pcdds$condition <- relevel(Pcdds$condition, ref = "w")
Pcdds <- DESeq(Pcdds)
pcldds <- DESeqDataSetFromHTSeqCount( data.frame(pcl), "~/rnaseq/PcGRNAi/htseq/", ~ condition )
pcldds$condition <- relevel(pcldds$condition, ref = "w")
pcldds <- DESeq(pcldds)
Jariddds <- DESeqDataSetFromHTSeqCount( data.frame(Jarid), "~/rnaseq/PcGRNAi/htseq/", ~ condition )
Jariddds$condition <- relevel(Jariddds$condition, ref = "w")
Jariddds <- DESeq(Jariddds)

Ezbam_res <- customDomains(Ezbamdds)
Ez_res <- customDomains(Ezdds)
scm_res <- customDomains(scmdds)
sce_res <- customDomains(scedds)
Pc_res <- customDomains(Pcdds)
pcl_res <- customDomains(pcldds)
Jarid_res <- customDomains(Jariddds)

plotScatter(Ezbam_res)
plotScatter(Ez_res)
plotScatter(scm_res)
plotScatter(sce_res)
plotScatter(Pc_res)
plotScatter(pcl_res)
plotScatter(Jarid_res)

plotBox(Ezbam_res)
plotBox(Ez_res)
plotBox(scm_res)
plotBox(sce_res)
plotBox(Pc_res)
plotBox(pcl_res)
plotBox(Jarid_res)