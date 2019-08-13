#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
library(ggplot2)
library(reshape)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(gplots)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(openxlsx)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Setup environment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set working directory to script location
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
# Update these if necessary. Paths are relative to the script location.
experiment.name <- "PC580"
protein.groups.path <- "../MaxQuant Results/proteinGroups.txt"
phosphosites.path <- "../MaxQuant Results/Phospho (STY)Sites.txt"
experimental.design.path <- "./ExperimentalDesign.csv"
phosphosites.literature.path <- "../../Phosphorylation_site_dataset.txt"
kinase.substrate.path <- "../../Kinase_Substrate_Dataset.txt"
regulatory.sites.path <- "../../Regulatory_sites.txt"
disease.associated.sites.path <- "../../Disease-associated_sites.txt"
species <- "mouse"
use.corrected = T

if (use.corrected)
{
  col.pattern <- "Reporter.intensity.corrected.*."
} else
{
  col.pattern <- "Reporter.intensity.\\d.*."
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#++++++++++++++++++++++++++++++++++++++
# Helper for pairwise correlation plot
#++++++++++++++++++++++++++++++++++++++
panel.cor <- function(x, y){
  points(x, y, pch=19, cex=0.5, col= rgb(red=0,green=0,blue=0,alpha=0.5))
  abline(a=0, b=1, col='green', lwd=1.5) # reference line for 1:1 correlation
  
  abline(a=-1, b=1, col='green', lwd=1)
  abline(a=1, b=1, col='green', lwd=1)
  
  abline(lm(y ~ x), col='red', lwd=1.5) # linear regression line for actual correlation
  r <- round(cor(x, y, use="pairwise"), digits=3)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt, cex = 1.0)
}

#++++++++++++++++++++++++++++++++++++++
# Pairwise correlation plot
#++++++++++++++++++++++++++++++++++++++
plotReporterIonPairwiseCorrelation <- function(data, design)
{
  createDir("figures")
  filename <- paste("figures/", experiment.name, "_pairwise_correlations_", design$Enriched[1], sep="")
  pdf(paste(filename, ".pdf", sep=""), width=7, height=7, compress=F)
  
  labels <- paste("Log2(intensity)\nLabel ", design$Label.Number, "\n", design$TMT.Label, "\n", design$Condition, " rep ", design$Replicate, "\n", design$Enriched, sep="")
  p <- pairs(data, labels=labels, cex.labels=0.9, lower.panel = NULL, upper.panel = panel.cor, main = paste("Pair-wise correlation of TMT reporter ion intensities -", experiment.name, design$Enriched[1]))
  #legend("bottomleft", legend=c("x=y reference line", "linear regression"), col=c("green", "red"), lty=c(1,1))
  
  dev.off()
}

#++++++++++++++++++++++++++++++++++++++
# Boxplots of reporter ion distributions
#++++++++++++++++++++++++++++++++++++++
plotReporterIonDistributions <- function(data, design)
{
  p <- (ggplot(data=data, aes(x=label, y=value, group=label, color=Condition))
        + geom_boxplot()
        + ylab("Log2(intensity)")
        + xlab("Label")
        + ggtitle(paste("Distributions of reporter ion intensities -", experiment.name, data$Enriched[1]))
        + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  )
  
  filename <- paste(experiment.name, "Boxplots", design$Enriched[1], sep="_")
  savefig(p, filename, 7, 7)
}

#++++++++++++++++++++++++++++++++++++++
# Principal Component Analysis (PCA)
#++++++++++++++++++++++++++++++++++++++
plotPCA <- function(data, design)
{
  data.t <- as.data.frame(t(data))
  pca <- prcomp(data.t)
  eigen <- pca$sdev^2
  variance <- (eigen/sum(eigen)) * 100
  variance <- format(round(variance, 2), nsmall=2) # show 2 digits after decimal
  
  plotting.data <- cbind(data.frame(pca$x), design)
  y.range <- max(plotting.data$PC2) - min(plotting.data$PC2)
  min.x <- min(plotting.data$PC1)
  max.x <- max(plotting.data$PC1)
  x.range <- max.x - min.x
  
  p <- (ggplot(plotting.data, aes(x=PC1, y=PC2, shape=as.factor(Replicate), color=Condition, label=rownames(plotting.data))) 
        + geom_point(size=3)
        + geom_text(size=3, nudge_y = y.range/20)
        
        + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
        
        + xlab(paste("PC1 (", variance[1], "% explained variance)", sep=""))
        + ylab(paste("PC2 (", variance[2], "% explained variance)", sep=""))
        + ggtitle(paste("Principal Component Analysis (PCA) -", experiment.name, design$Enriched[1]))
  )
  
  filename <- paste(experiment.name, "PCA", design$Enriched[1], sep="_")
  savefig(p, filename, 7, 4)
}

#++++++++++++++++++++++++++++++++++++++
# Histogram of fold change distribution
#++++++++++++++++++++++++++++++++++++++
plotFoldChangeDistribution <- function(data, design)
{
  for (c1 in levels(design$Condition))
  {
    abundance.c1 <- rowMeans(as.matrix(data[,design$Condition==c1]))
    for (c2 in levels(design$Condition))
    {
      if (c1 < c2)
      {
        abundance.c2 <- rowMeans(as.matrix(data[,design$Condition==c2]))
        ratios <- data.frame(ratio = abundance.c1-abundance.c2)
        p <- (ggplot(data=ratios, aes(x=ratio))
              + geom_histogram(bins=60)
              + scale_x_continuous(breaks=floor(min(ratios$ratio)):ceiling(max(ratios$ratio)))
              + xlab(paste("Average log2 ratio (", c1, " ", design$Enriched[1], " / ", c2, " ", design$Enriched[1], ")", sep=""))
              + ggtitle(paste("Distribution of log2 ratios -", experiment.name, design$Enriched[1]))
        )
        
        filename <- paste(experiment.name, "FoldChange_distribution", design$Enriched[1], c1, "vs", c2, sep="_")
        savefig(p, filename, 7, 7)
      }
    }
  }
}

#++++++++++++++++++++++++++++++++++++++
# Heatmaps
#++++++++++++++++++++++++++++++++++++++
plotHeatmap.protein <- function(data, 
                                all.data, 
                                design, 
                                bio.threshold=1, 
                                bio.method="log2.fc", # options are "log2.fc", "variance"
                                stat.threshold=0.05, 
                                stat.method="none", # options are "none", "p", "q"
                                stats=NULL,
                                cluster.cols=F,
                                show.rownames=T,
                                center=T,
                                scale=F
)
{
  # biological significance
  if (bio.method == "abs(log2(fc))")
  {
    changing <- get.changing(data, design, bio.threshold)
  }
  else if (bio.method == "variance")
  {
    changing <- rowVariances(data) >= bio.threshold
  }
  
  # statistical significance
  if (stat.method == "none")
  {
    sig.changing <- changing
  }
  else
  {
    sig.changing <- changing & stats <= stat.threshold
  }
  
  
  labels <- all.data$Labels[sig.changing]
  signficant.changing <- data[sig.changing, ]
  signficant.changing.all <- all.data[sig.changing, ]
  rownames(design) <- colnames(data)
  annotation_col <- select(design, Replicate, Condition)
  annotation_col$Replicate <- as.factor(annotation_col$Replicate)
  
  if (sum(sig.changing) == 0){
    warning(paste("Nothing passes thresholds. Heat map not created.", design$Enriched[1], bio.method, bio.threshold, stat.method, stat.threshold))
  } else
  {
    output_dir <- "heatmaps"
    createDir(output_dir)
    filename <- paste(experiment.name, "protein", bio.method, bio.threshold, stat.method, stat.threshold, sep="_")
    prefix <- paste(output_dir, filename, sep="/")
    write.table(cbind(signficant.changing.all, signficant.changing), file=paste(prefix,".txt", sep=""), sep="\t", row.names = F)
    output.protein.excel(signficant.changing.all, design, prefix)
    
    if (stat.method == "none")
    {
      selection.criteria <- paste("selection criteria:", bio.method, ">=", bio.threshold)
    } else 
    {
      selection.criteria <- paste("selection criteria:", bio.method, ">=", bio.threshold, stat.method, "<=", stat.threshold)
    }
    
    if (center & scale)
    {
      scale.name = "z-scored\nlog2(intensity)"
    } else if (center)
    {
      scale.name = "mean-centered\nlog2(intensity)"
    } else
    {
      scale.name = "log2(intensity)"
    }
    
    quant <- t(scale(t(as.matrix(signficant.changing)), center=center, scale=scale))
    rownames(quant) <- labels
    
    n.con <- nlevels(annotation_col$Condition)
    n.rep <- nlevels(annotation_col$Replicate)
    
    replicate.colors <- structure(brewer.pal(n.rep, "Set2")[1:n.rep], names=levels(annotation_col$Replicate))
    condition.colors <- structure(brewer.pal(n.con, "Set1")[1:n.con], names=levels(annotation_col$Condition))
    
    heatmap.quant.col.annot <- HeatmapAnnotation(annotation_col,
                                                 col = list(Condition = condition.colors,
                                                            Replicate = replicate.colors),
                                                 show_annotation_name = T,
                                                 annotation_name_gp = gpar(fontface = "bold", fontsize = 8),
                                                 annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface="bold"),
                                                                                labels_gp = gpar(fontsize = 8)))
    
    heatmap.quant <- Heatmap(quant,
                             # labels
                             column_title = paste(experiment.name, "protein", "\n", selection.criteria),
                             column_title_gp = gpar(fontface = "bold", fontsize=10),
                             show_row_names = show.rownames,
                             row_names_gp = gpar(fontsize = 5),
                             column_names_gp = gpar(fontsize = 8),
                             name = scale.name,
                             # colors
                             col = colorRamp2(seq(-max(abs(quant)), max(abs(quant)), length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                             # legends
                             heatmap_legend_param = list(color_bar = "continuous",
                                                         title_gp = gpar(fontsize = 8, fontface="bold"),
                                                         labels_gp = gpar(fontsize = 8)),
                             # clustering
                             cluster_columns = cluster.cols & sum(sig.changing)>1,
                             # annotations
                             top_annotation = heatmap.quant.col.annot,
                             # size
                             width=ncol(quant)*0.4,
                             clustering_distance_rows = "pearson"
    )
    
    filename <- paste(experiment.name, "HeatMap", "protein", bio.method, bio.threshold, stat.method, stat.threshold, sep="_")
    
    height <- min(7, 2.5 + (0.1*sum(sig.changing)))
    width <- 1.5 + (0.4*nrow(design))
    if (cluster.cols)
    {
      filename <- paste(filename, "clustering", sep="_")  
    }
    if (show.rownames)
    {
      filename <- paste(filename, "labeled", sep="_") 
      height <- 2.5 + (0.1*sum(sig.changing))
      width <- 2.5 + (0.4*ncol(quant))
    }
    
    savefig(heatmap.quant, filename, width, height)
  }
}


#++++++++++++++++++++++++++++++++++++++
# Heatmap for phospho
#++++++++++++++++++++++++++++++++++++++
plotHeatmap.phospho <- function(data, 
                                all.data, 
                                design, 
                                bio.threshold=1, 
                                bio.method="log2.fc", # options are "abs(log2(fc))", "variance"
                                stat.threshold=0.05, 
                                stat.method="none", # options are "none", "p", "q"
                                stats=NULL,
                                cluster.cols=F,
                                show.rownames=T,
                                center=T,
                                scale=F,
                                plot.proteins=T,
                                plot.kinase.subtrates=T
)
{
  # biological significance
  if (bio.method == "abs(log2(fc))")
  {
    changing <- get.changing(data, design, bio.threshold)
  }
  else if (bio.method == "variance")
  {
    changing <- rowVariances(data) >= bio.threshold
  }
  
  # statistical significance
  if (stat.method == "none")
  {
    sig.changing <- changing
  }
  else
  {
    sig.changing <- changing & stats <= stat.threshold
  }
  
  
  labels <- all.data$Labels[sig.changing]
  signficant.changing <- data[sig.changing, ]
  signficant.changing.all <- all.data[sig.changing, ]
  rownames(design) <- colnames(data)
  
  annotation_col <- select(design, Replicate, Condition)
  annotation_col$Replicate <- as.factor(annotation_col$Replicate)
  
  annotation_row <- select(all.data[sig.changing, ], Amino.acid, is.known.site)
  annotation_row$is.known.site <- as.factor(annotation_row$is.known.site)
  colnames(annotation_row) <- c("Amino acid", "Known site")
  
  # statistical significance
  if (stat.method == "none" || stat.method == "q")
  {
    protein.annotation_row <- select(signficant.changing.all, Protein.anova.q.value)
    protein.stat.label <- paste("Protein q-value <=", stat.threshold)
  }
  else
  {
    protein.annotation_row <- select(signficant.changing.all, Protein.anova.p.value)
    protein.stat.label <- paste("Protein p-value <=", stat.threshold)
  }
  colnames(protein.annotation_row) <- protein.stat.label
  
  
  if (sum(sig.changing) == 0){
    warning(paste("Nothing passes thresholds. Heat map not created.", design$Enriched[1], bio.method, bio.threshold, stat.method, stat.threshold))
  } else
  {
    output_dir <- "heatmaps"
    createDir(output_dir)
    filename <- paste(experiment.name, "phospho", bio.method, bio.threshold, stat.method, stat.threshold, sep="_")
    prefix <- paste(output_dir, filename, sep="/")
    write.table(signficant.changing.all, file=paste(prefix,".txt", sep=""), sep="\t", row.names = F)
    output.phospho.excel(signficant.changing.all, design, prefix)
    
    quant.protein <- subset(signficant.changing.all, select = grep(" protein ", colnames(signficant.changing.all)))
    
    if (stat.method == "none")
    {
      selection.criteria <- paste("selection criteria:", bio.method, ">=", bio.threshold)
    } else 
    {
      selection.criteria <- paste("selection criteria:", bio.method, ">=", bio.threshold, stat.method, "<=", stat.threshold)
    }
    
    if (center & scale)
    {
      scale.name = "z-scored\nlog2(intensity)"
    } else if (center)
    {
      scale.name = "mean-centered\nlog2(intensity)"
    } else
    {
      scale.name = "log2(intensity)"
    }
    
    quant <- t(scale(t(as.matrix(signficant.changing)), center=center, scale=scale))
    rownames(quant) <- labels
    
    n.con <- nlevels(annotation_col$Condition)
    n.rep <- nlevels(annotation_col$Replicate)
    
    replicate.colors <- structure(brewer.pal(n.rep, "Set2")[1:n.rep], names=levels(annotation_col$Replicate))
    condition.colors <- structure(brewer.pal(n.con, "Set1")[1:n.con], names=levels(annotation_col$Condition))
    
    protein.stat.colors <- list()
    protein.stat.colors[[protein.stat.label]] = c("TRUE" ="#80b1d3", "FALSE" = "#DDDDDD")
    
    heatmap.quant.col.annot <- HeatmapAnnotation(annotation_col,
                                                 col = list(Condition = condition.colors,
                                                            Replicate = replicate.colors),
                                                 show_annotation_name = T,
                                                 annotation_name_gp = gpar(fontface = "bold", fontsize = 8),
                                                 annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface="bold"),
                                                                                labels_gp = gpar(fontsize = 8)))
    
    heatmap.protein.quant.col.annot <- HeatmapAnnotation(annotation_col,
                                                         col = list(Condition = condition.colors,
                                                                    Replicate = replicate.colors),
                                                         show_annotation_name = F,
                                                         show_legend = F)
    
    
    
    heatmap.protein.quant.row.annot <- rowAnnotation(protein.annotation_row <= stat.threshold,
                                                     show_annotation_name = T,
                                                     annotation_name_gp = gpar(fontsize = 8),
                                                     annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface="bold"),
                                                                                    labels_gp = gpar(fontsize = 8)),
                                                     na_col = "#DDDDDD",
                                                     col = protein.stat.colors)
    
    heatmap.quant.row.annot <- rowAnnotation(annotation_row,
                                             col = list("Known site" = c("High-throughput" ="#BBBBBB", 
                                                                         "Low-throughput" = "#984ea3",
                                                                         "Novel" = "gold1"),
                                                        "Amino acid"=c("S" = "#DDDDDD", 
                                                                       "T" = "#80b1d3", 
                                                                       "Y" = "#fdb462")
                                             ),
                                             annotation_name_gp = gpar(fontsize = 8),
                                             annotation_legend_param = list(title_gp = gpar(fontsize = 8, fontface="bold"),
                                                                            labels_gp = gpar(fontsize = 8)),
                                             show_annotation_name = T)
    
    width.protein.quant <- 0
    if (plot.proteins)
    {
      width.protein.quant <- ncol(quant.protein)*0.4
      
      heatmap.protein.quant <- Heatmap(t(scale(t(as.matrix(quant.protein)), center=center, scale=scale)),
                                       cluster_columns=F,
                                       cluster_rows=F,
                                       column_title = paste(experiment.name, "protein"),
                                       column_title_gp = gpar(fontface = "bold", fontsize=10),
                                       column_names_gp = gpar(fontsize = 8),
                                       show_row_names=F,
                                       top_annotation = heatmap.protein.quant.col.annot,
                                       width=width.protein.quant,
                                       # colors
                                       col = colorRamp2(seq(-max(abs(quant)), max(abs(quant)), length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                                       na_col = "#BBBBBB",
                                       # legends
                                       show_heatmap_legend = F
      )
    }
    
    width.phospho.quant <- ncol(quant)*0.4
    
    heatmap.phospho.quant <- Heatmap(quant,
                                     # labels
                                     column_title = paste(experiment.name, "phospho", "\n", selection.criteria),
                                     column_title_gp = gpar(fontface = "bold", fontsize=10),
                                     show_row_names = show.rownames,
                                     row_names_gp = gpar(fontsize = 5),
                                     name = scale.name,
                                     column_names_gp = gpar(fontsize = 8),
                                     # colors
                                     col = colorRamp2(seq(-max(abs(quant)), max(abs(quant)), length=256), rev(colorRampPalette(brewer.pal(10, "RdBu"))(256))),
                                     # legends
                                     heatmap_legend_param = list(color_bar = "continuous",
                                                                 legend_direction = "horizontal",
                                                                 legend_width = unit(5, "cm"),
                                                                 title_position = "topcenter",
                                                                 title_gp = gpar(fontsize = 8, fontface="bold"),
                                                                 labels_gp = gpar(fontsize = 8)),
                                     # clustering
                                     cluster_columns = cluster.cols & sum(sig.changing)>1,
                                     # annotations
                                     top_annotation = heatmap.quant.col.annot,
                                     # size
                                     width=width.phospho.quant,
                                     clustering_distance_rows = "pearson"
    )
    
    ht.list <- heatmap.phospho.quant + heatmap.quant.row.annot
    
    kinases <- create.kinase.substrate.matrix(signficant.changing.all$kinases)
    
    width.kinases <- 0
    if (!is.null(kinases) && plot.kinase.subtrates)
    {
      width.kinases <- ncol(kinases)*0.15
      heatmap.kinase <- Heatmap(kinases,
                                name = "Kinase substrate",
                                col = c("Yes"="#fb8072", "No"="#DDDDDD"),
                                column_names_gp = gpar(fontsize = 7),
                                heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface="bold"),
                                                            labels_gp = gpar(fontsize = 8)),
                                width=width.kinases,
                                show_column_dend=F
      )
      ht.list <- ht.list + heatmap.kinase
    }
    
    if (plot.proteins)
    {
      ht.list <- ht.list + heatmap.protein.quant + heatmap.protein.quant.row.annot
    }
    
    fig <- make_layout(ht.list, heatmap_legend_side = "bottom")
    
    filename <- paste(experiment.name, "HeatMap", "phospho", bio.method, bio.threshold, stat.method, stat.threshold, sep="_")
    
    total.width <- width.protein.quant + width.phospho.quant + width.kinases
    height <- min(7, 3.5 + (0.1*sum(sig.changing)))
    width <- 2.5 + total.width
    if (cluster.cols)
    {
      filename <- paste(filename, "clustering", sep="_")  
    }
    if (show.rownames)
    {
      filename <- paste(filename, "labeled", sep="_") 
      height <- 3.5 + (0.1*sum(sig.changing))
      width <- 3.5 + total.width
    }
    
    savefig(fig, filename, width, height)
  }
}

#++++++++++++++++++++++++++++++++++++++
# Volcano plots
#++++++++++++++++++++++++++++++++++++++
plotVolcano <- function(data, 
                        all.data, 
                        design, 
                        bio.threshold=1, 
                        stat.threshold=0.05, 
                        stat.method="p", 
                        add.labels=T)
{
  is.phospho <- design$Enriched[1] == "phospho"
  # calc p-values using one-way anova
  stat.values <- anova.per.row(data, design)
  if (stat.method == "q")
  {
    # calc q-values based on p-values
    stat.values <- p.adjust(stat.values, method="fdr")
  }
  
  for (c1 in levels(design$Condition))
  {
    abundance.c1 <- rowMeans(as.matrix(data[,design$Condition==c1]))
    for (c2 in levels(design$Condition))
    {
      if (c1 > c2)
      {
        abundance.c2 <- rowMeans(as.matrix(data[,design$Condition==c2]))
        ratios <- abundance.c1-abundance.c2
        
        
        if (is.phospho)
        {
          ggdata <- data.frame(ratio=ratios, 
                               stat=-log10(stat.values), 
                               label=all.data$Labels,
                               significance="none",
                               amino.acid=all.data$Amino.acid,
                               is.known.site=as.factor(all.data$is.known.site))
        }
        else
        {
          ggdata <- data.frame(ratio=ratios, 
                               stat=-log10(stat.values), 
                               label=all.data$Labels,
                               significance="none")
        }
        
        
        levels(ggdata$significance) <- c("none", "sig fold-change", "sig statistic", "sig both")
        
        
        ggdata$significance[abs(ggdata$ratio) >= bio.threshold] = "sig fold-change"
        ggdata$significance[ggdata$stat >= -log10(stat.threshold)] = "sig statistic"
        ggdata$significance[abs(ggdata$ratio) >= bio.threshold & ggdata$stat >= -log10(stat.threshold)] = "sig both"
        
        if (!add.labels)
        {
          ggdata$label <- ""
        } 
        else
        {
          ggdata$label[ggdata$significance != "sig both"] <- ""
        }
        
        if (is.phospho)
        {
          p <- ggplot(data=ggdata, aes(x=ratio, y=stat, label=label, 
                                       color=significance, shape=amino.acid))#, 
          #size = is.known.site))
          quant.protein <- subset(all.data, select = grep(" protein ", colnames(all.data)))
          protein.fc <- rowMeans(as.matrix(quant.protein[,design$Condition==c1])) - rowMeans(as.matrix(quant.protein[,design$Condition==c2]))
          all.data[[paste("Protein Log2 Fold Change ", c1, "/", c2, sep="")]] <- protein.fc
        }
        else 
        {
          p <- ggplot(data=ggdata, aes(x=ratio, y=stat, label=label, 
                                       color=significance))
        }
        
        output_dir <- "volcano plots"
        createDir(output_dir)
        filename <- paste(experiment.name, design$Enriched[1], c1, c2, bio.threshold, stat.method, stat.threshold, sep="_")
        prefix <- paste(output_dir, filename, sep="/")
        write.table(cbind(all.data, ggdata), file=paste(prefix,".txt", sep=""), sep="\t", row.names = F)
        
        plotting.export.data <- subset(ggdata, select = c(ratio, stat, significance))
        colnames(plotting.export.data) <- c(paste(design$Enriched[1], " Log2 Fold Change ", c1, "/", c2, sep=""), 
                                            paste("-log10(ANOVA ", stat.method, "-value)", sep=""), 
                                            paste("Significance (",  design$Enriched[1], ")? (", stat.method, 
                                                  "-value <= ", stat.threshold, " and abs(FC) >= ", bio.threshold, ")", sep=""))
        if (is.phospho)
        {
          output.phospho.excel(cbind(all.data, plotting.export.data), design, prefix)
        } else {
          output.protein.excel(cbind(all.data, plotting.export.data), design, prefix)
        }
        
        
        
        min.x <- min(min(ggdata$ratio), -bio.threshold)
        max.x <- max(max(ggdata$ratio), bio.threshold)
        x.range <- max.x - min.x
        
        y.range <- max(ggdata$stat) - min(ggdata$stat)
        
        p <- (p 
              + geom_point(alpha=0.8) 
              + geom_text(nudge_y = 0.03 * y.range, size=3, color="#ff2020")
              
              + geom_vline(xintercept = -bio.threshold, linetype="dashed", color="#333333")
              + geom_vline(xintercept = bio.threshold, linetype="dashed", color="#333333")
              + geom_hline(yintercept = -log10(stat.threshold), linetype="dashed", color="#333333")
              
              + scale_color_manual(values=c("none"="#666666", "sig fold-change"="#f2c963", "sig statistic"="#f5aaaa", "sig both"="#ff2020"))
              #+ scale_size_manual(values=c("High-throughput"=2, "Low-throughput"=4, "No"=6))
              
              + scale_x_continuous(limits=c(min.x - x.range/10, max.x + x.range/10))
              
              + ylab(paste("-log10(",stat.method,"-value)", sep=""))
              + xlab(paste("log2(", c1, " / ", c2, ")", sep=""))
              + ggtitle(paste("Volcano plot -", experiment.name, design$Enriched[1], c1, "vs", c2))
        )
        filename <- paste(experiment.name, "Volcano", design$Enriched[1], c1, "vs", c2, bio.threshold, stat.method, stat.threshold, sep="_")
        
        if (add.labels)
        {
          filename <- paste(filename, "labeled", sep="_")  
        }
        
        savefig(p, filename, 9, 6.5)  
      }
    }
  }
}

savefig <- function(p, filename, width, height)
{
  createDir("figures")
  pdf(paste("figures/" ,filename, ".pdf", sep=""), width=width, height=height, compress=F)
  print(p)
  dev.off()
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Other functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
log2.intensity <- function(data)
{
  data[data==0] <- NA # Can't take the log of 0. Assume these are missing values.
  return(log2(data))
}

get.not.missing <- function(intensities, design)
{
  not.missing <- rep(F,nrow(intensities))
  for (condition in levels(design$Condition))
  {
    not.missing <- not.missing | rowSums(as.matrix(is.na(intensities[,design$Condition==condition]))) == 0
  }
  return(not.missing)
}

rowVariances <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2 , ...)/(dim(x)[2] - 1)
}

anova.per.row <- function(data, design)
{
  apply(data, 1, function(x) 
  {
    # x is a row from the data frame
    x.melted <- cbind(melt(x, variable_name="label", id.vars=c()), design)
    # export its p-value
    unlist(summary(aov(value ~ Condition, data=x.melted)))["Pr(>F)1"]
  })
}

t.test.per.row <- function(data, design, condition1, condition2)
{
  
  apply(data, 1, function(x) 
  {
    t.test(x=as.numeric(x[,design$Condition==condition1]),
           y=as.numeric(x[,design$Condition==condition2]))
  })
}

is.known.phosphosite.for.row <- function(proteins, positions)
{
  is.LT <- F
  is.HT <- F
  for (i in 1:length(proteins))
  {
    known.sites <- phosphosites.literature[.(proteins[i], as.numeric(positions[i])), nomatch=0]
    if (nrow(known.sites) > 0)
    {
      if (sum(known.sites$LT_LIT, na.rm = T) > 0)
      {
        is.LT <- T
      }
      if (sum(known.sites$MS_LIT, na.rm = T) > 0)
      {
        is.HT <- T
      }
    }
  }
  
  if (is.LT) return("Low-throughput")
  if (is.HT) return("High-throughput")
  return("Novel")
}

add.known.phosphosites <- function(phosphosites)
{
  apply(phosphosites, 1, function(x) 
  {
    # x is a row from the data frame
    is.known.phosphosite.for.row(unlist(strsplit(x[["Proteins"]], ";")), 
                                 unlist(strsplit(x[["Positions.within.proteins"]], ";")))
  })
}

create.kinase.substrate.matrix <- function(kinase.substrates)
{
  all.kinases <- c()
  for (kinases.for.row in kinase.substrates)
  {
    kinases <- unlist(strsplit(kinases.for.row, ";"))
    for (k in kinases)
    {
      all.kinases <- c(all.kinases, k)
    }
  }
  all.kinases <- as.factor(unique(all.kinases))
  
  kinase.matrix <- matrix(data="No", nrow=length(kinase.substrates), ncol=nlevels(all.kinases))
  
  if (length(all.kinases) == 0) return(NULL)
  colnames(kinase.matrix) <- 1:nlevels(all.kinases)
  
  for (i in 1:length(kinase.substrates))
  {
    kinases <- unlist(strsplit(as.character(kinase.substrates[i]),";"))
    for (k in kinases)
    {
      k.index <- match(k,levels(all.kinases))
      kinase.matrix[i, k.index] <- "Yes"
      colnames(kinase.matrix)[k.index] <- k
    }
  }
  
  return(kinase.matrix)
}

get.kinase.substrates.for.row <- function(proteins, positions)
{
  matched.kinases <- c()
  for (i in 1:length(proteins))
  {
    kinases <- kinase.substrates[.(proteins[i], as.numeric(positions[i])), nomatch=0]
    for (j in nrow(kinases))
    {
      matched.kinases <- c(matched.kinases, as.character(kinases$GENE[j]))
    }
  }
  if (length(matched.kinases) > 0)
  {
    return(paste(unique(matched.kinases), collapse=";"))
  } else 
  {
    return("")
  }
}

add.kinase.substrates <- function(phosphosites)
{
  apply(phosphosites, 1, function(x) 
  {
    # x is a row from the data frame
    get.kinase.substrates.for.row(unlist(strsplit(x[["Proteins"]], ";")), 
                                  unlist(strsplit(x[["Positions.within.proteins"]], ";")))
  })
}

get.regulatory.sites.for.row <- function(proteins, positions)
{
  matched.sites <- data.frame("ON_FUNCTION"=character(), "ON_PROCESS"=character(), 
                              "ON_PROT_INTERACT"=character(), "ON_OTHER_INTERACT"=character(), 
                              "NOTES"=character())
  
  for (i in 1:length(proteins))
  {
    sites <- regulatory.sites[.(proteins[i], as.numeric(positions[i])), nomatch=0]
    sites.subset <- subset(sites, select=c("ON_FUNCTION", "ON_PROCESS", 
                                           "ON_PROT_INTERACT", "ON_OTHER_INTERACT", "NOTES"))
    sites.subset[] <- lapply(sites.subset, as.character)
    matched.sites <- rbind(matched.sites, sites.subset)
  }
  
  if (nrow(matched.sites) == 0)
  {
    return(data.frame("ON_FUNCTION"="", "ON_PROCESS"="", 
                      "ON_PROT_INTERACT"="", "ON_OTHER_INTERACT"="", 
                      "NOTES"=""))
  }
  
  if (nrow(matched.sites) > 1)
  {
    matched.sites.collapsed <- apply(matched.sites, 2, function(x) {paste(x, collapse="|")})
    matched.sites.collapsed[which(matched.sites.collapsed == "|")] <- ""
    matched.sites[1,1]<- matched.sites.collapsed[["ON_FUNCTION"]]
    matched.sites[1,2] <- matched.sites.collapsed[["ON_PROCESS"]]
    matched.sites[1,3] <- matched.sites.collapsed[["ON_PROT_INTERACT"]]
    matched.sites[1,4] <- matched.sites.collapsed[["ON_OTHER_INTERACT"]]
    matched.sites[1,5] <- matched.sites.collapsed[["NOTES"]]
  }
  
  return(matched.sites[1,])
}

add.regulatory.sites <- function(phosphosites)
{
  matched.sites <- data.frame("ON_FUNCTION"=character(), "ON_PROCESS"=character(), 
                              "ON_PROT_INTERACT"=character(), "ON_OTHER_INTERACT"=character(), 
                              "NOTES"=character())
  apply(phosphosites, 1, function(x) 
  {
    # x is a row from the data frame
    matched.sites <<- rbind(matched.sites, get.regulatory.sites.for.row(unlist(strsplit(x[["Proteins"]], ";")), 
                                                                        unlist(strsplit(x[["Positions.within.proteins"]], ";"))))
  })
  
  return(matched.sites)
}

get.changing <- function(data, design, fold.change.threshold)
{
  # compute fold changes between all pairs of conditions to see if anything has a large effect size
  changing <- rep(F,nrow(data))
  for (c1 in levels(design$Condition))
  {
    abundance.c1 <- log2(rowMeans(as.matrix(2^data[,design$Condition==c1])))
    for (c2 in levels(design$Condition))
    {
      if (c1 > c2)
      {
        abundance.c2 <- log2(rowMeans(as.matrix(2^data[,design$Condition==c2])))
        changing <- changing | (abs(abundance.c1-abundance.c2) >= fold.change.threshold)
      }
    }
  }
  return(changing)
}

# emulates ggplot color scheme
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# create directory if it doesn't exist
createDir <- function(path)
{
  if (!file.exists(path))
  {
    dir.create(path)
  }
}

# get log2 fold change threshold
get.fold.change.threshold <- function(data, design, prob)
{
  ratios <- c()
  for (c1 in levels(design$Condition))
  {
    abundance.c1 <- log2(rowMeans(as.matrix((2^data[,design$Condition==c1]))))
    for (c2 in levels(design$Condition))
    {
      if (c1 > c2)
      {
        abundance.c2 <- log2(rowMeans(as.matrix((2^data[,design$Condition==c2]))))
        ratios <- c(ratios, abundance.c1-abundance.c2)
      }
    }
  }
  min(1, round(mean(abs(quantile(ratios, probs = c(prob, 1-prob)))),1))
}

get.variance.threshold <- function(data, prob)
{
  round(quantile(rowVariances(data), probs = prob),1)
}

create.labels <- function(data)
{
  is.phospho <- "Positions.within.proteins" %in% colnames(data)
  
  apply(data, 1, function(x) 
  {
    split.genes <- unlist(strsplit(x[["Gene.names"]], ";"))
    label <- split.genes[1]
    if (length(split.genes) > 1)
    {
      label <- paste(label, "*", sep="")
    }
    
    if(is.phospho)
    {
      mult <- ""
      if (x[["multiplicity"]] > 1)
      {
        mult <- paste(" ", x[["multiplicity"]], "ps", sep="")
      } 
      label <- paste(label, " p", x[["Amino.acid"]], unlist(strsplit(x[["Positions.within.proteins"]], ";"))[1], mult, sep="")
    }
    
    return(label)
  })
}

merge.unenriched <- function(phospho, proteins)
{
  # split protein into separate entries for each annotation
  proteins.long <- data.frame()
  for (i in 1:nrow(proteins))
  {
    protein.ids <- unlist(strsplit(as.character(proteins$Protein.IDs[i]), ";"))
    for (j in 1:length(protein.ids))
    {
      new.row <- proteins[i,]
      new.row$Protein.IDs <- protein.ids[j]
      proteins.long <- rbind(proteins.long, new.row)
    }
  }
  
  colnames(proteins.long) <- paste("Protein.",colnames(proteins.long), sep="")
  
  proteins.long <- data.table(proteins.long)
  setkey(proteins.long, Protein.Protein.IDs)
  
  
  phospho[,colnames(proteins.long)] <- NA
  
  for (i in 1:nrow(phospho))
  {
    protein.ids <- unlist(strsplit(as.character(phospho$Proteins[i]), ";"))
    most.abundant.row <- NULL
    
    for (j in 1:length(protein.ids))
    {
      protein.matches <- proteins.long[.(protein.ids[j]), nomatch=0]
      if (nrow(protein.matches) > 0)
      {
        for (k in nrow(protein.matches))
        {
          if (is.null(most.abundant.row) || protein.matches$Protein.Intensity[k] > most.abundant.row$Protein.Intensity)
          {
            most.abundant.row <- protein.matches[k,]
          }
        }
      }
    }
    
    if (!is.null(most.abundant.row))
    {
      phospho[i,colnames(proteins.long)] <- most.abundant.row
    }
  }
  return(phospho)
}

remove.unused.TMT.cols <- function(proteins, design)
{
  num.labels <- if (grep("TMT10", design$TMT.Label[1])) 10 else 6
  first.index <- min(grep(paste(col.pattern,"Global",sep=""), colnames(proteins)))
  label.indices <- 1:(num.labels)
  
  proteins[,-(first.index+label.indices[!(label.indices %in% design$Label.Number)]-1)]
}

remove.unused.TMT.phospho.cols <- function(phosphosites, design)
{
  num.labels <- if (grep("TMT10", design$TMT.Label[1])) 10 else 6
  first.index <- min(grep(paste(col.pattern,"Phospho",sep=""), colnames(phosphosites)))
  label.indices <- 1:(num.labels)
  to.remove <- first.index+(3*(label.indices[!(label.indices %in% design$Label.Number)]-1))
  to.remove <- c(to.remove, first.index+((3*(label.indices[!(label.indices %in% design$Label.Number)]-1))+1))
  to.remove <- c(to.remove, first.index+((3*(label.indices[!(label.indices %in% design$Label.Number)]-1))+2))
  
  phosphosites[,-to.remove]
}

#++++++++++++++++++++++++++++++++++++++
# Create minimal protein tables
#++++++++++++++++++++++++++++++++++++++
create.table.protein <- function(proteins, proteins.quant)
{
  columns.ID <- c("Protein.IDs", "Protein.names", "Gene.names", "Number.of.proteins", "Peptides", 
                  "Razor...unique.peptides", "Unique.peptides", "Sequence.coverage....", "Labels")
  columns.stats <- c("protein.anova.p.value", "protein.anova.q.value")
  
  table.out <- cbind(select(proteins, columns.ID), proteins.quant, select(proteins, columns.stats))
  
  return(table.out)
}

#++++++++++++++++++++++++++++++++++++++
# Create minimal phospho tables
#++++++++++++++++++++++++++++++++++++++
create.table.phospho <- function(phosphosites, phosphosites.quant)
{
  columns.ID <- c("Proteins", "Positions.within.proteins", "Protein.names", "Gene.names", "Localization.prob", 
                  "Number.of.Phospho..STY.", "Amino.acid", "Phospho..STY..Probabilities", "Labels",
                  "is.known.site", "kinases", "Protein.Protein.IDs", "Protein.Intensity")
  columns.quant.protein <-"Protein.*protein "
  columns.stats.protein <- c("Protein.anova.p.value", "Protein.anova.q.value")
  columns.stats.phospho <- c("phospho.anova.p.value", "phospho.anova.q.value")
  
  protein.quant <- subset(phosphosites, select = grep(columns.quant.protein, colnames(phosphosites)))
  
  table.out <- cbind(select(phosphosites, columns.ID), 
                     protein.quant,
                     select(phosphosites, columns.stats.protein),
                     phosphosites.quant,
                     select(phosphosites, columns.stats.phospho))
  
  return(table.out)
}

#++++++++++++++++++++++++++++++++++++++
# Output protein excel
#++++++++++++++++++++++++++++++++++++++
output.protein.excel <- function(proteins, design, filename)
{
  colnames(proteins)[1:9] <- c("Proteins", "Protein names", "Gene names", "Number of proteins", "Peptides",
                               "Razor + unique peptides", "Unique peptides", "Sequence coverage (%)", 
                               "Figure label")
  colnames(proteins)[(10+nrow(design)):(11+nrow(design))] <- c("Protein ANOVA p-value", "Protein ANOVA q-value")
  
  colnames(proteins)[10:(10+nrow(design)-1)] <- paste("log2(intensity)", colnames(proteins)[10:(10+nrow(design)-1)])
  
  wb <- createWorkbook("UNC Proteomics Core")
  addWorksheet(wb, "proteins") 
  
  writeData(wb, sheet = 1, proteins, rowNames = F)
  
  headerStyle.IDs <- createStyle(valign="center", halign="center", fgFill = "#CCCCCC", border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
  addStyle(wb, sheet = 1, headerStyle.IDs, rows = 1, cols = 1:9)
  
  n.con <- nlevels(design$Condition)
  condition.colors <- structure(brewer.pal(n.con, "Set1")[1:n.con], names=levels(design$Condition))
  
  i <- 10
  for (condition in levels(design$Condition))
  {
    num.condition <- sum(design$Condition == condition)
    headerStyle.quant <- createStyle(valign="center", halign="center", fgFill = condition.colors[[condition]], border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
    addStyle(wb, sheet = 1, headerStyle.quant, rows = 1, cols = i:(i+num.condition-1))
    i <- i + num.condition
  }
  
  headerStyle.stats <- createStyle(valign="center", halign="center", fgFill = "#CCCCCC", border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
  addStyle(wb, sheet = 1, headerStyle.stats, rows = 1, cols = (10+nrow(design)):ncol(proteins))
  
  geneStyle <- createStyle(numFmt="TEXT")
  addStyle(wb, sheet = 1, geneStyle, rows = 2:(nrow(proteins)+1), cols = c(3,9), gridExpand=T)
  
  setColWidths(wb, sheet=1, cols=1:ncol(proteins), widths=12)
  setRowHeights(wb, sheet=1, rows=1, heights=90)
  
  saveWorkbook(wb, paste(filename,".xlsx",sep=""), overwrite = T)
}


#++++++++++++++++++++++++++++++++++++++
# Output phospho excel
#++++++++++++++++++++++++++++++++++++++
output.phospho.excel <- function(phosphosites, design, filename)
{
  phosphosites <- subset(phosphosites, select=-Protein.Intensity)
  colnames(phosphosites)[1:12] <- c("Proteins", "Positions within proteins", "Protein names", 
                                    "Gene names", "Localiziation probability", "Number of Phospho(STY)",
                                    "Amino acid", "Phospho(STY) probabilities", "Figure label", 
                                    "Known site", "Kinases", "Matched Protein ID")
  
  colnames(phosphosites)[(13+nrow(design)):(14+nrow(design))] <- c("Protein ANOVA p-value", "Protein ANOVA q-value")
  
  colnames(phosphosites)[(15+2*nrow(design)):(16+2*nrow(design))] <- c("Phospho ANOVA p-value", "Phospho ANOVA q-value")
  
  wb <- createWorkbook("UNC Proteomics Core")
  addWorksheet(wb, "phosphosites") 
  
  writeData(wb, sheet = 1, phosphosites, rowNames = F)
  
  headerStyle.IDs <- createStyle(valign="center", halign="center", fgFill = "#CCCCCC", border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
  addStyle(wb, sheet = 1, headerStyle.IDs, rows = 1, cols = 1:12)
  
  n.con <- nlevels(design$Condition)
  condition.colors <- structure(brewer.pal(n.con, "Set1")[1:n.con], names=levels(design$Condition))
  
  # Protein quant columns
  i <- 13
  for (condition in levels(design$Condition))
  {
    num.condition <- sum(design$Condition == condition)
    headerStyle.quant <- createStyle(valign="center", halign="center", fgFill = condition.colors[[condition]], border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
    addStyle(wb, sheet = 1, headerStyle.quant, rows = 1, cols = i:(i+num.condition-1))
    i <- i + num.condition
  }
  
  # Protein stats
  headerStyle.stats <- createStyle(valign="center", halign="center", fgFill = "#CCCCCC", border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
  addStyle(wb, sheet = 1, headerStyle.stats, rows = 1, cols = (13+nrow(design)):ncol(phosphosites))
  
  # Phospho quant columns
  i <- 15 + nrow(design)
  for (condition in levels(design$Condition))
  {
    num.condition <- sum(design$Condition == condition)
    headerStyle.quant <- createStyle(valign="center", halign="center", fgFill = condition.colors[[condition]], border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
    addStyle(wb, sheet = 1, headerStyle.quant, rows = 1, cols = i:(i+num.condition-1))
    
    i <- i + num.condition
  }
  
  # Phospho stats and anything else
  headerStyle.stats <- createStyle(valign="center", halign="center", fgFill = "#CCCCCC", border="bottom", borderStyle = "thin", textDecoration="bold", wrapText=T)
  addStyle(wb, sheet = 1, headerStyle.stats, rows = 1, cols = (15+2*nrow(design)):ncol(phosphosites))
  
  
  geneStyle <- createStyle(numFmt="TEXT")
  addStyle(wb, sheet = 1, geneStyle, rows = 2:(nrow(phosphosites)+1), cols = 4)
  
  setColWidths(wb, sheet=1, cols=1:ncol(phosphosites), widths=12)
  setRowHeights(wb, sheet=1, rows=1, heights=90)
  
  saveWorkbook(wb, paste(filename,".xlsx",sep=""), overwrite = T)
}





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Phosphosites from the literature
phosphosites.literature <- read.table(phosphosites.literature.path, sep="\t", header=T, comment.char="", quote = "")#, skip=3)
phosphosites.literature <- subset(phosphosites.literature, phosphosites.literature$ORGANISM==species)
phosphosites.literature$MOD_RSD <- as.character(phosphosites.literature$MOD_RSD)
# extract residue index
phosphosites.literature$Residue.Index <- as.numeric(substr(phosphosites.literature$MOD_RSD, 2, nchar(phosphosites.literature$MOD_RSD)-2))
# convert to data.table for efficient lookup
phosphosites.literature <- data.table(phosphosites.literature)
setkey(phosphosites.literature, ACC_ID, Residue.Index)

# Kinase-substrate relationships
kinase.substrates <- read.table(kinase.substrate.path, sep="\t", header=T, comment.char="", quote = "")#, skip=3)
kinase.substrates <- subset(kinase.substrates, kinase.substrates$KIN_ORGANISM==species & kinase.substrates$SUB_ORGANISM==species)
kinase.substrates$SUB_MOD_RSD <- as.character(kinase.substrates$SUB_MOD_RSD)
kinase.substrates$Residue.Index <- as.numeric(substring(kinase.substrates$SUB_MOD_RSD, 2))
kinase.substrates <- data.table(kinase.substrates)
setkey(kinase.substrates, SUB_ACC_ID, Residue.Index)

# regulatory sites
regulatory.sites <- read.table(regulatory.sites.path, sep="\t", header=T, comment.char="", quote = "")#, skip=3)
regulatory.sites <- subset(regulatory.sites, regulatory.sites$ORGANISM==species)
regulatory.sites <- regulatory.sites[grep(".*p$", regulatory.sites$MOD_RSD), ]
regulatory.sites$MOD_RSD <- as.character(regulatory.sites$MOD_RSD)
regulatory.sites$Residue.Index <- as.numeric(substring(regulatory.sites$MOD_RSD, 2, nchar(regulatory.sites$MOD_RSD)-2))
regulatory.sites <- data.table(regulatory.sites)
setkey(regulatory.sites, ACC_ID, Residue.Index)

# disease associations
disease.associated.sites <- read.table(disease.associated.sites.path, sep="\t", header=T, comment.char="", quote="")#, skip=3)
disease.associated.sites <- subset(disease.associated.sites, disease.associated.sites$ORGANISM==species)
disease.associated.sites <- disease.associated.sites[grep(".*p$", disease.associated.sites$MOD_RSD), ]
disease.associated.sites$MOD_RSD <- as.character(disease.associated.sites$MOD_RSD)
disease.associated.sites$Residue.Index <- as.numeric(substring(disease.associated.sites$MOD_RSD, 2, nchar(disease.associated.sites$MOD_RSD)-2))
disease.associated.sites <- data.table(disease.associated.sites)
setkey(disease.associated.sites, ACC_ID, Residue.Index)

# Experimental design
design <- read.table(experimental.design.path, sep=",", header=T)
design.proteins <- subset(design, design$Enriched == "protein")
design.phospho <- subset(design, design$Enriched == "phospho")

column.names.protein <- paste(design.proteins$Condition, design.proteins$Enriched, design.proteins$Replicate, sep=" ")
column.names.phospho <- paste(design.phospho$Condition, design.phospho$Enriched, design.phospho$Replicate, sep=" ")

order.protein <- order(design.proteins$Condition, design.proteins$Replicate)
order.phospho <- order(design.phospho$Condition, design.phospho$Replicate)
order.phospho.individual <- rep(1:3, times=length(order.phospho), each=1) + rep(3*(order.phospho-1), times=1, each=3) # extra credit for figuring out how this works.

design.proteins <- design.proteins[order.protein, ]
design.phospho <- design.phospho[order.phospho, ]

# Original proteinGroups
proteins <- read.table(protein.groups.path, sep="\t", header=T, comment.char="", quote="")
# Remove unused TMT labels
proteins <- remove.unused.TMT.cols(proteins, design.proteins)
proteins$Labels <- create.labels(proteins)

# Original phosphoSites
phosphosites <- read.table(phosphosites.path, sep="\t", header=T, comment.char="", quote="")
# Remove unused TMT labels
phosphosites <- remove.unused.TMT.phospho.cols(phosphosites, design.phospho)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Protein filtering
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove potential contaminants, reverses, and identified only by site
proteins <- subset(proteins, proteins$Reverse=="" & 
                     proteins$Potential.contaminant=="" &
                     proteins$Only.identified.by.site=="")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Phosphosite filtering
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove potential contaminant and reverses
phosphosites <- subset(phosphosites, phosphosites$Reverse=="" & phosphosites$Potential.contaminant=="" )
# Localization
phosphosites <- subset(phosphosites, phosphosites$Localization.prob >= 0.75)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Add known phosphosites from phosphositePlus / literature
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
phosphosites$is.known.site <- as.factor(add.known.phosphosites(phosphosites))
phosphosites$kinases <- add.kinase.substrates(phosphosites)
phosphosites <- cbind(phosphosites, add.regulatory.sites(phosphosites))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Extract intensity columns
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins.quant <- subset(proteins, select = grep(paste(col.pattern,"Global",sep=""), colnames(proteins)))

colnames(proteins.quant) <- column.names.protein
proteins.quant <- proteins.quant[, order.protein]

# Reshape so that the separate cases for 1 phosphosite, 2 phosphosites, and 3+ phosphosites are in separate rows
phosphosites.long <- reshape(phosphosites, 
                             timevar="multiplicity",
                             idvar="id", 
                             varying=grep(paste(col.pattern,"Phospho",sep=""), colnames(phosphosites)),
                             sep="___",
                             direction='long')
phosphosites.long$Labels <- create.labels(phosphosites.long)

phosphosites.quant <- subset(phosphosites.long, select = grep(paste(col.pattern,"Phospho",sep=""), colnames(phosphosites.long)))

colnames(phosphosites.quant) <- column.names.phospho

phosphosites.quant <- phosphosites.quant[, order.phospho]


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Ratio correction
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
colSums.proteins <- colSums(proteins.quant)
ratios.proteins <- colSums.proteins / max(colSums.proteins)
print(ratios.proteins)

colSums.phospho <- colSums(phosphosites.quant)
ratios.phospho <- colSums.phospho / max(colSums.phospho)
print(ratios.phospho)

proteins.quant <- sweep(proteins.quant, 2, ratios.proteins, `/`)
phosphosites.quant <- sweep(phosphosites.quant, 2, ratios.proteins, `/`)

colSums.phospho <- colSums(phosphosites.quant)
ratios.phospho <- colSums.phospho / max(colSums.phospho)
print(ratios.phospho)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Log2 transformation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins.quant <- log2.intensity(proteins.quant)
phosphosites.quant <- log2.intensity(phosphosites.quant)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write tables so far
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output_dir <- "ratioCorrected_log2"
createDir(output_dir)
prefix <- paste(output_dir, experiment.name, sep="/")
write.table(cbind(proteins, proteins.quant), file=paste(prefix,"_proteinGroups_ratioCorrected_log2.txt", sep=""), sep="\t", row.names = F)
write.table(cbind(phosphosites.long, phosphosites.quant), file=paste(prefix,"_phosphoSites_ratioCorrected_log2.txt", sep=""), sep="\t", row.names = F)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Filter missing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# should be seen in all reps of at least 1 condition
not.missing.proteins <- get.not.missing(proteins.quant, design.proteins)
proteins.quant <- proteins.quant[not.missing.proteins, ]
proteins <- proteins[not.missing.proteins, ]

not.missing.phospho <- get.not.missing(phosphosites.quant, design.phospho)
phosphosites.quant <- phosphosites.quant[not.missing.phospho, ]
phosphosites.long <- phosphosites.long[not.missing.phospho, ]


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Imputation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins.quant[is.na(proteins.quant)] <- min(phosphosites.quant, na.rm=T)
phosphosites.quant[is.na(phosphosites.quant)] <- min(phosphosites.quant, na.rm=T)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Anova and protein alignment
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins$protein.anova.p.value <- anova.per.row(proteins.quant, design.proteins)
phosphosites.long$phospho.anova.p.value <- anova.per.row(phosphosites.quant, design.phospho)

proteins$protein.anova.q.value  <- p.adjust(proteins$protein.anova.p.value, method="fdr")
phosphosites.long$phospho.anova.q.value <- p.adjust(phosphosites.long$phospho.anova.p.value, method="fdr")

proteins.quant.annotated <- cbind(Protein.IDs = proteins$Protein.IDs, 
                                  anova.p.value = proteins$protein.anova.p.value, 
                                  anova.q.value = proteins$protein.anova.q.value, 
                                  Intensity = proteins$Intensity,
                                  proteins.quant)

phosphosites.long <- merge.unenriched(phosphosites.long, proteins.quant.annotated)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write tables so far
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
output_dir <- "filtered_imputed"
createDir(output_dir)
prefix <- paste(output_dir, experiment.name, sep="/")
write.table(cbind(proteins, proteins.quant), file=paste(prefix,"_proteinGroups_filtered_imputed.txt", sep=""), sep="\t", row.names = F)
write.table(cbind(phosphosites.long, phosphosites.quant) , file=paste(prefix,"_phosphoSites_filtered_imputed.txt", sep=""), sep="\t", row.names = F)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create minimal tables for client
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins.min <- create.table.protein(proteins, proteins.quant)
phosphosites.min <- create.table.phospho(phosphosites.long, phosphosites.quant)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Correlation plots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plotReporterIonPairwiseCorrelation(proteins.quant, design.proteins)
plotReporterIonPairwiseCorrelation(phosphosites.quant, design.phospho)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Distribution of intensities
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
proteins.quant.melted <- melt(proteins.quant, variable_name="label", id.vars=c())
repeated.design <- bind_rows(replicate(nrow(proteins.quant), design.proteins, simplify=F))
proteins.quant.melted <- cbind(proteins.quant.melted, repeated.design[order(repeated.design$Condition, repeated.design$Replicate),])

phosphosites.quant.melted <- melt(phosphosites.quant, variable_name="label", id.vars=c())
repeated.design <- bind_rows(replicate(nrow(phosphosites.quant), design.phospho, simplify=F))
phosphosites.quant.melted <- cbind(phosphosites.quant.melted, repeated.design[order(repeated.design$Condition, repeated.design$Replicate),])

plotReporterIonDistributions(proteins.quant.melted, design.proteins)
plotReporterIonDistributions(phosphosites.quant.melted, design.phospho)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# PCA
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plotPCA(proteins.quant, design.proteins)
plotPCA(phosphosites.quant, design.phospho)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Distribution of log2 fold changes
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plotFoldChangeDistribution(proteins.quant, design.proteins)
plotFoldChangeDistribution(phosphosites.quant, design.phospho)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# HeatMaps
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
variance.threshold.protein = get.variance.threshold(proteins.quant, 0.95)
variance.threshold.phospho = get.variance.threshold(phosphosites.quant, 0.95)
log.fc.threshold.protein = get.fold.change.threshold(proteins.quant, design.proteins, 0.975)
log.fc.threshold.phospho = get.fold.change.threshold(phosphosites.quant, design.phospho, 0.975)
p.value.threshold = 0.05
q.value.threshold = 0.05

for (show.rownames in c(T,F))
{
  # plot highly variant - column clustering
  plotHeatmap.protein(proteins.quant, proteins.min, design.proteins, 
                      bio.threshold=variance.threshold.protein, bio.method="variance",
                      cluster.cols=T, show.rownames=show.rownames, scale=T)
  
  plotHeatmap.phospho(phosphosites.quant, phosphosites.min, design.phospho, 
                      bio.threshold=variance.threshold.phospho, bio.method="variance",
                      cluster.cols=T, show.rownames=show.rownames, scale=T, 
                      plot.proteins=show.rownames, plot.kinase.subtrates=show.rownames)
  
  # plot signficantly different between conditions - no column clustering
  # perform one-way anova to see if anything is statistically changing between any condition
  plotHeatmap.protein(proteins.quant, proteins.min, design.proteins,
                      bio.threshold=log.fc.threshold.protein, bio.method="abs(log2(fc))",
                      stat.threshold=p.value.threshold, stat.method="p", stats=proteins$protein.anova.p.value,
                      show.rownames=show.rownames)
  
  plotHeatmap.phospho(phosphosites.quant, phosphosites.min, design.phospho, 
                      bio.threshold=log.fc.threshold.phospho, bio.method="abs(log2(fc))",
                      stat.threshold=p.value.threshold, stat.method="p", stats=phosphosites.long$phospho.anova.p.value,
                      show.rownames=show.rownames)
  
  
  
  plotHeatmap.protein(proteins.quant, proteins.min, design.proteins,
                      bio.threshold=log.fc.threshold.protein, bio.method="abs(log2(fc))",
                      stat.threshold=q.value.threshold, stat.method="q", stats=proteins$protein.anova.q.value,
                      show.rownames=show.rownames)
  
  plotHeatmap.phospho(phosphosites.quant, phosphosites.min, design.phospho, 
                      bio.threshold=log.fc.threshold.phospho, bio.method="abs(log2(fc))",
                      stat.threshold=q.value.threshold, stat.method="q",stats=phosphosites.long$phospho.anova.q.value,
                      show.rownames=show.rownames)
  
  
  
  plotHeatmap.protein(proteins.quant, proteins.min, design.proteins,
                      bio.threshold=1, bio.method="abs(log2(fc))",
                      stat.threshold=p.value.threshold, stat.method="p", stats=proteins$protein.anova.p.value,
                      show.rownames=show.rownames)
  
  plotHeatmap.phospho(phosphosites.quant, phosphosites.min, design.phospho, 
                      bio.threshold=1, bio.method="abs(log2(fc))",
                      stat.threshold=p.value.threshold, stat.method="p", stats=phosphosites.long$phospho.anova.p.value,
                      show.rownames=show.rownames)
  
  
  
  plotHeatmap.protein(proteins.quant, proteins.min, design.proteins,
                      bio.threshold=1, bio.method="abs(log2(fc))",
                      stat.threshold=q.value.threshold, stat.method="q", stats=proteins$protein.anova.q.value,
                      show.rownames=show.rownames)
  
  plotHeatmap.phospho(phosphosites.quant, phosphosites.min, design.phospho, 
                      bio.threshold=1, bio.method="abs(log2(fc))",
                      stat.threshold=q.value.threshold, stat.method="q", stats=phosphosites.long$phospho.anova.q.value,
                      show.rownames=show.rownames)  
}


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Volcano plots
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plotVolcano(proteins.quant, proteins.min, design.proteins, 
            bio.threshold=log.fc.threshold.protein, stat.threshold=p.value.threshold, 
            stat.method="p", add.labels=F)

plotVolcano(phosphosites.quant, phosphosites.min, design.phospho, 
            bio.threshold=log.fc.threshold.phospho, stat.threshold=p.value.threshold, 
            stat.method="p", add.labels=F)



plotVolcano(proteins.quant, proteins.min, design.proteins, 
            bio.threshold=1, stat.threshold=p.value.threshold, 
            stat.method="p", add.labels=F)

plotVolcano(phosphosites.quant, phosphosites.min, design.phospho, 
            bio.threshold=1, stat.threshold=p.value.threshold, 
            stat.method="p", add.labels=F)



plotVolcano(proteins.quant, proteins.min, design.proteins, 
            bio.threshold=log.fc.threshold.protein, stat.threshold=p.value.threshold, 
            stat.method="q", add.labels=F)

plotVolcano(phosphosites.quant, phosphosites.min, design.phospho, 
            bio.threshold=log.fc.threshold.phospho, stat.threshold=p.value.threshold, 
            stat.method="q", add.labels=F)


plotVolcano(proteins.quant, proteins.min, design.proteins, 
            bio.threshold=1, stat.threshold=q.value.threshold, 
            stat.method="q", add.labels=F)

plotVolcano(phosphosites.quant, phosphosites.min, design.phospho, 
            bio.threshold=1, stat.threshold=q.value.threshold, 
            stat.method="q", add.labels=F)


plotVolcano(proteins.quant, proteins.min, design.proteins, 
            bio.threshold=1, stat.threshold=q.value.threshold, 
            stat.method="q")

plotVolcano(phosphosites.quant, phosphosites.min, design.phospho, 
            bio.threshold=1, stat.threshold=q.value.threshold, 
            stat.method="q")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Write stats
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
