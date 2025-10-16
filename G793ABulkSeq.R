# Detaching all packages
lapply(names(sessionInfo()$otherPkgs), function(pkgs)
  detach(paste0('package:', pkgs),character.only = T,unload = T,force = T))

# Loading Libraries
suppressMessages({
  library(DESeq2)
  library(tximport)
  library(ggplot2)
  library(ggrepel)
  library(BiocParallel)
  library(pheatmap)
  library(reshape2)
  library(data.table)
  library(tidyverse)
  library(limma)
  library(gtable)
  library(grid)
  library(cowplot)
  library(threadr)
})

# =====================================================================================================================================  
# =============================================== Defining Custom Functions ===========================================================
# ===================================================================================================================================== 

# A modified PCA function that assigns shapes and colors by group
New_PCA <- function (object, var1 = NULL, var2 = NULL, ntop = NULL, labs = NULL){
  
  pcaData <- plotPCA(object, intgroup=c(var1, var2), ntop=ntop, returnData=TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  
  if (leg_pos == "none"){
    
    marg <- c(0.5, 0.5, 0, 0)
    
  } else {
      
    marg <- c(0,0,0,0)
    
    }
  
  if (!(labs)){
    
      if (is.null(var2)){
      
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), label=NULL)) + 
              geom_point(size=6) +
              theme_bw() +
              theme(axis.line=element_line(size=2), 
                    axis.ticks = element_line(color='black', size=2), 
                    axis.ticks.length = unit(0.2, 'cm'), 
                    axis.text = element_text(face='bold', color ='black', size=20),
                    text = element_text(face='bold', color ='black', size=20),
                    panel.border = element_blank(),
                    plot.margin = unit(marg, "cm"),
                    legend.position='right') +
              labs(color = var1) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    } else {
        
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), shape=get(var2), label=NULL)) + 
              geom_point(size=6) +
              theme_bw() +
              theme(axis.line=element_line(size=2), 
                    axis.ticks = element_line(color='black', size=2), 
                    axis.ticks.length = unit(0.2, 'cm'), 
                    axis.text = element_text(face='bold', color ='black', size=20),
                    text = element_text(face='bold', color ='black', size=20),
                    panel.border = element_blank(),
                    plot.margin = unit(marg, "cm"),
                    legend.position='right') +
              labs(color = var1, shape = var2) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
      }
    } else if (labs) {
      
    if (is.null(var2)){
      
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), label=name)) + 
              geom_point(size=6) +
              theme_bw() +
              geom_text_repel(size = 3) +
              theme(axis.line=element_line(size=2), 
                    axis.ticks = element_line(color='black', size=2), 
                    axis.ticks.length = unit(0.2, 'cm'), 
                    axis.text = element_text(face='bold', color ='black', size=20),
                    text = element_text(face='bold', color ='black', size=20),
                    panel.border = element_blank(),
                    plot.margin = unit(marg, "cm"),
                    legend.position='right') +
              labs(color = var1) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
    } else {
        
      pca <- ggplot(pcaData, aes(PC1, PC2, color=get(var1), shape=get(var2), label=name)) + 
              geom_point(size=6) +
              theme_bw() +
              geom_text_repel(size = 3) +
              theme(axis.line=element_line(size=2), 
                    axis.ticks = element_line(color='black', size=2), 
                    axis.ticks.length = unit(0.2, 'cm'), 
                    axis.text = element_text(face='bold', color ='black', size=20),
                    text = element_text(face='bold', color ='black', size=20),
                    panel.border = element_blank(),
                    plot.margin = unit(c(0,0,0,0), "cm"),
                    legend.position='right') +
              labs(color = var1, shape = var2) +
              xlab(paste0("PC1: ",percentVar[1],"% variance")) +
              ylab(paste0("PC2: ",percentVar[2],"% variance"))
      
      }
    }
  return(pca)
  
}

# A function that creates a list of FDR and LFC cutoff genes from the DESeq results
Make_Folded_Genes <- function (j, dat, FDR, LFC){
  
  # Extract the results from the DESeq object depending on the type of comparison
  if (length(dat[[j]]) == 1) {
    res <- results(dds_gene, name = dat[[j]], alpha = FDR)
  } else if (length(dat[[j]]) == 2){
    res <- results(dds_gene, list( dat[[j]] ), alpha = FDR)
    }
  # Remove all NAs from results
  nonNA_res <- res[!is.na(res$padj),]
  pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
  # Apply a log-fold change cutoff to the results
  folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
  folded_genes <- data.frame(folded_res)
  folded_genes$Gene <- rownames(folded_genes)
  folded_genes <- folded_genes[, c(7,1:6)]
  folded_genes <- folded_genes[rev(order(folded_genes$log2FoldChange)),]
  
  return(folded_genes)
  
}

# A function that generates and saves PCA plots and the corresponding VST files
# Will perform batch correction via limma if set to TRUE
savePCA <- function (wid=8, hei=6, reso=300) {
  # Setting the plot size based on whether there is a legend or not
  if (leg_pos == 'none'){
    
    wide <- 6
    tall <- 6
    
  }else {
      
    wide <- 8
    tall <- 6
    
    }
  
  
  # Saving the PCA plot
  if (!(Batch)){ # Batch correction == FALSE
    # Generate the VST file without batch correction and save
    vsdB_gene <- varianceStabilizingTransformation(dds_gene, blind = TRUE)
    VST <- as.data.frame(assay(vsdB_gene))
    VST$Gene <- rownames(VST)
    VST <- VST[,c(ncol(VST), 1:(ncol(VST)-1))]
    write.table(VST, file.path(outDir, paste0("VSTCounts_", Experiment, "_", out, "_", cells, ".tsv")),
                quote = F, sep = '\t', row.names = F)
    
    # Generating the PCA plot
    pca <- New_PCA(vsdB_gene, var1, var2, labs = labs, ntop = ntop) + theme(legend.position = leg_pos)
    #Displaying the PCA plot
    print(pca)
    
    if (!(labs)){ # Labels == FALSE
      if (is.null(var2)){ # var2 == NULL
        # Saving PCA plot without batch correction, without variable 2, and without labels
        
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Top", ntop, 
                                                 "Genes_Legend_", leg_pos, "_NoLabels.tiff")),
             res=300, units='in', width = wide, height=tall)
        
        print(pca)
        
        dev.off()
      } 
      else { # Batch == FALSE, Labels == FALSE, but var2 =/= NULL
        # Saving PCA plot without batch correction and without labels
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Shape_", var2, 
                                                 "_Top", ntop, "Genes_Legend_", leg_pos, "_NoLabels.tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
    } 
    else { # Batch == FALSE, Labels == TRUE
      if (is.null(var2)){ # var2 == NULL
        # Saving PCA plot without batch correction and without variable 2
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Top", ntop, 
                                                 "Genes_Legend_", leg_pos, ".tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
      else { # Batch == FALSE, Labels == TRUE, var2 =/= NULL
        # Saving PCA plot without batch correction
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Shape_", var2, 
                                                 "_Top", ntop, "Genes_Legend_", leg_pos, ".tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
    }
  }
  else { # Batch == TRUE
    # Performe the batch correction using limma and save the VST file
    assay(vsdB_gene) <- limma::removeBatchEffect(assay(vsdB_gene), vsdB_gene$Batch)
    VST <- as.data.frame(assay(vsdB_gene))
    VST$Gene <- rownames(VST)
    VST <- VST[,c(ncol(VST), 1:(ncol(VST)-1))]
    write.table(VST, file.path(outDir, paste0("VSTCounts_", Experiment, "_", out, "_", cells, "_BatchCorrected.tsv")),
                quote = F, sep = '\t', row.names = F)
    
    # Generating the PCA plot
    pca <- New_PCA(vsdB_gene, var1, var2, labs = labs, ntop = ntop) + theme(legend.position = leg_pos)
    #Displaying the PCA plot
    print(pca)
    
    if (!(labs)){ # Labels == FALSE
      if (is.null(var2)){ # var2 == NULL
        # Saving PCA plot with batch correction, without variable 2, and without labels
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_BatchCorrectedBy", 
                                                 batchVar, "_Top", ntop, "Genes_Legend_",
                                                 leg_pos, "_NoLabels.tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      } 
      else { # Batch == TRUE, Labels == FALSE, but var2 =/= NULL
        # Saving PCA plot with batch correction and without labels
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Shape_", var2, 
                                                 "_BatchCorrectedBy", batchVar, "_Top", ntop, 
                                                 "Genes_Legend_", leg_pos, "_NoLabels.tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
    } 
    else { # Batch == TRUE, Labels == TRUE
      if (is.null(var2)){ # var2 == NULL
        # Saving PCA plot with batch correction and without variable 2
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_BatchCorrectedBy", 
                                                 batchVar, "_Top", ntop, "Genes_Legend_", leg_pos, ".tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
      else { # Batch == TRUE, Labels == TRUE, var2 =/= NULL
        # Saving PCA plot with batch correction and with variable 2
        tiff(filename = file.path(outDir, paste0("PCA_", Experiment, "_Color_", var1, "_Shape_", var2, 
                                                 "_BatchCorrectedBy", batchVar, "_Top", ntop, "Genes_Legend_",
                                                 leg_pos, ".tiff")),
             res=300, units='in', width=wide, height=tall)
        
        print(pca)
        
        dev.off()
      }
    }
  }
}

# A function that creates and saves a JSON file that contains the details of the
# analysis
bulkJSON <- function(){
  
  etad <- format(Sys.time(), "%a %d %b %Y")
  emit <- format(Sys.time(), "%H:%M:%S %p %Z")
  verR <- strsplit(version$version.string, " ")[[1]][3]
  verDES <- packageDescription("DESeq2")$Version
  fastqs <- list(lapply(samples$Filename, function(x) {paste(as.character(x), sep="\n")}))
  filt <- paste("genes <=", rem, "total reads")
  if (Batch) {
    corr <- paste(Batch, "(", batchVar, ")")
  } else {
    corr <- Batch
  }
  comps <- list(lapply(names(comparisons), function(x) {paste(x, sep="\n")}))
  
  metrics <- c(etad, emit, verR, verDES, 
               fastqs, form, filt, corr, 
               FDR, LFC, cells, comps, outDir)
  
  names(metrics) <- c("Date: ", "Time: ", "R Version: ", "DESeq2 Version: ",
                      "Samples: ", "Design: ", "Genes Filtered: ", "Batch Correction: ",
                      "FDR: ", "LFC: ", "Groups: ", "Comparisons: ", "Output Directory: ")
  
  write_json(metrics, 
             file = file.path(outDir, "AnalysisInfo.json"), 
             auto_unbox = T)
}

# =====================================================================================================================================  
# ================================================== Main Code Runs Below =============================================================
# ===================================================================================================================================== 

# Setting the base.path
  base.path <- "/home/mbj-lab/Desktop/Analysis_Directory"
  person <- "JP"
  Exp <- c("G793A_PLX")

# Setting metadata file directory
  metaDir <- file.path(base.path, person, paste0("Metadata.tsv"))
  
# The metadata file is read into a dataframe
  samples1 <- read.table(metaDir, header = TRUE)

  Genome <- "GRCm39_107"
  Experiment <- Exp
  Treatment <- c("Control", "PLX3397", "PLX5622")
  Genotype <- c("WT", "G793A")
  Group <- c("WT_Control", "WT_3397", "WT_5622", "G793A_Control", "G793A_3397", "G793A_5622")
  
# Filter the metadata based on the above values
  samples <- samples1[(samples1$Genome %in% Genome),]
  #samples <- samples[(samples$Treatment %in% Treatment),]
  #samples <- samples[(samples$Genotype %in% Genotype),]
  samples <- samples[(samples$Group %in% Group),]
  
# Extracts file path information for each abundance file from the metadata
  files <- file.path(samples$Location, paste0("Kallisto_", samples$Genome), samples$Filename, "abundance.h5")
  names(files) <- paste0(samples$Sample_Name)

# Verifying that all specified files exist
  if (all(file.exists(files))==FALSE) {
    for (i in 1:length(files)){
      if (file.exists(files[i])==FALSE){
        print(paste0("File on line ", i, " is ", file.exists(files[i])))
      }
      }
    stop("Please ensure the metadata file matches the samples in analysis directory")
  } else {
    print("All files exits. Moving on.")
  }


  # Import Kallisto abundance files and summarize to the gene level using TxImport
  gene_list <- read.csv(file.path("/home/mbj-lab/Sequencing/References", samples$Species[1], "Kallisto", samples$Genome[1], paste0("tx2gene_", samples$Genome[1], ".tsv")), 
                        sep='\t', header=TRUE)
  txi.gene <- tximport(files, type = "kallisto", txIn = TRUE, txOut = FALSE, tx2gene=gene_list, 
                       countsFromAbundance = 'no' )
  
  
  # Creating DESeq object from sampleTable (only relevant data) and setting the design formula
  #form <- "~Treatment"
  #form <- "~Genotype"
  form <- "~Group"
  dds_gene <- DESeqDataSetFromTximport(txi.gene, colData=samples, design= eval(parse(text=form)))

  # Setting control group
  #dds_gene$Treatment <- relevel(dds_gene$Treatment, ref="Control")
  #dds_gene$Genotype <- relevel(dds_gene$Genotype, ref="WT")
  dds_gene$Group <- relevel(dds_gene$Group, ref="WT_Control")

# Filtering out genes with less than an average of 10 reads per sample. IF U CHANGED REMEMBER TYO CHNGE BACK
  rem <- 10*length(colnames(counts(dds_gene)))
  keep <- rowSums(counts(dds_gene)) >= rem
  print(paste0("Keeping ", nrow(dds_gene[keep,]), " out of ", nrow(dds_gene), " genes"))
  dds_gene <- dds_gene[keep,]
  remove(keep)

# Perform DESeq analysis
  dds_gene <- DESeq(dds_gene)
  
# =====================================================================================================================================  
# ===================================================== DESeq Results =================================================================
# ===================================================================================================================================== 
  
  FDR <- 0.05
  LFC <- 1
  resultsNames(dds_gene)
  
  #cells param below will name folder
  cells <- "allgroups_no372"
  out <- paste0(format(Sys.time(), "%m%d%y"), "_Analysis")
  outDir <- file.path(base.path, person, Exp, "DESeq_Results", out, cells)
  
  comparisons <- list("Group_G793A_Control_vs_WT_Control",
                      "Group_WT_3397_vs_WT_Control",
                      "Group_WT_5622_vs_WT_Control",
                      "Group_G793A_3397_vs_WT_Control",
                      "Group_G793A_5622_vs_WT_Control")
  
  
  names(comparisons) <- c("G793A_Control_vs_WT_Control",
                          "WT_3397_vs_WT_Control",
                          "WT_5622_vs_WT_Control",
                          "G793A_3397_vs_WT_Control",
                          "G793A_5622_vs_WT_Control")
  
  if (!(dir.exists(outDir))){
    if (!(dir.exists(file.path(base.path, person, Exp, "DESeq_Results", out)))){
      if (!(dir.exists(file.path(base.path, person, Exp, "DESeq_Results")))){
        dir.create(file.path(base.path, person, Exp, "DESeq_Results"))
      }
      dir.create(file.path(base.path, person, Exp, "DESeq_Results", out))
      }
    dir.create(file.path(base.path, person, Exp, "DESeq_Results", out, cells))
    }
  
  for (i in 1:length(comparisons)){
    
    # Create the output directory for the comparison and gene files
    if (!(dir.exists(file.path(outDir, names(comparisons)[i], "Genes")))){
      if (!(dir.exists(file.path(outDir, names(comparisons)[i])))){
        dir.create(file.path(outDir, names(comparisons)[i]))
        }
      dir.create(file.path(outDir, names(comparisons)[i], "Genes"))
    }
    if (!(dir.exists(file.path(outDir, names(comparisons)[i], "GSEA")))){
      dir.create(file.path(outDir, names(comparisons)[i], "GSEA"))
      }
    
    # Setting the output directory for the results files
    geneOut <- file.path(outDir, names(comparisons)[i], "Genes")
    
    print(paste0("Viewing results for: ", comparisons[[i]]))
    
    # Extract the results from the DESeq object depending on the type of comparison
    if (length(comparisons[[i]]) == 1) {
      res <- results(dds_gene, name = comparisons[[i]], alpha = FDR)
    } else if (length(comparisons[[i]]) == 2){
      res <- results(dds_gene, list( comparisons[[i]] ), alpha = FDR)
      }
    
    #resLFC <- lfcShrink(dds_gene, coef = resultsNames(dds_gene)[i], type = 'apeglm')
  # Remove all NAs from results
    nonNA_res <- res[!is.na(res$padj),]
  # Apply an FDR cutoff to the results
    pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
  # Apply a log-fold change cutoff to the results
    folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
    summary(pcutoff_res) 
    summary(folded_res)
    
    #DESeq2::plotMA(resLFC, ylim = c(-10, 10))
  
  # Saving the resulsts file with all genes in the dataset 
    All_genes <- data.frame(nonNA_res)
    All_genes$Gene <- rownames(All_genes)
    All_genes <- All_genes[, c(7,1:6)]
    All_genes <- All_genes[rev(order(All_genes$log2FoldChange)),]
    # Replacing any padj=0 with a value one order of magnitude lower than the lowest padj
    # This is so the -log10(padj) will not return Inf
    lowP <- min(subset(All_genes$padj, -log10(All_genes$padj) < Inf))
    All_genes$padj[All_genes$padj == 0] <- lowP*0.1
    All_genes$log10 <- -log10(All_genes$padj)
    All_genes$signed <- (All_genes$log2FoldChange * All_genes$log10)/abs(All_genes$log2FoldChange)
    write.table(All_genes, 
                file=file.path(geneOut, paste0("AllGenes_",  Experiment, "_", names(comparisons)[i], "_", cells, ".tsv")), 
                row.names = F, 
                quote=FALSE, 
                sep='\t')
  
  # Saving the gene list for use with GSEA
    GSEA_genes <- All_genes[,c(1,ncol(All_genes))]
    GSEA_genes <- GSEA_genes[order(GSEA_genes$signed, decreasing = T),]
    write.table(GSEA_genes, 
              file=file.path(outDir, names(comparisons)[i], "GSEA", paste0("GSEA_",  Experiment, "_", names(comparisons)[i], "_", cells, ".rnk")), 
              row.names = F, col.names = F, 
              quote=FALSE, 
              sep='\t')
  
  # Saving the results file with genes that passed the FDR cutoff
    diff_genes <- data.frame(pcutoff_res)
    diff_genes$Gene <- rownames(diff_genes)
    diff_genes <- diff_genes[, c(7,1:6)]
    diff_genes <- diff_genes[rev(order(diff_genes$log2FoldChange)),]
    write.table(diff_genes, 
                file=file.path(geneOut, paste0("SigGenes_",  Experiment, "_", names(comparisons)[i], "_FDR", FDR, "_", cells, ".tsv")),
                row.names = F, 
                quote=FALSE, 
                sep='\t')

  # Saving the genes that passed the FDR and LFC cutoffs
    folded_genes <- data.frame(folded_res)
    folded_genes$Gene <- rownames(folded_genes)
    folded_genes <- folded_genes[, c(7,1:6)]
    folded_genes <- folded_genes[rev(order(folded_genes$log2FoldChange)),]
    write.table(folded_genes, 
                file=file.path(geneOut, paste0("SigFoldGenes_",  Experiment, "_", names(comparisons)[i], "_FDR", FDR, "_LFC", LFC, "_", cells, ".tsv")), 
                row.names = F, 
                quote=FALSE, 
                sep='\t')

  }
  
  
  # Writing the raw counts into a table
  raw <- as.data.frame(txi.gene[['counts']])
  raw$Gene <- rownames(raw)
  raw <- raw[,c(ncol(raw), 1:(ncol(raw)-1))]
  write.table(raw, file.path(outDir, paste0("RawCounts_", Experiment, "_", out, "_", cells, ".tsv")),
              quote = F, sep = '\t', row.names = F)
  
# Calculating TPM and saving into a table
  RPK <- txi.gene[['counts']]/txi.gene[['length']]
  sums <- colSums(RPK)/1000000
  TPM <- matrix(ncol=ncol(txi.gene[['counts']]), nrow=nrow(RPK))
  for (i in 1:ncol(RPK)){
    
    TPM[,i] <- RPK[,i]/sums[i]
    
    }
  colnames(TPM) <- colnames(RPK)
  TPM <- as.data.frame(TPM)
  TPM$Gene <- rownames(RPK)
  TPM <- TPM[,c(ncol(TPM), 1:(ncol(TPM)-1))]
  write.table(TPM, file.path(outDir, paste0("TPMCounts_", Experiment, "_", out, "_", cells, ".tsv")),
              quote = F, sep = '\t', row.names = F)
  
  
  
# =====================================================================================================================================  
# ================================================ Dataset Visualizations =============================================================
# =====================================================================================================================================

  # Normalizing the data using DESeq2's Variance Stabilizing Transformation
    vsdB_gene <- varianceStabilizingTransformation(dds_gene, blind = TRUE)
    VST <- as.data.frame(assay(vsdB_gene))
    VST$Gene <- rownames(VST)
    VST <- VST[,c(ncol(VST), 1:(ncol(VST)-1))]
    write.table(VST, file.path(outDir, paste0("VSTCounts_", Experiment, "_", out, "_", cells, ".tsv")),
              quote = F, sep = '\t', row.names = F)
  
# =====================================================================================================================================  
# ======================================================== PCA Plots ==================================================================
# ===================================================================================================================================== 
  
  var1 <- "Genotype" # Must include for PCA generation (sets color)
  var2 <- "Treatment" # Make NULL if not using an additional variable (sets shape)
  ntop <- 1000 # Number of genes to include in PCA
  labs <- F # Make FALSE or F if not including labels
  leg_pos <- "right" # Legend position. Use "none" to remove legend
  Batch <- F # Make FALSE or F in not including batch correction
  batchVar <- "Batch" # Must include if you are including batch correction
  
  savePCA(vsdB_gene)
  
  # Saving a JSON file with the analysis details
  # Requires numerous variables above this point in the code (e.g., "Batch", "batchVar", "LFC)
  # so be sure to have those variables filled in as desired/relevant
  bulkJSON()
  
# =====================================================================================================================================  
# ======================================================== Heatmaps ===================================================================
# ===================================================================================================================================== 
  
  # Creating the Heatmap output directory
  if (!(dir.exists(file.path(outDir, "Heatmaps")))){
    dir.create(file.path(outDir, "Heatmaps"))
  }
  
  
  # Setting FDR and LFC cutoffs for the heatmap               
  FDR <- FDR
  LFC <- LFC
  
  # Setting the heatmap cell height
  height <- 8
  width <- 8
  
  # Creating and saving a median-centered heatmap for the "SigFold" genes from each comparison
  
  for (j in 1:length(comparisons)){
    
    # Creating color palette for heatmap
    colfunc <- colorRampPalette(c("midnightblue", "dodgerblue4", "lightskyblue1", "#1FA187FF", "darkgreen", "darkslategray"))
    #if want all genes, not just significant use line below and comment out folded_genes lines 607 and 609
    vsdB_filt <- vsdB_gene
     
    #####
    # Importing various gene lists
    #filt_genes <- read_tsv("/home/mbj-lab/Analysis_Directory/JP/G793A_PLX/heatmaplist.tsv")
    #filt_genes <- read_tsv("/home/mbj-lab/Analysis_Directory/JP/G793A_PLX/microgliasignature.tsv")
    
    # Filtering the data by the chosen gene list and assaying
    #vsdB_filt <- vsdB_filt[(rownames(vsdB_filt) %in% filt_genes$Gene),]
    
  
    #####
    # Make the list of genes to filter the heatmap
    #folded_genes <- Make_Folded_Genes(j, comparisons, FDR, LFC)
    
    
    #default below
   # vsdB_filt <- vsdB_gene[(rownames(vsdB_gene) %in% folded_genes$Gene),]
    
    #applies filtering criteria
    vsdB_filt_assay <- assay(vsdB_filt)
    
    # Median centering the filtered data
    vsdB_filt_assay <- t(as.data.frame(scale(t(vsdB_filt_assay), scale=F, center=T)))
    
    annot <- as.data.frame(colData(dds_gene)[,c("Genotype", "Treatment")])
    rownames(annot) <- colnames(vsdB_filt_assay)
    colnames(annot) <- c("Genotype", "Treatment")
    annotation_colors <- list(Genotype=c(WT="red", G793A="darkred"),
                              Treatment=c(Control="darkblue", PLX3397="cyan", PLX5622 = "green"))
    
    
    # Generate the heatmap
    map <- pheatmap(vsdB_filt_assay,
                    annotation_col = annot,
                    annotation_colors = annotation_colors,
                    border_color = NA,
                    cellwidth = width,
                    cellheight = height,
                    cluster_cols = T,
                    clustering_distance_cols = "euclidean",
                    cluster_rows = T,
                    clustering_distance_rows = "euclidean",
                    show_rownames = T,
                    show_colnames = T,
                    fontsize = height,
                    legend = T,
                    treeheight_row = 0,
                    treeheight_col = 10,
                    #cutree_cols=2, 
                    color = colfunc(100))
    #color = viridis(20))
    
    map$gtable$grobs[[1]]$gp <- gpar(lwd = 4)
    
    png(filename = file.path(outDir, "Heatmaps", 
                             paste0("Heatmap_", Experiment, "_SigFoldGenes_", names(comparisons)[j], "_FDR", FDR, "_LFC", LFC, "_", cells, "_Centered.png")), 
        height = (height/72 * nrow(vsdB_filt_assay))+2, width = (width/72 * ncol(vsdB_filt_assay)) + 2, units ='in', res = 300)
    
    grid.newpage()
    grid.draw(map$gtable)
    dev.off()
    
  }
  
  
  # =====================================================================================================================================  
  # ===================================================== Volcano Plots =================================================================
  # ===================================================================================================================================== 
  
  # Creating the Heatmap output directory
  if (!(dir.exists(file.path(outDir, "Volcano")))){
    dir.create(file.path(outDir, "Volcano"))
  }
  
  # Setting FDR cutoff for volcano plot                
  FDR <- FDR
  LFC <- LFC
  numLabs <- 30 # Number of genes to label
  
  for (k in 1:length(comparisons)){
    # Make the list of genes to filter the heatmap
    # Extract the results from the DESeq object depending on the type of comparison
    if (length(comparisons[[k]]) == 1) {
      res <- results(dds_gene, name = comparisons[[k]], alpha = FDR)
    } else if (length(comparisons[[k]]) == 2){
      res <- results(dds_gene, list( comparisons[[k]] ), alpha = FDR)
    }
    
    # Remove all NAs from results
    nonNA_res <- res[!is.na(res$padj),]
    # Replacing any padj=0 with a value one order of magnitude lower than the lowest padj
    lowP <- min(subset(nonNA_res$padj, -log10(nonNA_res$padj) < Inf))
    nonNA_res$padj[nonNA_res$padj == 0] <- as.numeric(lowP*0.1)
    # Apply an FDR cutoff to the results
    pcutoff_res <- nonNA_res[(nonNA_res$padj<FDR),]
    # Apply a log-fold change cutoff to the results
    folded_res <- pcutoff_res[(pcutoff_res$log2FoldChange >= LFC | pcutoff_res$log2FoldChange <= -LFC),]
    folded_genes <- data.frame(folded_res)
    folded_genes$Gene <- rownames(folded_genes)
    folded_genes <- folded_genes[, c(7,1:6)]
    
    # Calculating the ideal limit values for the x and y axes
    logMin <- min(nonNA_res$log2FoldChange)
    logMax <- max(nonNA_res$log2FoldChange)
    pMax <- max(-log10(nonNA_res$padj))
    
    print(paste0("logMin: ", logMin))
    print(paste0("logMax: ", logMax))
    print(paste0("pMax: ", pMax))
    
    logRound <- 2
    pRound <- 10
    
    xMin <- (logRound*round((logMin-(logRound/2))/logRound))
    xMax <- (logRound*round((logMax+(logRound/2))/logRound))
    xDivide <- 2
    yMax <- (pRound*round((pMax+(pRound/2))/pRound))
    yDivide <- 10
    
    if (TRUE) {
      # Generating and saving the volcano plot
      mutateddf <- mutate(as.data.frame(nonNA_res), 
                          Cutoff=rownames(nonNA_res) %in% rownames(folded_res),
                          Color=ifelse(nonNA_res$log2FoldChange>0, '#AD0505', '#0C0C80'),
                          Combo=nonNA_res$log2FoldChange*-log10(nonNA_res$padj)) #Will have different colors depending on significance
      mutateddf <- within(mutateddf, Color[Cutoff == "FALSE"] <- 'gray80')
      rownames(mutateddf) <- rownames(nonNA_res)
      
      input <- cbind(gene=rownames(mutateddf), mutateddf) #convert the rownames to a column
      input <- input[order(input$Combo, decreasing = T),]
      
      geneLabels <- subset(input, input$Cutoff == "TRUE")
      geneUp <- subset(geneLabels, geneLabels$log2FoldChange > 0) # Keeping only positive LFC genes
      geneDown <- subset(geneLabels, geneLabels$log2FoldChange < 0) # Keeping only negative LFC genes
      
      # Setting the color scale based on which subsets of genes are present
      if (nrow(geneUp)== 0 && nrow(geneDown) == 0){
        cols <- c("gray80")
      } else if (nrow(geneUp) == 0 && nrow(geneDown) > 0){
        cols <- c("#0C0C80", "gray80")
      } else if (nrow(geneUp) >0 && nrow(geneDown) == 0) {
        cols <- c("#AD0505", "gray80") 
      } else {
        cols <- c("#0C0C80", "#AD0505", "gray80")
      }
      
      # Checking if there are enough positive LFC genes to label the number set
      # by numLabs. If not, using all genes in the list
      if (nrow(geneUp) >= numLabs) {
        geneUp <- head(geneUp, numLabs)
      }
      
      # Checking if there are enough negative LFC genes to label the number set
      # by numLabs. If not, using all genes in the list
      if (nrow(geneDown) >= numLabs){
        geneDown <- tail(geneDown, numLabs)
      }
      
      volc <- ggplot(input, aes(log2FoldChange, -log10(padj))) + #volcano plot with log2Foldchange versus pvalue
        theme_classic() +
        theme(axis.line = element_line(linewidth=1), 
              axis.text = element_text(color = "black", size = 10, face = "bold"),
              axis.ticks = element_line(color = "black", linewidth = 1),
              axis.title = element_text(face = 'bold'),
              plot.title = element_text(face = 'bold'),
              legend.position = "NA") +
        geom_point(data=input, aes(col=Color)) + #add points colored by significance
        scale_color_manual(values = cols) + 
        xlim(xMin, xMax) +
        ylim(0, yMax) +
        scale_x_continuous(limits = c(-8, 8), breaks = seq(-8, 8, by=xDivide), expand = c(0, 0)) +
        #to adjust y axis, change c(0, "desired value") and seq(0, "desired value", by=yDivide
        scale_y_continuous(limits = c(0, 70), breaks = seq(0, 70, by=yDivide), expand = c(0, 0)) +
        geom_hline(yintercept = -log10(FDR), linetype = "dashed", size = 1, color = "gray30") +
        geom_vline(xintercept = c(LFC, -LFC), linetype = "dashed", size = 1, color = "gray30") +
        geom_label_repel(data=rbind(geneUp, geneDown), color="white", size=3, segment.color = 'black', 
                         aes(label=gene, fill = alpha(c(Color), 0.7))) +
        scale_fill_identity()
      
      # Saving the volcano plot  
      tiff(filename = file.path(outDir, "Volcano", 
                                paste0("Volcano_", Experiment, "_", names(comparisons)[k], "_FDR", FDR, "_LFC", LFC, "_", cells, ".tiff")),
           res=300, units='in', width=7.2, height=6, compression = 'lzw')
      
      plot(volc)
      
      dev.off()
      
    }
  }


# =====================================================================================================================================  
# =================================================== Correlation Matrix ==============================================================
# ===================================================================================================================================== 
  
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(ggcorrplot)
  
  rID_assay <- assay(vsdB_gene)
  
  cormat <- round(cor(rID_assay, method='pearson'), 2)
  reorder_cormat <- function(cormat){
    # Use correlation between variable as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <- cormat[hc$order, hc$order]
  }
  
  cormat <- reorder_cormat(cormat)
  melted_cormat <- melt(cormat)
  melted_cormat$squared <- melted_cormat[,"value"]^2
  
  ggplot(data=melted_cormat, aes(x=Var1, y=Var2, fill=squared)) + 
    geom_tile() +
    labs(x=NULL, y=NULL) +
    #scale_fill_distiller(type="seq", palette=colfunc(20), guide=FALSE) + 
    guides(fill=guide_colorbar(title='R^2', ticks=FALSE)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
  
  ggsave(filename=file.path(outDir, "Correlation", paste0("CorPlot_", Experiment, "_", out, ".tiff")), 
         width=5, height=3.8, dpi=300, unit='in', compression = 'lzw')
  


  
  
  
  
