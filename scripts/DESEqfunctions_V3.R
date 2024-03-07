# Function to create breaks and colour vectors scaled to 0 (middle) and a specified upper and lower bound value
scaleColors <- function(data = input_scale, # data to use
                        maxvalue = NULL # value at which the color is fully red / blue
){
  if(is.null(maxvalue)){
    maxvalue <- floor(min(abs(min(data)), max(data)))
  }
  if(max(data) > abs(min(data))){
    if(ceiling(max(data)) == maxvalue){
      myBreaks <- c(floor(-max(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(max(data)))
    } else{
      myBreaks <- c(floor(-max(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(max(data)))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  } else {
    if(-floor(min(data)) == maxvalue){
      myBreaks <- c(floor(min(data)), seq(-maxvalue+0.2, maxvalue-0.2, 0.2),  ceiling(min(data)))
    } else{
      myBreaks <- c(floor(min(data)), seq(-maxvalue, maxvalue, 0.2),  ceiling(abs(min(data))))
    }
    paletteLength <- length(myBreaks)
    myColor <- colorRampPalette(c("blue", "white", "red"))(paletteLength)
  }
  return(list(breaks = myBreaks, color = myColor))
}

# Plot BoxPlot of highest expressed genes

highestGenes <- function(numGenes = 10, data = norm_anno, plot_title = NULL)
{
  tmp <- data[, colnames(data) %in% sample_table$ID]
  tmp <- tmp[order(rowMeans(tmp), decreasing = T), ] # order according to maximal mean expression value
  tmp <- tmp[1:numGenes, ]
  tmp <- melt(t(tmp))
  colnames(tmp) <- c("sample", "gene", "value")
  idx <- match(tmp$gene, data$GENEID)
  tmp$symbol <- as.factor(data$SYMBOL[idx])
  tmp$symbol <- factor(tmp$symbol, levels = rev(unique(tmp$symbol)))
  tmp$description <- data$DESCRIPTION[idx]
  tmp$genetype <- data$GENETYPE[idx]
  
  if(is.null(plot_title)){
    title <- paste("Expression of", numGenes, "highest expressed genes")
  } else {
    title <- plot_title
  }
  
  ggplot(tmp, aes(x = tmp$symbol, y = value)) +
    geom_quasirandom(size = 0.8) +
    geom_boxplot(aes(fill = genetype), alpha = 0.6, outlier.shape = NA) +
    scale_fill_brewer(palette = "Paired", name = "Gene Type") +
    xlab("") +
    ylab("Normalized Expression") +
    ggtitle(paste(title)) +
    theme_bw() +
    coord_flip() +
    theme(axis.text.x = element_text(size = 18, angle = 90, hjust = 1),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          plot.title = element_text(size = 18, face = "bold")) +
    guides(fill = guide_legend(override.aes = list(alpha = 1)))
}


# Function to plot heatmaps based on pheatmap
plotHeatmap <- function(input = norm_anno,
                        geneset = "all",
                        title = "",
                        keyType = "Ensembl",
                        show_rownames = FALSE,
                        cluster_cols = FALSE,
                        column_annotation = c("donor", "treatment", "condition"),
                        sample_annotation = sample_table,
                        plot_mean = FALSE,
                        plot_mean_column = "condition")
{
  if (geneset[1] != "all") {
    if (keyType == "Ensembl") {
      input <- input[input$GENEID %in% geneset, ]
    }
    else if (keyType == "Symbol") {
      input <- input[input$SYMBOL %in% geneset, ]
    }
    else {
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep = "") # collapse GeneID and SYMBOL for rownames. Changes to the gene annotation should be done here
  
  input <- input[, colnames(input) %in% sample_annotation$ID]
  
  if(plot_mean == FALSE ){
    
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[, order(sample_annotation[[plot_order]], decreasing = FALSE)]
    column_annotation <- sample_table[,c(column_annotation), drop = F]
    
  } else {
    
    if(!plot_mean_column %in% names(sample_table)){
      stop("Error: plot_mean_cloumn needs to be a column name in your column_annotation")
    }
    input$gene <- rownames(input)
    tmp <- gather(input, key = "sample", value = "expression", -gene) %>%
      merge(., sample_table[,c("ID", plot_mean_column)], by.x = "sample", by.y = "ID") %>%
      group_by_(plot_mean_column, "gene")
    tmp$expression <- as.numeric(tmp$expression)
    tmp <- summarise(tmp, mean = mean(expression))
    tmp <- as.data.frame(spread_(tmp, key = plot_mean_column, value = "mean"))
    rownames(tmp) <- tmp$gene
    tmp$gene <- NULL
    
    input_scale <- t(scale(t(tmp)))
    column_annotation <- data.frame(row.names = unique(sample_table[,plot_mean_column]),
                                    plot_mean_column = unique(sample_table[,plot_mean_column]))
    names(column_annotation)[1] <- plot_mean_column
  }
  
  pheatmap(input_scale, main = title,
           show_rownames = show_rownames,
           show_colnames = TRUE,
           cluster_cols = cluster_cols,
           fontsize = 7,
           annotation_col = column_annotation[colnames(input_scale), ,drop = F],
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]],
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]], 
           border_color = NA, use_raster = TRUE)
}



# Function to plot heatmaps based on pheatmap
plotMeanHeatmap <- function(geneset,
                        title="",
                        keyType = "Ensembl",
                        show_rownames = FALSE,
                        cluster_cols = FALSE){
  if(geneset[1] =="all"){
    input <- norm_anno
  }else{
    if(keyType == "Ensembl"){
      input <- norm_anno[norm_anno$GENEID %in% geneset,]
    } else if(keyType == "Symbol"){
      input <- norm_anno[norm_anno$SYMBOL %in% geneset,]
    } else{
      print("Wrong keyType. Choose Ensembl or Symbol!")
    }
  }
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
  input <- input[,colnames(input) %in% sample_table$ID]
  input$gene <- rownames(input)
  group_gene <- "gene"
  tmp <- melt(input)
  tmp <- merge(tmp, sample_table, by.x = "variable", by.y = "ID")
  tmp <- tmp %>%  group_by(gene, Genotype_Age) %>% summarise(mean = mean(value))
  tmp <- spread(tmp, key = Genotype_Age, value = mean)
  tmp <- as.data.frame(tmp)
  rownames(tmp) <- tmp$gene
  tmp$gene <- NULL
  input <- tmp
  input_scale <- t(scale(t(input)))
  
  pheatmap(input_scale,
           main=title,
           show_rownames=show_rownames,
           show_colnames=TRUE,
           cluster_cols = cluster_cols, 
           fontsize = 7,
           breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]], 
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
}



# Heatmap of genes of specified GO, KEGG or HALLMARK gene sets
plotGeneSetHeatmap <- function(input. = norm_anno,
                               sample_annotation. = sample_table,
                               cat,
                               term,
                               organism,
                               show_rownames =TRUE,
                               cluster_cols = FALSE,
                               column_annotation. =  plot_annotation){
  if(organism == "mouse"){
    GO <- GO_mm
    KEGG <- KEGG_mm
  } else if(organism == "human"){
    GO <- GO_hs
    KEGG <- KEGG_hs
  } else (stop("Wrong organism specified!"))
  
  xterm <- paste("^", term, "$", sep="")
  if(cat=="GO"){
    genes <- unique(GO[grep(xterm,GO$TERM),"SYMBOL"])
  }
  if(cat=="KEGG"){
    genes <- unique(KEGG[grep(xterm,KEGG$PATHWAY),"SYMBOL"])
  }
  if(cat=="HALLMARK"){
    genes <- unique(hallmark_genes[grep(xterm,hallmark_genes$ont),"gene"])
    if(organism == "mouse"){
      genes <- getLDS(attributes = c("entrezgene_id"),
                      filters = "entrezgene_id",
                      values = genes,
                      mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                      attributesL = c("mgi_symbol"),
                      martL = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                      uniqueRows=T)[,2]
    }
  }
  
  plotHeatmap(input = input.,
              sample_annotation = sample_annotation.,
              geneset = genes,
              keyType = "Symbol",
              title = paste("Heatmap of present genes annotated to: ",
                            term, sep=""),
              show_rownames = show_rownames,
              cluster_cols = cluster_cols,
              column_annotation = column_annotation.)
}

# Function to plot a PCA
plotPCA <- function(pca_input = dds_vst,
                    pca_sample_table = sample_table,
                    ntop=500,
                    xPC=1,
                    yPC=2,
                    color,
                    anno_colour,
                    shape="NULL",
                    point_size=3,
                    title="PCA",
                    label = NULL,
                    label_subset = NULL){
  
  if(!is.data.frame(pca_input)){
    vst_matrix <-as.matrix(assay(pca_input))
  }else{
    vst_matrix <- pca_input
  }
  
  if(ntop=="all"){
    pca <- prcomp(t(vst_matrix))
  }else{
    # select the ntop genes by variance
    select <- order(rowVars(vst_matrix), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(vst_matrix[select,]))
  }
  
  #calculate explained variance per PC
  explVar <- pca$sdev^2/sum(pca$sdev^2)
  # transform variance to percent
  percentVar <- round(100 * explVar[c(xPC,yPC)], digits=1)
  
  # Define data for plotting
  pcaData <- data.frame(xPC=pca$x[,xPC],
                        yPC=pca$x[,yPC],
                        color = pca_sample_table[[color]],
                        name= as.character(pca_sample_table$ID),
                        stringsAsFactors = F)
  
  #plot PCA
  if(is.factor(pcaData$color) || is.character(pcaData$color)|| is.integer(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size)
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_shape_discrete(name=shape)
    }
    
    if(anno_colour[1] == "NULL"){
      pca_plot <- pca_plot + scale_color_discrete(name=color)
    }else{
      pca_plot <- pca_plot + scale_color_manual(values=anno_colour, name=color)
    }
    
  }else if(is.numeric(pcaData$color)){
    if(shape == "NULL"){
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color)
    }else{
      pcaData$shape = pca_sample_table[[shape]]
      pca_plot <- ggplot(pcaData, aes(x = xPC, y = yPC, colour=color, shape=shape)) +
        geom_point(size =point_size) +
        scale_color_gradientn(colours = bluered(100),name=color)+
        scale_shape_discrete(name=shape)
    }
  }
  
  # adds a label to the plot. To label only specific points, put them in the arument label_subset
  if (!is.null(label) == TRUE){
    pcaData$label <- pca_sample_table[[label]]
    if(!is.null(label_subset) == TRUE){
      pcaData_labeled <- pcaData[pcaData$label %in% label_subset,]
    } else {
      pcaData_labeled <- pcaData
    }
    pca_plot <- pca_plot +
      geom_text_repel(data = pcaData_labeled, aes(label = label), nudge_x = 2, nudge_y = 2, colour = "black")
  }
  
  pca_plot <- pca_plot+
    xlab(paste0("PC ",xPC, ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC ",yPC,": ", percentVar[2], "% variance")) +
    coord_fixed()+
    theme_classic()+
    theme(aspect.ratio = 1)+
    ggtitle(title)
  
  pca_plot
}

# Function to plot heatmaps of PC loadings
plotLoadings <- function(pca_input = dds_vst, 
                         heatmap_input = norm_anno, 
                         sample_annotation = sample_table, 
                         PC, 
                         ntop, 
                         col_annotation = NULL){
  if(ntop=="all"){
    pca <- prcomp(t(assay(pca_input)))
  }else{
    select <- order(rowVars(assay(pca_input)), decreasing=TRUE)[c(1:ntop)]
    pca <- prcomp(t(assay(pca_input)[select,]))
  }
  
  Loadings <- pca$rotation[,PC]
  Loadings <- Loadings[order(Loadings, decreasing = T)]
  Loadings <- names(Loadings[c(1:20,(length(Loadings)-19):length(Loadings))])
  
  heatmap <- heatmap_input[heatmap_input$GENEID %in% Loadings,]
  rownames(heatmap) <- paste(heatmap$GENEID,": ",heatmap$SYMBOL,sep="")
  heatmap <- heatmap[,colnames(heatmap) %in% sample_annotation$ID]
  heatmap_scale <- as.matrix(t(scale(t(heatmap))))
  
  column_annotation <- sample_table[colnames(heatmap_scale),]
  if(!is.null(col_annotation)){
    column_annotation <- column_annotation[,col_annotation, drop = F]
  }
  
  # Heatmap
  pheatmap(heatmap_scale,
           main=paste("Hierarchical Clustering of top20 ",PC, " loadings in both directions",sep=""),
           show_rownames=TRUE,
           show_colnames = TRUE,
           annotation_col = column_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(heatmap_scale, 2)[["breaks"]],
           color = scaleColors(heatmap_scale, 2)[["color"]],
           cluster_cols = T,
           fontsize=6)
}


# Function for plotting multiple plots in grid
multiplot<-function(plots=plots,
                    cols=1){
  
  layout <- matrix(seq(1, cols * length(plots)/cols),
                   ncol = cols, 
                   nrow = length(plots)/cols)
  
  
  if (length(plots)==1) {
    print(plots[[1]])
  }else{
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:length(plots)){
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Function to plot the normalized or batch-corrected counts for a single gene in a box plot
plotSingleGene <-function(data=norm_anno, 
                          symbol, 
                          condition="Genotype_Age", 
                          anno_colour=col_genotype_age,
                          shape = NULL) {
  
  input<-as.data.frame(data)
  rownames(input)<- input$GENEID
  
  if(sum(input$SYMBOL == symbol) == 0){
    stop("Gene not present")
  }else{
    plots<-list()
    for (i in 1:sum(input$SYMBOL == symbol)) {
      geneCounts <- as.data.frame(t(input[input$SYMBOL == symbol, colnames(input) %in% sample_table$ID]))
      geneCounts$condition <- sample_table[[condition]]
      GENEID<-colnames(geneCounts)[i]
      colnames(geneCounts)[i]<-"y"
      
      if(!is.null(anno_colour)){
        if (is.null(shape)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_manual(values=anno_colour)
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_manual(values=anno_colour)+
            geom_beeswarm(cex = 3, na.rm=T, aes(shape=shape)) +
            scale_shape(name=legend_shape)
        }
      }else{
        if (is.null(shape_opt)){
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_brewer(palette = "Spectral")
        }else{
          geneCounts$shape <- sample_table[[shape]]
          legend_shape<-paste0(shape)
          plot<-ggplot(geneCounts, aes(x = condition, y = y, colour=condition)) +
            scale_color_brewer(palette = "Spectral")+
            geom_beeswarm(cex = 3, na.rm=T, aes(shape=shape)) +
            scale_shape(name=legend_shape)
        }
      }
      plots[[i]]<-plot+
        geom_boxplot(width=.5,alpha=0) + 
        stat_boxplot(geom ='errorbar',width=.25) +
        ylab("Normalized counts") +
        scale_y_continuous(expand=c(0.05,0.25)) +
        expand_limits(y=0) +
        labs(title=paste(symbol, GENEID, sep=": "),colour=condition)+
        theme_classic()+
        theme(plot.title = element_text(hjust=0.5))
    }
    if(sum(input$SYMBOL== symbol)>1){
      print("Selected gene symbol assigned to more than one gene (Ensembl ID)")
      multiplot(plots)
    }else{
      print("Selected gene symbol assigned to one gene (Ensembl ID)")
      multiplot(plots)
    }
  }
}

# Function to remove potential batch effects using the removeBatchEffect function from limma
limmaBatchEffectRemoval <- function(input=dds_vst,
                                    batchfactor, # name of batch effect column in sample_table
                                    batchfactor_2=NULL,
                                    modelfactor){ # name of model effect column in sample_table
  
  # rlog-transformed input
  x <- as.matrix(assay(input)) 
  
  # design matrix
  model <- model.matrix(~sample_table[,c(modelfactor)])
  
  # run batch remocal function
  if(is.numeric(sample_table[,colnames(sample_table) == batchfactor[1]])==T){
    as.data.frame(removeBatchEffect(x,
                                    covariates = sample_table[,colnames(sample_table) %in% batchfactor],
                                    design = model))
  }else{
    if(is.null(batchfactor_2)){
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      design = model))
    }else{
      as.data.frame(removeBatchEffect(x=x,
                                      batch = sample_table[,colnames(sample_table) == batchfactor],
                                      batch2 = sample_table[,colnames(sample_table) == batchfactor_2],
                                      design = model))
    }
  }
}


# Specify structure of DESeq2_analysis_object
setClass(Class = "DESeq2_analysis_object",
         slots = c(results="data.frame", DE_genes="list", Number_DE_genes="list"))


# Wrapper Function to perform DESeq2 differential testing
DEAnalysis <- function(input = dds,
                       condition,
                       alpha = 0.05,
                       lfcThreshold = 0,
                       sigFC = 2,
                       multiple_testing = "IHW",
                       independentFiltering= TRUE,
                       shrinkage = TRUE,
                       shrinkType = "normal"){
  # create results_list
  results_list <- list()
  # print parameters
  results_list$parameters <-list(multiple_testing = multiple_testing,
                                 p_value_threshold = alpha,
                                 log2_FC_threshold = lfcThreshold,
                                 shrinkage = shrinkage,
                                 shrinkage_type = shrinkType)
  # Run results() function on comparisons defined in comparison table
  for (i in 1:nrow(comparison_table)){
    # create DE_object
    DE_object <- new(Class = "DESeq2_analysis_object")
    # IHW
    if (multiple_testing=="IHW") {
      res_deseq_lfc <- results(input,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               filterFun = ihw,
                               altHypothesis = "greaterAbs")
      # Independent Filtering
    }else {
      res_deseq_lfc <- results(input,
                               contrast = c(condition,
                                            paste(comparison_table$comparison[i]),
                                            paste(comparison_table$control[i])),
                               lfcThreshold = lfcThreshold,
                               alpha = alpha,
                               independentFiltering = independentFiltering,
                               altHypothesis = "greaterAbs",
                               pAdjustMethod= multiple_testing)
    }
    if(shrinkage == TRUE){
      res_deseq_lfc <- lfcShrink(input,
                                 contrast = c(condition,
                                              paste(comparison_table$comparison[i]),
                                              paste(comparison_table$control[i])),
                                 res=res_deseq_lfc,
                                 type = shrinkType)
    }
    res_deseq_lfc <- as.data.frame(res_deseq_lfc)
    # indicate significant DE genes
    res_deseq_lfc$regulation <- ifelse(!is.na(res_deseq_lfc$padj)&
                                         res_deseq_lfc$padj <= alpha&
                                         res_deseq_lfc$log2FoldChange > log(sigFC,2),
                                       "up",
                                       ifelse(!is.na(res_deseq_lfc$padj)&
                                                res_deseq_lfc$padj <= alpha&
                                                res_deseq_lfc$log2FoldChange < -log(sigFC,2),
                                              "down",
                                              "n.s."))
    # add gene annotation to results table
    res_deseq_lfc$GENEID <- row.names(res_deseq_lfc) # ensembl-IDs as row names
    res_deseq_lfc <- merge(res_deseq_lfc,
                           norm_anno[,c("GENEID",
                                        "SYMBOL",
                                        "GENETYPE",
                                        "DESCRIPTION",
                                        "CHR")],
                           by = "GENEID")
    row.names(res_deseq_lfc) <- res_deseq_lfc$GENEID
    res_deseq_lfc$comparison<-paste(comparison_table$comparison[i]," vs ",comparison_table$control[i],
                                    sep="")
    # re-order results table
    if (multiple_testing=="IHW") {
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj",
                                      "weight")]
    }else{
      res_deseq_lfc<-res_deseq_lfc[,c("GENEID",
                                      "SYMBOL",
                                      "GENETYPE",
                                      "DESCRIPTION",
                                      "CHR",
                                      "comparison",
                                      "regulation",
                                      "baseMean",
                                      "log2FoldChange",
                                      "lfcSE",
                                      "stat",
                                      "pvalue",
                                      "padj")]
    }
    # print result table
    DE_object@results <- res_deseq_lfc
    # print DE genes in seperate tables
    DE_object@DE_genes <- list(up_regulated_Genes = res_deseq_lfc[res_deseq_lfc$regulation =="up",],
                               down_regulated_Genes= res_deseq_lfc[res_deseq_lfc$regulation =="down",])
    # print the numbers of DE genes
    DE_object@Number_DE_genes <- list(up_regulated_Genes = nrow(DE_object@DE_genes$up_regulated_Genes),
                                      down_regulated_Genes= nrow(DE_object@DE_genes$down_regulated_Genes))
    # write DE_object into results_list
    results_list[[paste(comparison_table$comparison[i], "vs", comparison_table$control[i], sep="_")]] <- DE_object
  }
  return(results_list)
}

# Function: Union of DEgenes
uDEG <- function(input = DEresults, keyType = "Ensembl", comparisons){
  uDEGs <- NULL
  tmp <- input[names(input) %in% comparisons]
  for(i in 1:length(comparisons)){
    DEGs <- as.data.frame(tmp[[i]]@results[tmp[[i]]@results$regulation %in% c("up","down"),])
    if(!keyType %in% c("Ensembl", "Symbol")){stop("keyType should be one of Ensembl or Symbol")}
    if(keyType == "Ensembl"){
      uDEGs <- unique(c(uDEGs, DEGs$GENEID))
    } else if (keyType == "Symbol"){
      uDEGs <- unique(c(uDEGs, DEGs$SYMBOL))
    }
  }
  uDEGs
}


# Function: Venn Diagram
plotVenn <- function(comparisons,
                     regulation=NULL){
  venn <- NULL
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    comp <- as.data.frame(res[[i]]@results)
    if(is.null(regulation)){
      DE <- ifelse(comp$regulation %in% c("up","down"), 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], "up&down", sep=": ")
    } else {
      DE <- ifelse(comp$regulation == regulation, 1, 0)
      venn <- cbind(venn, DE)
      colnames(venn)[i]<- paste(names(res)[[i]], regulation, sep=": ")
    }
    
  }
  vennDiagram(venn,cex = 1, counts.col = "blue")
}

# Ratio plot function
plotRatios <- function(comp1, comp2){
  U <- NULL
  c <- c(comp1,comp2)
  U <- uDEG(c)
  Ratio <- NULL
  for(i in 1:length(c)){
    tmp <- DEresults[names(DEresults) %in% c]
    comp <- as.data.frame(tmp[[i]]@results)
    DE <- as.data.frame(comp[rownames(comp) %in% U,])
    Ratio <- as.data.frame(cbind(Ratio,DE$log2FoldChange))
  }
  colnames(Ratio)<- c
  rownames(Ratio) <- U
  ggplot(Ratio, aes(x=Ratio[,1], y=Ratio[,2])) +
    geom_point(colour = "grey", size = 1.5) + 
    theme_bw() +
    xlab(comp1)+
    ylab(comp2) +
    geom_abline(slope = c(-1,1),intercept = 0, colour="grey") +
    geom_hline(yintercept = c(0,log(2,2),-(log(2,2))))+
    geom_hline(yintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))))+
    geom_vline(xintercept = c(log(2,2),-(log(2,2))), colour="firebrick1")+
    theme(text = element_text(size=10))+
    ggtitle(paste(comp1," vs ",comp2,": ",length(U)," DE genes",sep=""))
}

# GO & KEGG enrichment across comparisons
compareGSEA <- function(comparisons, 
                        organism, 
                        GeneSets =c("GO","KEGG"),
                        ontology= "BP",
                        pCorrection = "bonferroni", # choose the p-value adjustment method
                        pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                        qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                        showMax = 20){
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {stop("Wrong Organism. Select mouse or human.")}
  
  ENTREZlist <-  list()
  for(i in 1:length(comparisons)){
    res <- DEresults[names(DEresults) %in% comparisons]
    DE_up <- as.data.frame(res[[i]]@DE_genes$up_regulated_Genes)$SYMBOL
    entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    DE_down <- as.data.frame(res[[i]]@DE_genes$down_regulated_Genes)$SYMBOL
    entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID  
    x <- setNames(list(entrez_up, entrez_down),
                  c(paste(names(res[i]),"_up",sep=""), 
                    paste(names(res[i]),"_down",sep="")))
    ENTREZlist <- c(ENTREZlist,x)
  }
  
  list <- list()
  
  # Compare the Clusters regarding their GO enrichment  
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    CompareClusters_GO <- compareCluster(geneCluster = ENTREZlist, 
                                         fun = "enrichGO",  
                                         universe = universe_Entrez,
                                         OrgDb = OrgDb,
                                         ont = ontology, 
                                         pvalueCutoff  = pvalueCutoff, 
                                         pAdjustMethod = pCorrection, 
                                         qvalueCutoff  = pvalueCutoff,  
                                         readable      = T)
    list$GOresults <- as.data.frame(CompareClusters_GO)
    list$GOplot <- clusterProfiler::dotplot(CompareClusters_GO, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse"){org = "mmu"} 
    if(organism == "human"){org = "hsa"}
    
    # Compare the Clusters regarding their KEGG enrichment  
    CompareClusters_KEGG <- compareCluster(geneCluster = ENTREZlist, 
                                           fun = "enrichKEGG",  
                                           universe = universe_Entrez,
                                           organism = org, 
                                           pvalueCutoff  = pvalueCutoff, 
                                           pAdjustMethod = pCorrection, 
                                           qvalueCutoff  = pvalueCutoff)
    list$KEGGresults <- as.data.frame(CompareClusters_KEGG)
    list$KEGGplot <- clusterProfiler::dotplot(CompareClusters_KEGG, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  list
}

# Custom function to plot baseMean versus fold change
plotMA <- function(comparison, 
                   ylim=c(-2,2),  
                   padjThreshold=0.05,
                   xlab = "mean of normalized counts", 
                   ylab = expression(log[2]~fold~change),
                   log = "x", 
                   cex=0.45){
  x <- as.data.frame(DEresults[[comparison]]@results)
  if (!(is.data.frame(x) && all(c("baseMean", "log2FoldChange") %in% colnames(x)))){
    stop("'x' must be a data frame with columns named 'baseMean', 'log2FoldChange'.")
  }
  col = ifelse(x$padj>=padjThreshold, "gray32", "red3")
  py = x$log2FoldChange
  if(missing(ylim)){
    ylim = c(-1,1) * quantile(abs(py[is.finite(py)]), probs=0.99) * 1.1
  }
  plot(x=x$baseMean, 
       y=pmax(ylim[1], pmin(ylim[2], py)),
       log=log, 
       pch=ifelse(py<ylim[1], 6, ifelse(py>ylim[2], 2, 16)),
       cex=cex, 
       col=col, 
       xlab=xlab, 
       ylab=ylab, 
       ylim=ylim,
       main=comparison)
  abline(h=0, lwd=4, col="#ff000080")
  abline(h=c(-1,1), lwd=2, col="dodgerblue")
}

# Plot p value distribution  
plotPvalues <- function(comparison){
  res <- as.data.frame(DEresults[[comparison]]@results)
  ggplot(na.omit(res), aes(x=pvalue)) + 
    geom_histogram(aes(y=..count..),               
                   binwidth = 0.01) +
    theme_bw()+
    ggtitle(paste("p value histogram of: ",comparison,sep=""))
}

# Function to plot heatmaps of DE genes based on pheatmap
plotDEHeatmap_onlyDEsubset <- function(comparison,
                                       factor,
                                       conditions="all",
                                       show_rownames = FALSE,
                                       cluster_cols = FALSE){
  geneset <- c(DEresults[[comparison]]@DE_genes$up_regulated_Genes$GENEID, 
               DEresults[[comparison]]@DE_genes$down_regulated_Genes$GENEID)
  
  input <- norm_anno[norm_anno$GENEID %in% geneset,]
  rownames(input) <- paste(input$GENEID, ":", input$SYMBOL, sep="")
  
  if(conditions[1] == "all"){
    input <- input[,colnames(input) %in% sample_table$ID]
    input_scale <- t(scale(t(input)))
    input_scale <- input_scale[,order(sample_table[[plot_order]], decreasing = FALSE)]
  } else {
    input <- input[,colnames(input) %in% sample_table[as.vector(sample_table[[factor]]) %in% conditions,]$ID,]
    input_scale <- t(scale(t(input)))
  }
  
  pheatmap(input_scale,
           main=paste("Heatmap of significant DE genes in: ",comparison,sep=""),
           show_rownames=show_rownames,
           show_colnames=TRUE,
           cluster_cols = cluster_cols, 
           fontsize = 7,
           annotation_col = plot_annotation,
           annotation_colors = ann_colors,
           breaks = scaleColors(data = input_scale, maxvalue = 2)[["breaks"]], 
           color = scaleColors(data = input_scale, maxvalue = 2)[["color"]])
}

# Volcano Plot function
plotVolcano <-  function(comparison,
                         labelnum=20, input = DEresults){
  
  # specify labeling
  upDE <-  as.data.frame(input[[comparison]]@results[input[[comparison]]@results$regulation =="up",])
  FClabel_up <- upDE[order(abs(upDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_up)>labelnum){
    FClabel_up <- as.character(FClabel_up[c(1:labelnum),"GENEID"])
  } else {
    FClabel_up <- as.character(FClabel_up$GENEID)}
  plabel_up <- upDE[order(upDE$padj, decreasing = FALSE),]
  if(nrow(plabel_up)>labelnum){
    plabel_up <- as.character(plabel_up[c(1:labelnum),"GENEID"])
  } else {
    plabel_up <- as.character(plabel_up$GENEID)}
  
  downDE <-  as.data.frame(input[[comparison]]@results[input[[comparison]]@results$regulation =="down",])
  FClabel_down <- downDE[order(abs(downDE$log2FoldChange), decreasing = TRUE),]
  if(nrow(FClabel_down)>labelnum){
    FClabel_down <- as.character(FClabel_down[c(1:labelnum),"GENEID"])
  } else {
    FClabel_down <- as.character(FClabel_down$GENEID)}
  plabel_down <- downDE[order(downDE$padj, decreasing = FALSE),]
  if(nrow(plabel_down)>labelnum){
    plabel_down <- as.character(plabel_down[c(1:labelnum),"GENEID"])
  } else {
    plabel_down <- as.character(plabel_down$GENEID)}
  
  
  label<- unique(c(FClabel_up, plabel_up, FClabel_down, plabel_down))
  
  data <- input[[comparison]]@results
  data$label<- ifelse(data$GENEID %in% label == "TRUE",as.character(data$SYMBOL), "")
  data <- data[,colnames(data) %in% c("label", "log2FoldChange", "padj", "regulation")]
  
  # Volcano Plot
  ggplot(data=na.omit(data), aes(x=log2FoldChange, y=-log10(padj), colour=regulation)) +
    geom_point(alpha=0.4, size=1.75) +
    scale_color_manual(values=c("cornflowerblue","grey", "firebrick"))+
    scale_x_continuous() +
    scale_y_continuous() +
    xlab("log2(FoldChange)") +
    ylab("-log10(padj)") +
    geom_vline(xintercept = 0, colour="black")+
    geom_vline(xintercept = c(-log(2,2),log(2,2)), colour="red")+
    geom_hline(yintercept=-log(0.05,10),colour="red")+
    geom_text_repel(data=na.omit(data[!data$label =="",]),aes(label=label), size=3)+
    guides(colour=FALSE) +
    ggtitle(paste("Volcano Plot of: ",comparison,sep="")) +
    theme_bw()
}

# GSEA function
GSEA <-  function(comparison,
                  organism,
                  
                  
                  DE_results = DEresults,
                  GeneSets =c("GO","KEGG","DO","Hallmark","canonicalPathways","Motifs","ImmunoSignatures"),
                  
                  GOntology = "BP",
                  pCorrection = "bonferroni", # choose the p-value adjustment method
                  pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                  qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                  showMax = 20,
                  font.size = 8){
  
  results <- list()
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {print("Wrong Organism. Select mouse or human.")}
  
  res <- DE_results[[comparison]]
  DE_up <- as.data.frame(res@DE_genes$up_regulated_Genes)$SYMBOL
  entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
  DE_down <- as.data.frame(res@DE_genes$down_regulated_Genes)$SYMBOL
  entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
  
  # GO enrichment ###############################
  
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    if(length(entrez_up)<20){
      print("Too few upregulated genes for GO enrichment (<20)")
      results$GO_up <- "Too few upregulated genes for GO enrichment (<20)"
    }else{
      eGO_up <- enrichGO(gene = entrez_up,
                         universe = universe_Entrez,
                         OrgDb = OrgDb,
                         ont = GOntology,
                         pAdjustMethod = pCorrection,
                         pvalueCutoff  = pvalueCutoff,
                         qvalueCutoff  = qvalueCutoff,
                         readable      = T)
      
      results$GOup <- as.data.frame(eGO_up)
      if(nrow(results$GOup)<1){
        results$GOup_plot <- "No GO enrichment for upregulated genes"
      }else{
        results$GOup_plot <- clusterProfiler::dotplot(eGO_up, 
                                                      showCategory = showMax, 
                                                      font.size= font.size, 
                                                      title = paste("GO enrichment for genes upregulated in: ", comparison,sep="")
        )
      }
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for GO enrichment (<20)")
      results$GO_down <- "Too few downregulated genes for GO enrichment (<20)"
    }else{
      eGO_down <- enrichGO(gene = entrez_down,
                           universe = universe_Entrez,
                           OrgDb = OrgDb,
                           ont = GOntology,
                           pAdjustMethod = pCorrection,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff  = qvalueCutoff,
                           readable      = T)
      
      results$GOdown <- as.data.frame(eGO_down)
      if(nrow(results$GOdown)<1){
        results$GOdown_plot <- "No GO enrichment for downregulated genes"
      }else{
        results$GOdown_plot <- clusterProfiler::dotplot(eGO_down, 
                                                        showCategory = showMax, 
                                                        font.size= font.size, 
                                                        title = paste("GO enrichment for genes downregulated in: ", comparison,sep="")
        )
      }
    }
  }
  
  
  # KEGG enrichment ##########################################
  
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse") {org = "mmu"} 
    if(organism == "human"){org = "hsa"}
    
    if(length(entrez_up)<20){
      print("Too few upregulated genes for KEGG enrichment (<20)")
      results$KEGG_up <- "Too few upregulated genes for KEGG enrichment (<20)"
    }else{
      eKEGG_up <- enrichKEGG(gene = entrez_up, 
                             organism = org,
                             universe = universe_Entrez, 
                             pAdjustMethod = pCorrection,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff = qvalueCutoff)
      
      results$KEGGup <- as.data.frame(eKEGG_up)
      if(nrow(results$KEGGup)<1){
        results$KEGGup_plot <- "No KEGG enrichment for upregulated genes"
      }else{
        results$KEGGup_plot <- clusterProfiler::dotplot(eKEGG_up,  
                                                        showCategory = showMax, 
                                                        font.size= font.size, 
                                                        title = paste("KEGG enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for KEGG enrichment (<20)")
      results$KEGG_down <- "Too few downregulated genes for KEGG enrichment (<20)"
    } else{
      eKEGG_down <- enrichKEGG(gene = entrez_down, 
                               organism = org,
                               universe = universe_Entrez, 
                               pAdjustMethod = pCorrection,
                               pvalueCutoff  = pvalueCutoff,
                               qvalueCutoff = qvalueCutoff)
      
      results$KEGGdown <- as.data.frame(eKEGG_down)
      if(nrow(results$KEGGdown)<1){
        results$KEGGdown_plot <- "No KEGG enrichment for downregulated genes"
      }else{
        results$KEGGdown_plot <- clusterProfiler::dotplot(eKEGG_down,
                                                          showCategory = showMax,
                                                          font.size= font.size,
                                                          title = paste("KEGG enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  if("Hallmark" %in% GeneSets |
     "DO" %in% GeneSets |
     "canonicalPathways" %in% GeneSets|
     "ImmunoSignatures" %in% GeneSets |
     "Motifs" %in% GeneSets){
    if(organism == "mouse"){
      
      entrez_up_hsa <- as.character(getLDS(attributes = c("mgi_symbol"),
                                           filters = "mgi_symbol",
                                           values = DE_up,
                                           mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                                           attributesL = c("entrezgene_id"),
                                           martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                           uniqueRows=T)[,2])
      entrez_down_hsa <- getLDS(attributes = c("mgi_symbol"),
                                filters = "mgi_symbol",
                                values = DE_down,
                                mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl"),
                                attributesL = c("entrezgene_id"),
                                martL = useMart("ensembl", dataset = "hsapiens_gene_ensembl"),
                                uniqueRows=T)[,2]
      
    } else if(organism == "human"){
      entrez_up_hsa <- entrez_up
      entrez_down_hsa <- entrez_down
    }
  }
  
  # DO enrichment ########################################
  if("DO" %in% GeneSets){
    print("Performing Disease Ontology enrichment")
    
    if(length(entrez_up)<20){
      print("Too few upregulated genes for DO enrichment (<20)")
      results$DOup <- "Too few upregulated genes for DO enrichment (<20)"
    }else{
      results$DOup <- as.data.frame(enrichDO(gene = entrez_up_hsa,
                                             universe = universe_mouse2human_Entrez,
                                             pAdjustMethod = pCorrection,
                                             pvalueCutoff  = pvalueCutoff,
                                             qvalueCutoff = qvalueCutoff,
                                             minGSSize     = 5,
                                             maxGSSize     = 500,
                                             readable=TRUE))
      if(nrow(results$DOup)>0){results$DOup$Enrichment <- paste("DO enrichment for genes upregulated in ",comparison,sep="")}
    }
    if(length(entrez_down)<20){
      print("Too few downregulated genes for DO enrichment (<20)")
      results$DOdown <- "Too few downregulated genes for DO enrichment (<20)"
    } else{
      results$DOdown <- as.data.frame(enrichDO(gene = entrez_down_hsa,
                                               universe = universe_mouse2human_Entrez,
                                               pAdjustMethod = pCorrection,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff = qvalueCutoff,
                                               minGSSize     = 5,
                                               maxGSSize     = 500,
                                               readable=TRUE))
      if(nrow(results$DOdown)>0){results$DOdown$Enrichment <- paste("DO enrichment for genes downregulated in ",comparison,sep="")}
    }
  }
  
  # Hallmark enrichment ################################
  
  if("Hallmark" %in% GeneSets){
    print("Performing Hallmark enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Hallmark enrichment (<20)")
      results$Hallmark_up <- "Too few upregulated genes for Hallmark enrichment (<20)"
    }else{
      Hallmark_up <- enricher(entrez_up_hsa,
                              TERM2GENE=hallmark_genes,
                              universe = universe_mouse2human_Entrez,  
                              pAdjustMethod = pCorrection,
                              pvalueCutoff  = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
      
      results$HALLMARKup <- as.data.frame(Hallmark_up)
      if(nrow(results$HALLMARKup)<1){
        results$HALLMARKup_plot <- "No Hallmark enrichment for upregulated genes"
      }else{
        results$HALLMARKup_plot <- clusterProfiler::dotplot(Hallmark_up,
                                                            showCategory = showMax,
                                                            font.size= font.size,
                                                            title = paste("Hallmark enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Hallmark enrichment (<20)")
      results$Hallmark_down <- "Too few downregulated genes for Hallmark enrichment (<20)"
    }else{
      Hallmark_down <- enricher(entrez_down_hsa,
                                TERM2GENE=hallmark_genes,
                                universe = universe_mouse2human_Entrez,  
                                pAdjustMethod = pCorrection,
                                pvalueCutoff  = pvalueCutoff,
                                qvalueCutoff = qvalueCutoff)
      
      results$HALLMARKdown <- as.data.frame(Hallmark_down)
      if(nrow(results$HALLMARKdown)<1){
        results$HALLMARKdown_plot <-"No Hallmark enrichment for downregulated genes"
      }else{
        results$HALLMARKdown_plot <- clusterProfiler::dotplot(Hallmark_down,
                                                              showCategory = showMax,
                                                              font.size= font.size,
                                                              title = paste("Hallmark enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  # Canonical Pathway enrichment #############################
  if("canonicalPathways" %in% GeneSets){
    print("Performing Canonical Pathway (C2) enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Canonical Pathway enrichment (<20)")
      results$canonicalPathwaysup <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      results$canonicalPathwaysup <- as.data.frame(enricher(entrez_up_hsa,
                                                            TERM2GENE=canonicalPathway_genes,
                                                            universe = universe_mouse2human_Entrez,
                                                            pAdjustMethod = pCorrection,
                                                            pvalueCutoff  = pvalueCutoff,
                                                            qvalueCutoff = qvalueCutoff))
      if(nrow(results$canonicalPathwaysup)>0){results$canonicalPathwaysup$Enrichment <- paste("Canonical pathway enrichment for genes upregulated in ",comparison,sep="")}
      
    }
    
    if(length(entrez_down_hsa)<20){
      
      print("Too few downregulated genes for canonical pathway  enrichment (<20)")
      results$canonicalPathwaysdown <- "Too few downregulated genes for canonical pathway enrichment (<20)"
    }else{
      results$canonicalPathwaysdown <- as.data.frame(enricher(entrez_down_hsa,
                                                              TERM2GENE=canonicalPathway_genes,
                                                              universe = universe_mouse2human_Entrez,
                                                              pAdjustMethod = pCorrection,
                                                              pvalueCutoff  = pvalueCutoff,
                                                              qvalueCutoff = qvalueCutoff))
      if(nrow(results$canonicalPathwaysdown)>0){results$canonicalPathwaysdown$Enrichment <- paste("Canonical pathway enrichment for genes downregulated in ",comparison,sep="")}
      
    }
  }
  
  # Motif enrichment
  if("Motifs" %in% GeneSets){
    print("Performing Motif enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Motif enrichment (<20)")
      results$Motif_up <- "Too few upregulated genes for Motif enrichment (<20)"
    }else{
      Motif_up <- enricher(entrez_up_hsa,
                           TERM2GENE=motifs,
                           universe = universe_mouse2human_Entrez,  
                           pAdjustMethod = pCorrection,
                           pvalueCutoff  = pvalueCutoff,
                           qvalueCutoff = qvalueCutoff)
      
      results$Motifup <- as.data.frame(Motif_up)
      if(nrow(results$Motifup)<1){
        results$Motifup_plot <- "No Motif enrichment for upregulated genes"
      }else{
        results$Motifup_plot <- clusterProfiler::dotplot(Motif_up,
                                                         showCategory = showMax,
                                                         font.size= font.size,
                                                         title = paste("Motif enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Motif enrichment (<20)")
      results$Motif_down <- "Too few downregulated genes for Motif enrichment (<20)"
    }else{
      Motif_down <- enricher(entrez_down_hsa,
                             TERM2GENE=motifs,
                             universe = universe_mouse2human_Entrez,  
                             pAdjustMethod = pCorrection,
                             pvalueCutoff  = pvalueCutoff,
                             qvalueCutoff = qvalueCutoff)
      
      results$Motifdown <- as.data.frame(Motif_down)
      if(nrow(results$Motifdown)<1){
        results$Motifdown_plot <-"No Motif enrichment for downregulated genes"
      }else{
        results$Motifdown_plot <- clusterProfiler::dotplot(Motif_down,
                                                           showCategory = showMax,
                                                           font.size= font.size,
                                                           title = paste("Motif enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  
  # Immunosignatures enrichment
  if("ImmunoSignatures" %in% GeneSets){
    print("Performing Immunosignature enrichment")
    if(length(entrez_up_hsa)<20){
      print("Too few upregulated genes for Immunosignature enrichment (<20)")
      results$ImmSig_up <- "Too few upregulated genes for Immunosignature enrichment (<20)"
    }else{
      ImmSig_up <- enricher(entrez_up_hsa,
                            TERM2GENE=immuno_genes,
                            universe = universe_mouse2human_Entrez,  
                            pAdjustMethod = pCorrection,
                            pvalueCutoff  = pvalueCutoff,
                            qvalueCutoff = qvalueCutoff)
      
      results$ImmSigup <- as.data.frame(ImmSig_up)
      if(nrow(results$ImmSigup)<1){
        results$ImmSigup_plot <- "No Immunosignature enrichment for upregulated genes"
      }else{
        results$ImmSigup_plot <- clusterProfiler::dotplot(ImmSig_up,
                                                          showCategory = showMax,
                                                          font.size= font.size,
                                                          title = paste("Immunosignature enrichment for genes upregulated in: ",comparison,sep="")
        )
      }
    }
    if(length(entrez_down_hsa)<20){
      print("Too few downregulated genes for Immunosignature enrichment (<20)")
      results$ImmSig_down <- "Too few downregulated genes for Immunosignature enrichment (<20)"
    }else{
      ImmSig_down <- enricher(entrez_down_hsa,
                              TERM2GENE=immuno_genes,
                              universe = universe_mouse2human_Entrez,  
                              pAdjustMethod = pCorrection,
                              pvalueCutoff  = pvalueCutoff,
                              qvalueCutoff = qvalueCutoff)
      
      results$ImmSigdown <- as.data.frame(ImmSig_down)
      if(nrow(results$ImmSigdown)<1){
        results$ImmSigdown_plot <- "No Immunosignature enrichment for downregulated genes"
      }else{
        results$ImmSigdown_plot <- clusterProfiler::dotplot(ImmSig_down,
                                                            showCategory = showMax,
                                                            font.size= font.size,
                                                            title = paste("Immunosignature enrichment for genes downregulated in: ",comparison,sep="")
        )
      }
    }
  }
  results
}

# Function to plot a heatmap of genes responsible for gene set enrichment
plotGSEAHeatmap<-function(GSEA_result,
                          GeneSet, 
                          term,
                          regulation){
  xterm <- paste("^", term, "$", sep="") 
  tmp <- GSEA_result[grep(xterm,GSEA_result$Description),]
  gene.list <- unique(unlist(strsplit(tmp$geneID, split = "/")))
  
  if(GeneSet == "KEGG"){
    gene.list <- bitr(gene.list, 
                      fromType = "ENTREZID", 
                      toType="SYMBOL", 
                      OrgDb="org.Mm.eg.db")[,2]
  }
  
  if(GeneSet == "HALLMARK" | GeneSet == "ImmunoSignatures" | GeneSet == "Motifs"){
    gene.list <- getLDS(attributes = c("hgnc_symbol"), 
                        filters = "hgnc_symbol", 
                        values = gene.list, 
                        mart = human, 
                        attributesL = c("mgi_symbol"), 
                        martL = mouse, 
                        uniqueRows=T)[,2]
  }
  
  plotHeatmap(geneset = gene.list,
              keyType = "Symbol",
              title = paste("Heatmap of genes responsible for enrichment of term: ",term,", in ",deparse(substitute(GSEA_result)),sep=""),
              show_rownames = TRUE,
              cluster_cols = FALSE)
}


compareGSEA <- function(comparisons,
                        DE_results = DEresults,
                        organism,
                        GeneSets =c("GO","KEGG"),
                        ontology= "BP",
                        pCorrection = "bonferroni", # choose the p-value adjustment method
                        pvalueCutoff = 0.05, # set the unadj. or adj. p-value cutoff (depending on correction method)
                        qvalueCutoff = 0.05, # set the q-value cutoff (FDR corrected)
                        showMax = 20){
  
  if(organism == "mouse") {
    OrgDb = org.Mm.eg.db
  } else if(organism == "human"){
    OrgDb = org.Hs.eg.db
  } else {stop("Wrong Organism. Select mouse or human.")}
  
  ENTREZlist <-  list()
  for(i in 1:length(comparisons)){
    res <- DE_results[names(DE_results) %in% comparisons]
    DE_up <- as.data.frame(res[[i]]@DE_genes$up_regulated_Genes)$SYMBOL
    entrez_up <- bitr(DE_up, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    DE_down <- as.data.frame(res[[i]]@DE_genes$down_regulated_Genes)$SYMBOL
    entrez_down <- bitr(DE_down, fromType = "SYMBOL", toType="ENTREZID", OrgDb=OrgDb)$ENTREZID
    x <- setNames(list(entrez_up, entrez_down),
                  c(paste(names(res[i]),"_up",sep=""),
                    paste(names(res[i]),"_down",sep="")))
    ENTREZlist <- c(ENTREZlist,x)
  }
  
  list <- list()
  
  # Compare the Clusters regarding their GO enrichment
  if("GO" %in% GeneSets){
    print("Performing GO enrichment")
    CompareClusters_GO <- compareCluster(geneCluster = ENTREZlist,
                                         fun = "enrichGO",
                                         universe = universe_Entrez,
                                         OrgDb = OrgDb,
                                         ont = ontology,
                                         pvalueCutoff  = pvalueCutoff,
                                         pAdjustMethod = pCorrection,
                                         qvalueCutoff  = pvalueCutoff,
                                         readable      = T)
    list$GOresults <- as.data.frame(CompareClusters_GO)
    list$GOplot <- clusterProfiler::dotplot(CompareClusters_GO, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  
  if("KEGG" %in% GeneSets){
    print("Performing KEGG enrichment")
    
    if(organism == "mouse"){org = "mmu"}
    if(organism == "human"){org = "hsa"}
    
    # Compare the Clusters regarding their KEGG enrichment
    CompareClusters_KEGG <- compareCluster(geneCluster = ENTREZlist,
                                           fun = "enrichKEGG",
                                           universe = universe_Entrez,
                                           organism = org,
                                           pvalueCutoff  = pvalueCutoff,
                                           pAdjustMethod = pCorrection,
                                           qvalueCutoff  = pvalueCutoff)
    list$KEGGresults <- as.data.frame(CompareClusters_KEGG)
    list$KEGGplot <- clusterProfiler::dotplot(CompareClusters_KEGG, showCategory = showMax, by = "geneRatio", font.size=10)
  }
  list
}

dotplotGSEA <- function(x,
                        show=25,
                        font.size=10,
                        title.size=10,
                        title.width=100,
                        order="count"){
  if(nrow(x)<1){
    print("No enrichment found.")
  }else{
    x <- if(nrow(x)>show){x[c(1:show),]}else{x}
    if(order=="padj"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
      x <- x[order(x$p.adjust,decreasing=TRUE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
    }
    if(order=="count"){
      x <- x[order(x$Count,decreasing=FALSE),]
      x$Description <- factor(x$Description, levels = unique(x$Description))
      x$GeneRatio <- factor(x$GeneRatio, levels = unique(x$GeneRatio))
    }
    ggplot(x, aes(x = GeneRatio, y = Description, color = p.adjust)) +
      geom_point(aes(size = Count)) +
      scale_colour_gradientn(colours=c('red',
                                       'orange',
                                       'darkblue',
                                       'darkblue'),
                             limits=c(0,1),
                             values   = c(0,0.05,0.2,0.5,1),
                             breaks   = c(0.05,0.2,1),
                             labels = format(c(0.05,0.2,1))) +
      ylab(NULL) +
      ggtitle(paste(strwrap(unique(x$Enrichment), width=title.width), collapse = "\n"))+
      theme_bw() +
      theme(text = element_text(size=font.size),
            plot.title = element_text(size=title.size))
  }
}


