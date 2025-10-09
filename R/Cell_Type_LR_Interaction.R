#' @export
Ligand_Receptor_Interaction_Scores <- function(single_cell_metadata,
                                               single_cell_gex,
                                               sample_colname,
                                               celltype_colname,
                                               subset_data = NULL,
                                               lr_database = NULL,
                                               celltype_to_exclude = NULL){


  require(limma)
  require(usethis)
  require(reshape2)
  require(dplyr)
  require(Matrix)

  ## Check Validity of Parameters and Objects
  #########################################################################
  # Check for presence of required parameters
  if(is.null(single_cell_metadata)){
    warning("Please provide the a single-cell metadata object with the sample and cell type annotation under the `single_cell_metadata` parameter as explained in the vignette")
    return()
  }
  if(is.null(single_cell_gex)){
    warning("Please provide the a single-cell gene expression matrix with genes in rows and cells in columns under the `single_cell_gex` parameter as explained in the vignette")
    return()
  }

  # Check that cells are identical between single_cell_metadata and single_cell_gex
  if(!identical(rownames(single_cell_metadata), colnames(single_cell_gex))){
    common_cells <- intersect(rownames(single_cell_metadata), colnames(single_cell_gex))
    if(length(common_cells)==0){
      warning("No common cells found between the row.names of `single_cell_metadata` and the col.names of `single_cell_gex`")
      return()
    } else {
      warning(paste0("Cells and/or order of cells were not matching between `single_cell_metadata` and `single_cell_gex`. Subsetting down to ",length(common_cells)," common cells"))
      # Subset and reorder both metadata and gene expression matrix
      single_cell_metadata <- single_cell_metadata[common_cells, , drop = FALSE]
      single_cell_gex <- single_cell_gex[, common_cells, drop = FALSE]
    }
  }

  if(is.null(sample_colname)){
    warning("Please provide the name of the column containing the sample identifier in the `single_cell_metadata` object under the `sample_colname` parameter")
    return()
  } else {
    if(length(which(colnames(single_cell_metadata) %in% c(sample_colname))) != 1){
      warning("The `sample_colname` identifier was not found in the column names of the `single_cell_metadata` object")
      return()
    }
  }

  if(is.null(celltype_colname)){
    warning("Please provide the name of the column containing the cell type identifier in the `single_cell_metadata` object under the `celltype_colname` parameter")
    return()
  } else {
    if(length(which(colnames(single_cell_metadata) %in% c(celltype_colname))) != 1){
      warning("The `celltype_colname` identifier was not found in the column names of the `single_cell_metadata` object")
      return()
    }
  }
  # Load or check format of lr_database
  if(is.null(lr_database)){
    #data("LR_database", package = "scAvengeR")
    #lr_database <- get("LR_database", envir = asNamespace("scAvengeR"))
    lr_database <- scAvengeR::lr_database
  } else {
    if(length(which(c("receptor","ligand") %in% colnames(lr_database))) != 2){
      warning("The lr_database data.frame should contain two columns named `receptor` and `ligand` containing the pairs of ligand-receptors")
      return()
    } else {
      lr_database <- lr_database[,c("receptor","ligand")]
    }
  }

  # Check that ligand and Receptor are found in Gene Expression data
  all_genes <- unique(c(lr_database$receptor,lr_database$ligand))
  common_genes <- intersect(all_genes,rownames(single_cell_gex))
  warning(paste0(round((length(common_genes) / length(all_genes)) * 100,digits = 1)," % (",length(common_genes),"/",length(all_genes)," genes) of the LR genes were found in the `single_cell_gex` object"))

  if(length(common_genes)==0){
    warning("No common genes found between the `single_cell_gex` and `lr_database` objects")
    return()
  } else {
    single_cell_gex <- single_cell_gex[which(rownames(single_cell_gex) %in% common_genes),]
  }
  #########################################################################

  ## Subset Data
  #########################################################################
  if(!is.null(subset_data)){
    if (!is.list(subset_data)) {
      warning("The `subset_data` parameter should be a list")
      return()
    }
    if(length(which(names(subset_data)  %in% colnames(single_cell_metadata) ==F))>0){
      warning("The selected subset_data parameters were not found in col.names of the single-cell object")
      return()
    }
    idx_to_keep <- 1:nrow(single_cell_metadata)
    for(l in 1:length(subset_data)){
      idx_to_keep <- intersect(idx_to_keep,which(single_cell_metadata[,names(subset_data)[l]] %in% subset_data[[l]]))
    }
    single_cell_metadata <- single_cell_metadata[idx_to_keep,]
    single_cell_gex <- single_cell_gex[,idx_to_keep]
  }
  #########################################################################

  ## Compute interaction scores
  #########################################################################
  idx <- which(lr_database$receptor %in% rownames(single_cell_gex)[which(rownames(single_cell_gex) %in% unique(c(lr_database$ligand,lr_database$receptor)))])
  lr_database <- lr_database[idx,]
  idx <- which(lr_database$ligand %in% rownames(single_cell_gex)[which(rownames(single_cell_gex) %in% unique(c(lr_database$ligand,lr_database$receptor)))])
  lr_database <- lr_database[idx,]

  # Filter-out non expressed genes
  suppressWarnings(idx <- which(lr_database$receptor %in% names(which(apply(single_cell_gex[unique(c(lr_database$ligand,lr_database$receptor)),],1,sum)==0))))
  if(length(idx) > 0){
    lr_database <- lr_database[-idx,]
  }
  suppressWarnings(idx <- which(lr_database$ligand %in% names(which(apply(single_cell_gex[unique(c(lr_database$ligand,lr_database$receptor)),],1,sum)==0))))
  if(length(idx) > 0){
    lr_database <- lr_database[-idx,]
  }

  comb_LR <- c(paste(lr_database$ligand,lr_database$receptor,sep=":"),paste(lr_database$receptor,lr_database$ligand,sep=":"))

  cell_types <- unique(as.character(single_cell_metadata[,celltype_colname]))
  if(!is.null(celltype_to_exclude)){
    if(length(which(cell_types %in% celltype_to_exclude))>0){
      cell_types <- cell_types[-which(cell_types %in% celltype_to_exclude)]
    }
  }
  comb_CT <- as.data.frame(utils::combn(cell_types,2))

  # Create the data.frame to collect the stats for all possibilities
  gene_split <- do.call(rbind, strsplit(comb_LR, split = ":"))
  colnames(gene_split) <- c("gene1", "gene2")
  poss_comb_list <- lapply(1:ncol(comb_CT), function(i) {
    data.frame(
      celltype1 = comb_CT[1, i],
      gene1 = gene_split[, 1],
      celltype2 = comb_CT[2, i],
      gene2 = gene_split[, 2]
    )
  })
  poss_comb <- do.call(rbind, poss_comb_list)
  rm(poss_comb_list);invisible(gc())

  # Compute Proportions of cell types per sample
  proportion_cell_type <- t(as.data.frame(unclass(prop.table(table(single_cell_metadata[,which(colnames(single_cell_metadata) %in% c(sample_colname,celltype_colname))]),2))))

  # Filter
  samples <- unique(single_cell_metadata[,sample_colname])
  genes <- intersect(unique(c(poss_comb$gene1,poss_comb$gene2)),unique(c(lr_database$ligand,lr_database$receptor)))
  single_cell_gex <- single_cell_gex[genes,]
  poss_comb <- poss_comb[which((poss_comb$gene1 %in% genes)&(poss_comb$gene2 %in% genes)),]

  # Create unique LR ID
  poss_comb$unique_id <- paste0(poss_comb$celltype1,"_",poss_comb$gene1,":",
                                poss_comb$celltype2,"_",poss_comb$gene2)

  interaction_score <- as.data.frame(matrix(nrow=nrow(poss_comb),ncol=length(samples)))
  rownames(interaction_score) <- poss_comb$unique_id
  colnames(interaction_score) <- samples

  # Create 4 matrices with Cell Type Proportion and Averaged Gene expression per patient and cell type for every possible combination

  # Prop_CT1
  ##########
  Prop_CT1 <- as.data.frame(t(proportion_cell_type[match(colnames(interaction_score),rownames(proportion_cell_type)),
                                                   match(poss_comb$celltype1,colnames(proportion_cell_type))]))
  rownames(Prop_CT1) <- rownames(interaction_score)

  # Prop_CT2
  ##########
  Prop_CT2 <- as.data.frame(t(proportion_cell_type[match(colnames(interaction_score),rownames(proportion_cell_type)),
                                                   match(poss_comb$celltype2,colnames(proportion_cell_type))]))
  rownames(Prop_CT2) <- rownames(interaction_score)

  # GEX_CT1 & GEX_CT2
  ##########
  GEX_CT1 <- interaction_score
  GEX_CT2 <- interaction_score

  template <- matrix(nrow=nrow(single_cell_gex),ncol=length(cell_types),0)
  rownames(template) <- rownames(single_cell_gex)
  colnames(template) <- cell_types

  for(i in 1:ncol(interaction_score)){
    idx <- which(single_cell_metadata[,sample_colname]==colnames(interaction_score)[i])

    tmp <- as.data.frame(single_cell_gex[,idx])
    tmp <- as.data.frame(t(tmp))

    tmp <- cbind(tmp,"Cell_Type"=single_cell_metadata[,celltype_colname][idx])
    tmp <- suppressWarnings(reshape2::melt(tmp, id.vars = "Cell_Type"))

    avg <- tmp %>%
      group_by(variable, Cell_Type) %>%
      summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

    for(c in 1:ncol(template)){
      tmp <- avg[which(avg$Cell_Type==colnames(template)[c]),]
      if(nrow(tmp) > 0){
        template[,c] <- tmp$value[match(rownames(template),tmp$variable)]
        idx1 <- which(poss_comb$celltype1==colnames(template)[c])
        GEX_CT1[idx1,i] <- template[match(poss_comb$gene1[idx1],rownames(template)),c]
        idx2 <- which(poss_comb$celltype2==colnames(template)[c])
        GEX_CT2[idx2,i] <- template[match(poss_comb$gene2[idx2],rownames(template)),c]
      }
    }

    # Reinitialize template matrix
    template <- matrix(nrow=nrow(single_cell_gex),ncol=length(cell_types),0)
    rownames(template) <- rownames(single_cell_gex)
    colnames(template) <- cell_types
  }

  # Compute Interaction Score by multiplying the 4 matrices
  interaction_score <- Prop_CT1*GEX_CT1*Prop_CT2*GEX_CT2

  # When the cell type was not present, the pseudobulk gene expression is obviously NA, so change it to 0
  interaction_score[is.na(interaction_score)] <- 0

  #save(interaction_score,poss_comb,GEX_CT1,GEX_CT2,Prop_CT1,Prop_CT2,file="./Interaction_objects_NeoTIL_v2.rda")
  returned_object <- list("interaction_score"=interaction_score,
                          "possible_combination"=poss_comb,
                          "proportion_cell_type"=proportion_cell_type)
  return(returned_object)
}

#' @export
Ligand_Receptor_Differential_Analysis <- function(interaction_score_object,
                                                  group_up,
                                                  label_up=NULL,
                                                  group_dn,
                                                  label_dn=NULL,
                                                  statistics = c("ttest","wilcoxon","lm"),
                                                  lr_pathways_to_plot = NULL,
                                                  color_pathways = NULL,
                                                  cell_type_groups){

  require(limma)
  require(usethis)
  require(reshape2)
  require(dplyr)
  require(Matrix)
  require(pbapply)
  require(ggplot2)
  require(patchwork)
  require(RColorBrewer)
  require(circlize)
  require(scales)
  require(colorspace)
  require(ComplexHeatmap)

  ## Check Validity of Parameters and Objects
  #########################################################################
  # Check for presence of required parameters
  if(is.null(interaction_score_object)){
    warning("Please provide the a `interaction_score_object` object that is the output of the `Ligand_Receptor_Interaction_Scores` function in this package.")
    return()
  } else {
    if(length(which(names(interaction_score_object) %in% c("interaction_score","possible_combination"))) != 2){
      warning("Please provide a valid `interaction_score_object` object. The `interaction_score` and/or `possible_combination` objects were not found in `interaction_score_object`.")
      return()
    }
    interaction_score <- interaction_score_object$interaction_score
    poss_comb <- interaction_score_object$possible_combination
  }

  if(is.null(group_up)){
    warning("Please provide the `group_up` parameter which is a vector of sample names found in the colnames of the `interaction_score` data.frame of the `interaction_score_object` object.")
    return()
  } else {
    group_up <- group_up[!is.na(which(group_up %in% colnames(interaction_score)))]
    if(length(group_up)==0){
      warning("The values in `group_up` do not intersect with the colnames of the `interaction_score` data.frame of the `interaction_score_object` object.")
      return()
    }
    if(length(group_up)==1){
      warning("Only one value in `group_up` found to intersect with the colnames of the `interaction_score` data.frame of the `interaction_score_object` object. No possible differential analysis.")
      return()
    }
  }

  if(is.null(group_dn)){
    warning("Please provide the `group_dn` parameter which is a vector of sample names found in the colnames of the `interaction_score` data.frame of the `interaction_score_object` object.")
    return()
  } else {
    group_dn <- group_dn[!is.na(which(group_dn %in% colnames(interaction_score)))]
    if(length(group_dn)==0){
      warning("The values in `group_dn` do not intersect with the colnames of the `interaction_score` data.frame of the `interaction_score_object` object.")
      return()
    }
    if(length(group_dn)==1){
      warning("Only one value in `group_dn` found to intersect with the colnames of the `interaction_score` data.frame of the `interaction_score_object` object. No possible differential analysis.")
      return()
    }
  }

  if(is.null(label_up)){label_up<-"up"}
  if(is.null(label_dn)){label_up<-"down"}

  statistics <- match.arg(statistics)
  cat("The selected statistics mode is:", statistics, "\n")

  # Load or check format of lr_pathways_to_plot
  if(is.null(lr_pathways_to_plot)){
    #data("LR_main_pathways", package = "scAvengeR")
    lr_pathways_to_plot <- scAvengeR::lr_pathways_to_plot
    if(is.null(color_pathways)){
      color_pathways <- c("#D52626","#320EE6","#1B7727","#E6890E","#169A9A")
    } else {
      if(length(color_pathways) != length(lr_pathways_to_plot)){
        warning(paste0("The `color_pathways` should be a color vector of size ",length(lr_pathways_to_plot)))
        return()
      }
    }
  } else {
    if (!is.list(lr_pathways_to_plot)) stop("'lr_pathways_to_plot' must be a list.")
    for (i in seq_along(lr_pathways_to_plot))
      if (!is.data.frame(lr_pathways_to_plot[[i]]) ||
          !all(c("receptor", "ligand") %in% colnames(lr_pathways_to_plot[[i]])))
        stop(paste("Element", i, "of 'lr_pathways_to_plot' must be a data.frame with 'receptor' and 'ligand' columns."))

    if(is.null(color_pathways)){
      # Generate a color vector of size length(lr_pathways_to_plot)
      color_pathways <- brewer.pal(min(length(lr_pathways_to_plot), 12), "Set3") # Up to 12 distinct colors
      # If more than 12 elements, generate more colors using colorRampPalette
      if (length(lr_pathways_to_plot) > 12) {
        color_pathways <- colorRampPalette(brewer.pal(12, "Set3"))(length(lr_pathways_to_plot))
      }
    } else {
      if(length(color_pathways) != length(lr_pathways_to_plot)){
        warning(paste0("The `color_pathways` should be a color vector of size ",length(lr_pathways_to_plot)))
        return()
      }
    }
  }

  # Check for cell_type_groups
  if(is.null(cell_type_groups)){
    warning("Please provide the `cell_type_groups` parameter which is a named list containing the vectors of cell type subgroups")
    return()
  } else {
    detect_invalid <- match(as.character(unlist(cell_type_groups)),colnames(interaction_score_object$proportion_cell_type))
    if(length(which(is.na(detect_invalid))) > 0){
      warning(paste0("Some cell types were not found in `the interaction_score_object` : "),paste0(as.character(unlist(cell_type_groups))[which(is.na(detect_invalid))],collapse=", "))
      return()
    }
  }
  #########################################################################

  ## Compute statistics
  #########################################################################
  idx_up <- match(group_up,colnames(interaction_score))
  idx_dn <- match(group_dn,colnames(interaction_score))

  if(statistics == "ttest"){

    perform_test <- function(row){
      group1 <- as.numeric(row[idx_up])
      group2 <- as.numeric(row[idx_dn])
      t <- t.test(group1, group2)
      return(c(t$statistic,t$p.value))
    }
    message("Computing statistical tests")
    DE_LR <- as.data.frame(t(pbapply(interaction_score, 1, perform_test)))
    colnames(DE_LR) <- c("directionality","p.value")
    poss_comb$directionality <- NA
    poss_comb$p.value <- NA
    poss_comb$directionality[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$directionality
    poss_comb$p.value[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$p.value
  }

  if(statistics == "wilcoxon"){
    perform_test <- function(row){
      group1 <- as.numeric(row[idx_up])
      group2 <- as.numeric(row[idx_dn])

      w <- suppressWarnings(wilcox.test(group1, group2))
      wilcox_test <- w$p.value
      medianx <- median(group1)
      mediany <- median(group2)
      wilcox_dir <- medianx - mediany

      return(c(wilcox_dir,w$p.value))
    }
    message("Computing statistical tests")
    DE_LR <- as.data.frame(t(pbapply(interaction_score, 1, perform_test)))
    colnames(DE_LR) <- c("directionality","p.value")
    if(length(which(is.nan(DE_LR$p.value))) > 0){
      DE_LR$directionality[which(is.nan(DE_LR$p.value))] <- NaN
    }
    poss_comb$directionality <- NA
    poss_comb$p.value <- NA
    poss_comb$directionality[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$directionality
    poss_comb$p.value[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$p.value
  }

  if(statistics == "lm"){
    # Then perform differential expression
    design <- data.frame("Sample"=colnames(interaction_score),
                         "DEG"=NA)
    design$DEG[which(design$Sample %in% group_up)]<-"up"
    design$DEG[which(design$Sample %in% group_dn)]<-"dn"
    design$DEG <- factor(design$DEG,levels=c("dn","up"))
    design <- model.matrix(~ DEG, data = design)

    idx_zero <- which(apply(interaction_score,1,sum)==0)
    #if(length(idx_zero) > 0){interaction_score <- interaction_score[-idx_zero,]}
    message("Computing statistical tests")
    DE_LR <- eBayes(lmFit(interaction_score[-idx_zero,], design))
    DE_LR <- topTable(DE_LR, coef = colnames(design)[ncol(design)], adjust = "fdr",number=nrow(interaction_score))

    poss_comb$directionality <- NA
    poss_comb$p.value <- NA
    poss_comb$directionality[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$logFC
    poss_comb$p.value[match(rownames(DE_LR),poss_comb$unique_id)] <- DE_LR$P.Value
    poss_comb[is.na(poss_comb)] <- NaN
  }
  DE_LR <- poss_comb
  DE_LR_up <- DE_LR[which((DE_LR$directionality > 0)&(DE_LR$p.value < 0.05)),]
  DE_LR_dn <- DE_LR[which((DE_LR$directionality < 0)&(DE_LR$p.value < 0.05)),]
  #########################################################################

  ## Make Barplot of Total and Significant Interactions
  #########################################################################
  # Chi-square test
  dat <- data.frame(
    "Sig" = c(nrow(DE_LR_up), nrow(DE_LR_dn)),
    "Non_Sig" = c((nrow(DE_LR)-nrow(DE_LR_up)), (nrow(DE_LR)-nrow(DE_LR_dn))),
    row.names = c(label_up, label_dn),
    stringsAsFactors = FALSE
  )
  dat <- as.matrix(t(dat))
  chisquare <- chisq.test(dat)

  # Convert to a long data frame
  df <- as.data.frame(as.table(dat))
  colnames(df) <- c("Significance", "Condition", "Count")

  df$Significance <- factor(df$Significance, levels=c("Non_Sig","Sig"))
  df$Condition <- factor(df$Condition, levels= c(label_dn,label_up))

  # Subset for “Significant” only – must match your row label in ‘dat’
  df_sig <- subset(df, Significance == "Sig")

  # First (stacked) barplot: total interactions
  p1 <- ggplot(df, aes(x=Condition, y=Count, fill=Significance)) +
    geom_col(color="black", position="stack", linewidth = 0.2) +
    scale_fill_manual(values=c("gray","red")) +
    labs(title="# total interactions",
         y="Number of total interactions",
         x=NULL) +
    theme_minimal(base_size=10) +  # Make labels smaller
    theme(legend.title=element_text(size=9),
          legend.text=element_text(size=8))

  # Second (single-color) barplot: significant only
  p2 <- ggplot(df_sig, aes(x=Condition, y=Count)) +
    geom_col(fill="red", color="black",linewidth = 0.2) +
    labs(title="# significant interactions",
         y="Number of significant interactions",
         x=NULL) +
    theme_minimal(base_size=10) +
    annotate(
      "text",
      x=1.2,
      y=max(df_sig$Count)*0.9,
      label=paste0("pval chi-sq: ",
                   format(signif(chisquare$p.value, 2), scientific=TRUE)),
      size=3
    )

  # Combine side-by-side
  total_interaction_plot <- p1 + p2
  #########################################################################

  ## Circos for pathways of interest
  #########################################################################
  cell_types <- as.character(unlist(cell_type_groups))

  # Make a function for plotting the circos and call it for up and down groups
  plot_circos <- function(group){

    # Parameters for circos plot
    LR_types <- names(lr_pathways_to_plot)
    alpha <- 1
    plot_intensity <- T
    generate_ranges <- function(n, step = 2) {
      starts <- seq(1, by = step, length.out = n)
      ends <- starts + step
      Map(c, starts, ends)
    }
    range <- generate_ranges(length(lr_pathways_to_plot))
    pval <- 0.05
    collect_sig_genes <- data.frame()
    collect_links <- data.frame()

    # Define Sectors
    sector <- data.frame("seq"=rep(seq(1,max(unlist(range)),1),length(cell_types)),
                         "sector"=rep(cell_types,each=max(unlist(range))))
    sector$sector <- factor(sector$sector,levels=cell_types)

    # Define space between sectors
    idx_space <- c(rep(2,length(cell_types)))
    idx_space[cumsum(lapply(cell_type_groups,length))] <- 8

    circos.clear()
    circos.par("track.height" = 0.1, gap.after =idx_space,
               track.margin = c(0.2,0.2))
    #plot.new()
    # Circos initialization and plotting
    circos.initialize(sector$sector,x = sector$seq)

    group_vector <- as.character(unlist(mapply(rep, names(cell_type_groups), lengths(cell_type_groups))))
    color_mapping <- setNames(qualitative_hcl(length(unique(group_vector)), palette = "Pastel1"), unique(group_vector))
    group_colors <- color_mapping[group_vector]
    circos.track(sector$sector, ylim = c(0, 1),
                 track.height = 0.05,        # Controls thickness of the track
                 bg.col = group_colors, # Box color
                 bg.border = "black")

    suppressMessages(
      for(cat in unique(sector$sector)){
        circos.text(
          x = 9,  # Increase this value if labels still touch boxes
          y = 3,  # Adjust this to move labels outward
          labels = cat,
          facing = "clockwise",
          niceFacing = TRUE,
          adj = c(0.1, 0),
          cex = 0.4,
          col = "black",
          font = 1,
          sector.index = cat
        )
      }
    )

    title(paste0(group," LR sig interations"), cex.main = 1, font.main = 1.5)

    melt_Sig_Int_sum <- data.frame()

    for(type in 1:length(LR_types)){

      DE_LR_pathway <- DE_LR[which((DE_LR$celltype1 %in% cell_types)&(DE_LR$celltype2 %in% cell_types)),]
      lr_pathways_to_plot[[type]]$pair <- paste(lr_pathways_to_plot[[type]]$receptor,lr_pathways_to_plot[[type]]$ligand,sep="_")
      DE_LR_pathway$comb1<-paste(DE_LR_pathway$gene1,DE_LR_pathway$gene2,sep="_")
      DE_LR_pathway$comb2<-paste(DE_LR_pathway$gene2,DE_LR_pathway$gene1,sep="_")

      idx1 <- which(DE_LR_pathway$comb1 %in% lr_pathways_to_plot[[type]]$pair)
      idx2 <- which(DE_LR_pathway$comb2 %in% lr_pathways_to_plot[[type]]$pair)
      DE_LR_pathway <- DE_LR_pathway[unique(c(idx1,idx2)),]
      DE_LR_pathway <- DE_LR_pathway[complete.cases(DE_LR_pathway),]
      DE_LR_pathway$comb1 <- NULL
      DE_LR_pathway$comb2 <- NULL

      # Collect number of significant interactions in both groups
      comb <- as.data.frame(combn(unique(as.character(unique(c(DE_LR_pathway$celltype1,DE_LR_pathway$celltype2)))),2))
      Sig_Int <- as.data.frame(matrix(nrow=2,ncol=ncol(comb)))
      rownames(Sig_Int) <- c("group_up","group_dn")

      for(i in 1:ncol(comb)){

        idx1 <- which((DE_LR_pathway$celltype1==comb[1,i])&(DE_LR_pathway$celltype2==comb[2,i]))
        idx2 <- which((DE_LR_pathway$celltype1==comb[2,i])&(DE_LR_pathway$celltype2==comb[1,i]))
        poss_comb_pos <- DE_LR_pathway[c(idx1,idx2),]
        poss_comb_pos <- poss_comb_pos[which((poss_comb_pos$p.value < pval)),]

        Sig_Int[1,i] <- nrow(poss_comb_pos[which(poss_comb_pos$directionality>0),])
        Sig_Int[2,i] <- nrow(poss_comb_pos[which(poss_comb_pos$directionality<0),])

        colnames(Sig_Int)[i] <- paste0(comb[1,i],":",comb[2,i])

        collect_sig_genes <- rbind(collect_sig_genes,
                                   rbind(poss_comb_pos[which(poss_comb_pos$directionality > 0),],
                                         poss_comb_pos[which(poss_comb_pos$directionality < 0),]))

      }

      Sig_Int$Groups<-rownames(Sig_Int)

      melt_Sig_Int <- reshape2::melt(Sig_Int,id.vars = "Groups")

      # Reorder LR interactions. It's important for putting the high intensity link in the front
      melt_Sig_Int$variable <- factor(melt_Sig_Int$variable,
                                      levels = names(apply(Sig_Int[,-ncol(Sig_Int)],2,sum)[order(apply(Sig_Int[,-ncol(Sig_Int)],2,sum),decreasing=T)]))

      # Filter out LR interaction having 0 interactions for both up and down groups
      toremove<-c()
      for(i in 1:length(unique(melt_Sig_Int$variable))){
        if(sum(melt_Sig_Int$value[which(melt_Sig_Int$variable %in% as.character(unique(melt_Sig_Int$variable)[i]))]) == 0){
          toremove <- c(toremove,as.character(unique(melt_Sig_Int$variable)[i]))
        }
      }
      if(length(toremove) > 0){melt_Sig_Int <- melt_Sig_Int[-which(melt_Sig_Int$variable %in% toremove),]}

      # Circos
      melt_Sig_Int$cell_type1 <- strsplit2(melt_Sig_Int$variable,split=":")[,1]
      melt_Sig_Int$cell_type2 <- strsplit2(melt_Sig_Int$variable,split=":")[,2]
      melt_Sig_Int$Pathway <- LR_types[type]

      melt_Sig_Int_sum <- rbind(melt_Sig_Int_sum,melt_Sig_Int)

      if(plot_intensity==T){
        melt_Sig_Int$color<-colorRampPalette(c("white",color_pathways[type]))(101)[round(scales::rescale(melt_Sig_Int$value,to=c(0,100)))+1]
      }
      if(plot_intensity==F){
        melt_Sig_Int$color<-color_pathways[type]
      }

      circos.par(track.margin=c(-0.18,-0.18))

      # Select Appropriate Group
      tmp <- melt_Sig_Int[which(melt_Sig_Int$Groups==group),]
      tmp <- tmp[order(tmp$value,decreasing=F),]

      for(i in 1:nrow(tmp)){
        if(tmp$value[i] != 0){
          circos.link(tmp$cell_type1[i], range[[type]],
                      tmp$cell_type2[i], range[[type]],
                      col = alpha(tmp$color[i],alpha),lwd=0.01,border = "black",
                      h.ratio = sample(seq(0.2,0.7,0.1),size = 1))
        }
      }
      collect_links <- rbind(collect_links,melt_Sig_Int)
    }

    # Format collect_links
    collect_links <- collect_links[which(collect_links$Groups == group),]
    colnames(collect_links)[which(colnames(collect_links)=="variable")] <- "CellType_Pair"
    collect_links$color <- NULL
    colnames(collect_links)[which(colnames(collect_links)=="value")] <- "Interaction_Number"

    # Capture the plot and return
    p <- recordPlot()
    return(list("circos_plot"=p,"Number_Significant_Interactions"=collect_links))
  }

  circos_up <- plot_circos("group_up")
  circos_dn <- plot_circos("group_dn")
  Significant_Interaction_Number <- rbind(circos_up$Number_Significant_Interactions,
                                          circos_dn$Number_Significant_Interactions)
  #########################################################################

  ## Barplot Interaction Main Cell Type
  #########################################################################
  Significant_Interaction_Number$cell_type1_main <- result <- sapply(Significant_Interaction_Number$cell_type1,
                                                                     function(x) {
                                                                       names(cell_type_groups)[sapply(cell_type_groups, function(y) x %in% y)]})

  Significant_Interaction_Number$cell_type2_main <- result <- sapply(Significant_Interaction_Number$cell_type2,
                                                                     function(x) {
                                                                       names(cell_type_groups)[sapply(cell_type_groups, function(y) x %in% y)]})
  Significant_Interaction_Number$CellType_Main_Pair <- paste(Significant_Interaction_Number$cell_type1_main,
                                                             Significant_Interaction_Number$cell_type2_main,
                                                             sep=":")

  # For Fine Cell Type
  ordered_pairs <- Significant_Interaction_Number %>%
    group_by(CellType_Pair) %>%
    summarise(avg_interaction = mean(Interaction_Number, na.rm = TRUE)) %>%
    arrange(desc(avg_interaction)) %>%
    pull(CellType_Pair)
  Significant_Interaction_Number$CellType_Pair <- factor(Significant_Interaction_Number$CellType_Pair,levels=ordered_pairs)

  Barplot_Fine_Per_Pathway <- ggplot(data=Significant_Interaction_Number, aes(x=CellType_Pair, y=Interaction_Number, fill=Pathway)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    ylab("Significant Interaction Number") + ggtitle("Number of Significant Interactions") +
    scale_fill_manual(values=color_pathways) + facet_wrap(~Groups+Pathway,nrow=2) +
    theme(legend.position = "right",
          title =element_text(size=6, face='bold'),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=1,angle=90),axis.title.y = element_text(size=8),
          axis.title.x.bottom = element_blank())

  # For Main Cell Type
  ordered_pairs <- Significant_Interaction_Number %>%
    group_by(CellType_Main_Pair) %>%
    summarise(avg_interaction = mean(Interaction_Number, na.rm = TRUE)) %>%
    arrange(desc(avg_interaction)) %>%
    pull(CellType_Main_Pair)
  Significant_Interaction_Number$CellType_Main_Pair <- factor(Significant_Interaction_Number$CellType_Main_Pair,levels=ordered_pairs)

  Barplot_Main_Per_Pathway <- ggplot(data=Significant_Interaction_Number, aes(x=CellType_Main_Pair, y=Interaction_Number, fill=Pathway)) +
    geom_bar(stat="identity", color="black", position=position_dodge()) +
    ylab("Significant Interaction Number") + ggtitle("Number of Significant Interactions") +
    scale_fill_manual(values=color_pathways) + facet_wrap(~Groups+Pathway,nrow=2) +
    theme(legend.position = "right",
          title =element_text(size=6, face='bold'),
          axis.text.y = element_text(size=8),
          axis.text.x = element_text(size=4,angle=90),axis.title.y = element_text(size=8),
          axis.title.x.bottom = element_blank())
  #########################################################################

  # Heatmap Interaction for the five pathways
  #########################################################################
  template <- as.data.frame(matrix(nrow=length(cell_types),ncol=length(cell_types),0))
  rownames(template) <- cell_types
  colnames(template) <- cell_types

  groups <- c("group_up","group_dn")
  for(group in groups){
    tmp <- Significant_Interaction_Number[which(Significant_Interaction_Number$Groups == group),]

    for(i in 1:nrow(template)){
      for(j in 1:ncol(template)){
        idx1 <- which((tmp$cell_type1 == rownames(template)[i])&(tmp$cell_type2 == colnames(template)[j]))
        idx2 <- which((tmp$cell_type1 == colnames(template)[j])&(tmp$cell_type2 == rownames(template)[i]))
        idx <- c(idx1,idx2)
        template[i,j] <- sum(tmp$Interaction_Number[idx])

      }
    }
    assign(paste0("heat_",group),template)
    rm(tmp)
  }

  breaks <- seq(0,max(c(max(heat_group_dn[cell_types,cell_types]),max(heat_group_up[cell_types,cell_types]))))

  for(i in 1:nrow(heat_group_up)){
    for(j in i:ncol(heat_group_up)){
      heat_group_dn[i,j]<-heat_group_up[i,j]
    }
  }
  for(i in 1:nrow(heat_group_dn)){heat_group_dn[i,i]<-NA}

  Heat_Pathway <- ComplexHeatmap::pheatmap(as.matrix(heat_group_dn),cluster_rows = F,cluster_cols = F,breaks = breaks,
                           color = colorRampPalette(c("white","darkblue","orange","red"))(length(breaks)),
                           main="Interaction Pathways - group_up (top-right) versus group_dn (bottom left)", scale="none",border_color="black",
                           treeheight_col =4,treeheight_row = 4,
                           show_rownames=T,show_colnames=T,
                           fontsize_row=7,fontsize_col = 7,
                           cellwidth = 8,cellheight = 8,
                           gaps_row = as.numeric(cumsum(lapply(cell_type_groups,length))),
                           gaps_col = as.numeric(cumsum(lapply(cell_type_groups,length))),
                           name="# Interactions")
  #########################################################################

  # Heatmap Interaction all LRs all subsets
  #########################################################################
  template <- as.data.frame(matrix(nrow=length(cell_types),ncol=length(cell_types),0))
  rownames(template) <- cell_types
  colnames(template) <- cell_types

  groups <- c("group_up","group_dn")
  for(group in groups){
    if(group == "group_up"){
      tmp <- DE_LR[which((DE_LR$directionality > 0) & (DE_LR$p.value < 0.05)),]
    }
    if(group == "group_dn"){
      tmp <- DE_LR[which((DE_LR$directionality < 0) & (DE_LR$p.value < 0.05)),]
    }

    for(i in 1:nrow(template)){
      for(j in 1:ncol(template)){
        idx1 <- which((tmp$celltype1 == rownames(template)[i])&(tmp$celltype2 == colnames(template)[j]))
        idx2 <- which((tmp$celltype1 == colnames(template)[j])&(tmp$celltype2 == rownames(template)[i]))
        idx <- c(idx1,idx2)
        template[i,j] <- length(idx)
      }
    }
    assign(paste0("heat_",group),template)
    rm(tmp)
  }

  breaks <- seq(0,max(c(max(heat_group_dn[cell_types,cell_types]),max(heat_group_up[cell_types,cell_types]))))

  for(i in 1:nrow(heat_group_up)){
    for(j in i:ncol(heat_group_up)){
      heat_group_dn[i,j]<-heat_group_up[i,j]
    }
  }
  for(i in 1:nrow(heat_group_dn)){heat_group_dn[i,i]<-NA}

  Heat_All_LR <- ComplexHeatmap::pheatmap(as.matrix(heat_group_dn),cluster_rows = F,cluster_cols = F,breaks = breaks,
                                           color = colorRampPalette(c("white","darkblue","orange","red"))(length(breaks)),
                                           main="Interaction All LR - group_up (top-right) versus group_dn (bottom left)", scale="none",border_color="black",
                                           treeheight_col =4,treeheight_row = 4,
                                           show_rownames=T,show_colnames=T,
                                           fontsize_row=7,fontsize_col = 7,
                                           cellwidth = 8,cellheight = 8,
                                           gaps_row = as.numeric(cumsum(lapply(cell_type_groups,length))),
                                           gaps_col = as.numeric(cumsum(lapply(cell_type_groups,length))),
                                           name="# Interactions")
  #########################################################################

  returned_object <- list("Differential_Analysis" = DE_LR,
                          "Barplot_Total_Interactions" = total_interaction_plot,
                          "circos_group_up" = circos_up$circos_plot,
                          "circos_group_dn" = circos_dn$circos_plot,
                          "Significant_Interaction_Number" = Significant_Interaction_Number,
                          "Barplot_Interaction_Pathway_Fine"= Barplot_Fine_Per_Pathway,
                          "Barplot_Interaction_Pathway_Main"= Barplot_Main_Per_Pathway,
                          "Heatmap_Interacion_Pathway" = Heat_Pathway,
                          "Heatmap_Interacion_All_LR" = Heat_All_LR)


  return(returned_object)
}

