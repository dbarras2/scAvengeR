#' @export
CellType_Proportion_Heatmap <- function(single_cell_data,
                                        sample_colname,
                                        subset_data = NULL,
                                        stratification,
                                        group_up,
                                        group_dn,
                                        group_up_col = NULL,
                                        group_dn_col = NULL,
                                        statistics = "wilcoxon",
                                        include_ratio = TRUE,
                                        stratification_names = NULL,
                                        pval_range = NULL,
                                        stratification_for_ratio = NULL){
  require(ComplexHeatmap)
  require(limma)
  require(reshape2)
  require(crayon)
  require(circlize)
  require(usethis)
  require(scales)
  require(dplyr)

  meta <- single_cell_data

  if(is.null(meta)){
    warning("Please provide a single-cell metadata object with cell type annotation under the `single_cell_data` parameter as explained in the vignette")
    return()
  }
  if(is.null(sample_colname)){
    warning("Please provide the name of the column containing the sample identifier in the single-cell metadata object under the `sample_colname` parameter")
    return()
  }
  if(is.na(match(statistics,c("wilcoxon","ttest")))){
    warning("The selected statistics should be `wilcoxon` or `ttest`")
    return()
  }

  ## Subset Data
  #########################################################################
  if(!is.null(subset_data)){
    if (!is.list(subset_data)) {
      warning("The `subset_data` parameter should be a list")
      return()
    }
    if(length(which(names(subset_data)  %in% colnames(meta) ==F))>0){
      warning("The selected subset_data parameters were not found in col.names of the single-cell object")
      return()
    }
    idx_to_keep <- 1:nrow(meta)
    for(l in 1:length(subset_data)){
      idx_to_keep <- intersect(idx_to_keep,which(meta[,names(subset_data)[l]] %in% subset_data[[l]]))
    }
    meta <- meta[idx_to_keep,]
  }
  #########################################################################

  ## Define categorizing factor and groups
  #########################################################################
  # Check that all parameters were filled correctly and are available
  if(1){
    if(is.null(stratification)){
      warning("Please provide the name of the column containing the stratification group that will be compared in the single-cell metadata object under the `stratification` parameter")
      return()
    }
    if(length(which(stratification  %in% colnames(meta) ==F))>0){
      warning("The selected stratification name parameter was not found in col.names of the single-cell object")
      return()
    }
    if(is.null(group_up)){
      warning("Please provide the stratification factor group down `group_up` that will be used for comparison with `group_dn` contained in the `stratification` column")
      return()
    }
    idx<-which(meta[,stratification]==group_up)
    if(length(idx) == 0){
      warning(paste0("Stratification factor not found for `group_up`; available ones are: ",paste(unique(meta[,stratification]),collapse=", ")
      ))
      return()
    }
    if(is.null(group_dn)){
      warning("Please provide the stratification factor group down `group_dn` that will be used for comparison with `group_up` contained in the `stratification` column")
      return()
    }
    idx<-which(meta[,stratification]==group_dn)
    if(length(idx) == 0){
      warning(paste0("Stratification factor not found for `group_dn`; available ones are: ",paste(unique(meta[,stratification]),collapse=", ")
      ))
      return()
    }
  }

  # Creating Sample annotation file
  Sample_annot <- meta %>%
    group_by(across(sample_colname)) %>%
    summarise(across(all_of(stratification), unique))

  Sample_annot <- as.data.frame(Sample_annot)
  Sample_annot$Stratif <- Sample_annot[,stratification]

  idx_up <- which(Sample_annot$Stratif %in% group_up)
  idx_dn <- which(Sample_annot$Stratif %in% group_dn)
  if(length(intersect(idx_up,idx_dn)) > 0){
    warning("Group.up and Group_dn are overlapping - Please provide non-overlapping categories")
    return()
  }

  Sample_annot$Stratif[idx_up] <- "up"
  Sample_annot$Stratif[idx_dn] <- "dn"
  Sample_annot <- Sample_annot[which(Sample_annot$Stratif %in% c("up","dn")),]
  Sample_annot$Stratif <- factor(Sample_annot$Stratif, levels = c("up","dn"))

  # Define patient levels
  sample_levels <- levels(factor(as.character(Sample_annot[,sample_colname])))
  meta[,sample_colname] <- factor(meta[,sample_colname], levels = sample_levels)

  # Define sample order
  sample_order_1 <- Sample_annot[which(Sample_annot$Stratif=="up"),sample_colname]
  sample_order_2 <- Sample_annot[which(Sample_annot$Stratif=="dn"),sample_colname]
  #########################################################################

  # List and stratification of all cell types
  #########################################################################
  order_strat <- order(as.numeric(gsub("Cell_Type_Strat","",colnames(meta)[grep("Cell_Type_Strat",colnames(meta))])))
  cell_type_strats <- colnames(meta)[grep("Cell_Type_Strat",colnames(meta))][order_strat]

  if (!is.factor(meta[,cell_type_strats[length(cell_type_strats)]])) {
    warning(paste0("Please provide the last cell type stratification (",cell_type_strats[length(cell_type_strats)],") as factor column and define the order of cell types (how it will appear in the heatmap) with levels using factor(x,levels=...)"))
    return()
  }

  column_cell_type <- cell_type_strats[length(cell_type_strats)]

  # Vector of granular cell types
  cell_subtypes <- levels(droplevels(meta[,cell_type_strats[length(cell_type_strats)]]))

  # Vector of bulked categories
  cell_types <- c()
  for(strata in rev(cell_type_strats[-length(cell_type_strats)])){
    tmp <- meta[order(match(meta[,column_cell_type],cell_subtypes)),strata]
    for(str in unique(tmp)){
      if(length(unique(meta[,column_cell_type][which(meta[,strata] == str)]))>1){
        cell_types <- c(cell_types,str)
      }
    }
  }
  if(length(which(duplicated(cell_types)))>0){
    cell_types <- cell_types[-which(duplicated(cell_types))]
  }

  # Rename the bulk cell types (the parent cell type that have more than one daughter cell type)
  meta_backup <- meta
  cell_types_bulk <- paste0(cell_types,"_bulk")
  for(strata in rev(cell_type_strats[-length(cell_type_strats)])){
    tmp <- meta[,strata]
    idx <- which(!is.na(match(tmp,cell_types)))
    tmp[idx] <- cell_types_bulk[match(tmp,cell_types)[idx]]
    meta[,strata] <- tmp
  }
  cell_types <- cell_types_bulk

  # Layers of stratification
  row_norm <- c("All",cell_type_strats[-length(cell_type_strats)])

  collect_stats <- as.data.frame(matrix(nrow=length(row_norm),ncol=length(c(cell_subtypes,cell_types)),NA))
  rownames(collect_stats) <- row_norm
  colnames(collect_stats) <- c(cell_subtypes,cell_types)

  collect_stats_ratio <- as.data.frame(matrix(nrow=length(c(cell_subtypes,cell_types)),
                                              ncol=length(c(cell_subtypes,cell_types)),
                                              NA))
  rownames(collect_stats_ratio) <- c(cell_subtypes,cell_types)
  colnames(collect_stats_ratio) <- c(cell_subtypes,cell_types)
  #########################################################################

  # Decide Statistical Test
  ########################################################################
  makeStatisticalTest <- function(Imm1,Imm2,test){
    result<-list()
    if((length(Imm1[!is.na(Imm1)])<1)|(length(Imm2[!is.na(Imm2)])<1)){
      result$statistics<-NA
      result$p.value<-NA
      return(result)
    }
    if(test=="ttest"){
      t<-suppressWarnings(t.test(Imm1,Imm2))
      result$statistics<-t$statistic
      result$p.value<-t$p.value
      return(result)
    }
    if(test=="wilcoxon"){
      t<-suppressWarnings(wilcox.test(Imm1,Imm2))
      result$p.value<-t$p.value
      val<-ifelse((median(Imm1)-median(Imm2))==0,(mean(Imm1)-mean(Imm2)),(median(Imm1)-median(Imm2)))
      result$statistics<-val
      return(result)
    }
  }
  ########################################################################

  # Proportions out of total cells
  ########################################################################
  Proportions <- prop.table(table(meta[,which(colnames(meta) %in% c(sample_colname,column_cell_type))]),1)
  Proportions <- reshape2::melt(Proportions)
  Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
  colnames(Proportions)[2] <- "Cell_Type"

  for(i in 1:length(cell_subtypes)){
    tmp1 <- Proportions[which(Proportions$Cell_Type==cell_subtypes[i]),]
    if(nrow(tmp1)>0){
      Imm1 <- tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1[,sample_colname][which(tmp1$Stratif=="up")])
      Imm2 <- tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1[,sample_colname][which(tmp1$Stratif=="dn")])


      t <- makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats["All",cell_subtypes[i]] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
    }
  }

  strats <- cell_type_strats[-length(cell_type_strats)]
  for(strat in strats){
    Proportions <- prop.table(table(meta[,which(colnames(meta) %in% c(sample_colname,strat))]),1)
    Proportions <- reshape2::melt(Proportions)
    Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
    colnames(Proportions)[2] <- "Cell_Type"

    cts <- as.character(unique(Proportions$Cell_Type))
    cts <- cts[which(cts %in% colnames(collect_stats))]
    for(ct in cts){
      tmp1 <- Proportions[which(Proportions$Cell_Type==ct),]
      Imm1 <- tmp1$value[which(tmp1$Stratif=="up")];names(Imm1)<-as.character(tmp1[,sample_colname][which(tmp1$Stratif=="up")])
      Imm2 <- tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2)<-as.character(tmp1[,sample_colname][which(tmp1$Stratif=="dn")])

      t <- makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
      collect_stats["All",ct] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
    }
  }
  ########################################################################

  # Proportions out of stratification layers
  ########################################################################

  # Interrogate Fine cell subtypes out of all strats
  strats <- cell_type_strats[-length(cell_type_strats)]
  for(strat in strats){
    outofs <- unique(meta[,strat])
    for(outof in outofs){
      idx_outof <- which(meta[,strat] == outof)

      Proportions <- prop.table(table(meta[idx_outof,which(colnames(meta) %in% c(sample_colname,column_cell_type))]),1)
      Proportions <- reshape2::melt(Proportions)
      Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions)[2] <- "Cell_Type"

      cts <- as.character(unique(meta[idx_outof,column_cell_type]))
      cts <- cts[which(cts %in% colnames(collect_stats))]
      for(ct in cts){
        tmp1 <- Proportions[which(Proportions$Cell_Type==ct),]
        Imm1 <- tmp1$value[which(tmp1$Stratif=="up")];names(Imm1) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="up")])
        Imm2 <- tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="dn")])

        t <- makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
        collect_stats[strat,ct] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
      }
    }
  }

  # Interrogate Bulk populations in lower strats
  strats <- cell_type_strats[-length(cell_type_strats)]
  if(length(strats) > 0){
    for(s in 1:length(strats)){
      strat <- strats[s]
      for(o in 1:(s-1)){
        outof_strat <- strats[o]

        outofs <- unique(meta[,outof_strat])
        for(outof in outofs){
          idx_outof <- which(meta[,outof_strat] == outof)

          Proportions <- prop.table(table(meta[idx_outof,which(colnames(meta) %in% c(sample_colname,strat))]),1)
          Proportions <- reshape2::melt(Proportions)
          Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
          colnames(Proportions)[2] <- "Cell_Type"

          cts <- as.character(unique(meta[idx_outof,strat]))
          cts <- cts[which(cts %in% colnames(collect_stats))]
          for(ct in cts){
            tmp1 <- Proportions[which(Proportions$Cell_Type==ct),]
            Imm1 <- tmp1$value[which(tmp1$Stratif=="up")];names(Imm1) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="up")])
            Imm2 <- tmp1$value[which(tmp1$Stratif=="dn")];names(Imm2) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="dn")])

            t <- makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
            collect_stats[outof_strat,ct] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
          }
        }
      }
    }
  }
  ########################################################################

  # Ratio of cell types out of total cells
  ########################################################################

  # Subset Data For Ratio Computation
  if(!is.null(stratification_for_ratio)){
    if(!is.list(stratification_for_ratio)) {
      warning("The `stratification_for_ratio` parameter should be a list")
      return()
    }
    if(length(which(names(stratification_for_ratio)  %in% colnames(meta) ==F))>0){
      warning("The selected stratification_for_ratio parameters were not found in col.names of the single-cell object")
      return()
    }
    idx_to_keep <- 1:nrow(meta_backup)
    for(l in 1:length(stratification_for_ratio)){
      idx_to_keep <- intersect(idx_to_keep,which(meta_backup[,names(stratification_for_ratio)[l]] %in% stratification_for_ratio[[l]]))
    }
    meta_ratio <- meta[idx_to_keep,]
  } else {
    meta_ratio <- meta
    }

  # Compute Ratio for all Cell Type at all stratification and for all patients
  strats <- cell_type_strats
  for(s in 1:length(strats)){
    strat <- strats[s]
    if(s==1){
      Proportions <- prop.table(table(meta_ratio[,which(colnames(meta_ratio) %in% c(sample_colname,strat))]),1)
      Proportions <- reshape2::melt(Proportions)
      Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions)[2] <- "Cell_Type"
    } else {
      Proportions_tmp <- prop.table(table(meta_ratio[,which(colnames(meta_ratio) %in% c(sample_colname,strat))]),1)
      Proportions_tmp <- reshape2::melt(Proportions_tmp)
      Proportions_tmp$Stratif <- Sample_annot$Stratif[match(Proportions_tmp[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions_tmp)[2] <- "Cell_Type"

      idx <- which(unique(Proportions_tmp$Cell_Type) %in% unique(Proportions$Cell_Type))
      if(length(idx) > 0){
        Proportions_tmp <- Proportions_tmp[-which(Proportions_tmp$Cell_Type %in% unique(Proportions_tmp$Cell_Type)[idx]),]
      }
      Proportions <- rbind(Proportions,Proportions_tmp)
      rm(Proportions_tmp)
    }
  }

  # Compute and collect statistics for Ratio comparisons
  for(ct_row in rownames(collect_stats_ratio)){
    for(ct_col in colnames(collect_stats_ratio)){

      tmp1<-Proportions[which(Proportions$Cell_Type==ct_col),]
      tmp2<-Proportions[which(Proportions$Cell_Type==ct_row),]
      tmp2<-tmp2[match(tmp1[,sample_colname],tmp2[,sample_colname]),]

      tmp1$Ratio<-log2(tmp1$value/tmp2$value)

      # Remove Inf values if
      if(statistics == "ttest"){
        warning("Infinite values for ratio will be removed as statistics has been set to `ttest`")
        if(length(which(tmp1$Ratio==(-Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(-Inf))]<-NA}
        if(length(which(tmp1$Ratio==(Inf)))>0){tmp1$Ratio[which(tmp1$Ratio==(Inf))]<-NA}
      }

      Imm1 <- tmp1$Ratio[which(tmp1$Stratif=="up")];names(Imm1) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="up")]);Imm1<-Imm1[complete.cases(Imm1)]
      Imm2 <- tmp1$Ratio[which(tmp1$Stratif=="dn")];names(Imm2) <- as.character(tmp1[,sample_colname][which(tmp1$Stratif=="dn")]);Imm2<-Imm2[complete.cases(Imm2)]

      if((length(Imm1) > 0)&(length(Imm2) > 0)){
        t <- makeStatisticalTest(Imm1[!is.na(Imm1)],Imm2[!is.na(Imm2)],statistics)
        collect_stats_ratio[ct_row,ct_col] <- ifelse(t$statistic>0,(-log10(t$p.value)),-(-log10(t$p.value)))
      }
    }
  }
  ########################################################################

  # Set up heatmap parameters
  ########################################################################
  if(length(row_norm) == 1){
    include_ratio <- TRUE
  }
  if(include_ratio){
    collect_stats <- rbind(collect_stats,collect_stats_ratio)
  }

  signif <- abs(collect_stats)
  signif_annot <- signif
  signif_annot[] <- lapply(signif_annot,as.character)

  for(i in 1:ncol(signif_annot)){
    signif_annot[,i] <- " "
    signif_annot[which(signif[,i] > (-log10(0.05))),i] <- "*"
    signif_annot[which(signif[,i] > (-log10(0.01))),i] <- "**"
    signif_annot[which(signif[,i] > (-log10(0.001))),i] <- "***"
  }

  lim <- max(abs(collect_stats),na.rm=T)

  if(!is.null(pval_range)){
    lim <- max(abs(pval_range),na.rm=T)
  }

  # Set up colors
  if(!is.null(group_up_col)){
    col_fun = colorRamp2(c(-lim, 0, lim), c("#3A34DF","white",group_up_col))
  }
  if(!is.null(group_dn_col)){
    col_fun = colorRamp2(c(-lim, 0, lim), c(group_dn_col,"white","#DF3434"))
  }
  if((!is.null(group_dn_col))&(!is.null(group_up_col))){
    col_fun = colorRamp2(c(-lim, 0, lim), c(group_dn_col,"white",group_up_col))
  }
  if((is.null(group_dn_col))&(is.null(group_up_col))){
    col_fun = colorRamp2(c(-lim, 0, lim), c("#3A34DF","white","#DF3434"))
  }
  ########################################################################

  if(length(row_norm) == 1){
    rownames(collect_stats)[1] <- "Out of total cells"
    heatmap <- Heatmap(as.matrix(collect_stats), name = "directed -log10(pval)", na_col = "lightgray",
                       col=col_fun,cluster_rows=F,
                       #row_order = 1:nrow(collect_stats),
                       column_order = 1:ncol(collect_stats),
                       column_title = "Cell Type Proportion and Ratio",
                       column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                       row_names_side = "left",column_names_side = "top",
                       row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),
                       row_labels=gsub("_"," ",rownames(collect_stats)),
                       row_split=c("A",rep("B",nrow(collect_stats)-1)),
                       row_gap = unit(3, "mm"),
                       border = TRUE,
                       cell_fun = function(j, i, x, y, w, h, col) {
                         grid.text(signif_annot[i, j], x, y,gp = gpar(fontsize = 8))
                       },
                       row_title = c("Proportion","Ratio between cell type proportion (out of total)"),
                       row_title_gp = gpar(fontsize = 12),
                       width = unit(14, "cm"), height = unit(16, "cm"))
  } else {

    # Create unique categories in alphabetical order
    Gap_template <- data.frame("CT"=colnames(collect_stats))
    Gap_template$strat <- meta[match(Gap_template$CT,meta[,column_cell_type]),cell_type_strats[length(cell_type_strats)-1]]
    idx <- which(is.na(Gap_template$strat))
    if(length(idx) > 0){
      Gap_template$strat[idx]<-"Bulk"
    }
    new_vector <- character(length(Gap_template$strat))

    current_label <- "B1"
    increment_label <- function(label) {
      num <- as.numeric(substring(label, 2))
      letter <- substring(label, 1, 1)

      if (num < 9) {
        return(paste0(letter, num + 1))
      } else {
        return(paste0(LETTERS[which(LETTERS == letter) + 1], 1))
      }
    }
    new_vector[1] <- current_label
    for(i in 2:length(Gap_template$strat)){
      if (Gap_template$strat[i] != Gap_template$strat[i - 1]) {
        current_label <- increment_label(current_label)
      }
      new_vector[i] <- current_label
    }
    Gap_template$alphabetical <- new_vector

    gaps_row <- c(rep("A",length(row_norm)),Gap_template$alphabetical)
    gaps_col <- gaps_row[-c(1:length(row_norm))]

    # Change names for stratification
    if(!is.null(stratification_names)){
      if(length(stratification_names) != length(row_norm)){
        warning("The vector of names for stratification layers (`stratification_names` parameter) should be the same length than the number of cell type stratification (= ",length(row_norm),")")
        return()
      }

      if(length(stratification_names) == 1){
        rownames(collect_stats)[1] <- "total cells"
      }
      if(length(stratification_names) > 1){
        names <- c("total cells",stratification_names[-length(stratification_names)])
        rownames(collect_stats)[1:length(names)] <- names
      }
    }
    rownames(collect_stats)[1:length(row_norm)] <- paste0("Out of ",rownames(collect_stats)[1:length(row_norm)])

    cat <- rle(Gap_template$strat)$values
    ########################################################################

    if(!include_ratio){
      heatmap <- Heatmap(as.matrix(collect_stats), name = "directed -log10(pval)", na_col = "lightgray",
                         col=col_fun,cluster_rows=F,
                         #row_order = 1:nrow(collect_stats),
                         column_order = 1:ncol(collect_stats),
                         column_title = "Cell Type Proportion",
                         column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                         row_names_side = "left",column_names_side = "top",
                         row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),
                         row_labels=gsub("_"," ",rownames(collect_stats)),
                         column_labels=gsub("_"," ",colnames(collect_stats)),
                         column_split=gaps_col,column_gap = unit(c(rep(1,length(unique(gaps_row))-1)), "mm"),
                         border = TRUE,
                         cell_fun = function(j, i, x, y, w, h, col) {
                           grid.text(signif_annot[i, j], x, y,gp = gpar(fontsize = 8))
                         },
                         row_title = c("Proportion"),
                         row_title_gp = gpar(fontsize = 12),
                         width = unit(14, "cm"), height = unit(1, "cm"),
                         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white"),
                                                                             labels = cat,
                                                                             labels_gp = gpar(col = "black", fontsize = 5))))

    }

    if(include_ratio){
      heatmap <- Heatmap(as.matrix(collect_stats), name = "directed -log10(pval)", na_col = "lightgray",
                         col=col_fun,cluster_rows=F,
                         #row_order = 1:nrow(collect_stats),
                         column_order = 1:ncol(collect_stats),
                         column_title = "Cell Type Proportion and Ratio",
                         column_title_gp = gpar(fontsize = 15, fontface = "bold"),
                         row_names_side = "left",column_names_side = "top",
                         row_names_gp = gpar(fontsize = 6),column_names_gp = gpar(fontsize = 6),
                         row_labels=gsub("_"," ",rownames(collect_stats)),
                         column_labels=gsub("_"," ",colnames(collect_stats)),
                         row_split=gaps_row,row_gap = unit(c(3,rep(1,length(unique(gaps_row))-1)), "mm"),
                         column_split=gaps_col,column_gap = unit(c(rep(1,length(unique(gaps_row))-1)), "mm"),
                         border = TRUE,
                         cell_fun = function(j, i, x, y, w, h, col) {
                           grid.text(signif_annot[i, j], x, y,gp = gpar(fontsize = 8))
                         },
                         row_title = c("Proportion",
                                       rep(" ",(length(unique(gaps_row))-2)/2),
                                       "Ratio between cell type proportion (out of total)",
                                       rep(" ",length(unique(gaps_row))-2-length(rep(" ",(length(unique(gaps_row))-2)/2)))),
                         row_title_gp = gpar(fontsize = 12),
                         width = unit(14, "cm"), height = unit(16, "cm"),
                         top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = "white"),
                                                                             labels = cat,
                                                                             labels_gp = gpar(col = "black", fontsize = 5))))
    }
  }

  row_names <- rownames(collect_stats)  # Save row names
  collect_stats[] <- lapply(collect_stats, function(x) {
    if (is.numeric(x)) x[is.nan(x)] <- NA
    return(x)
  })
  rownames(collect_stats) <- row_names

  return(list("heatmap"=heatmap,"statistics"=collect_stats))

}

#' @export
Compute_Proportions_Ratios <- function(single_cell_data,
                                       sample_colname,
                                       stratification,
                                       group_up,
                                       group_dn,
                                       stratification_for_ratio = NULL){

  meta <- single_cell_data

  ## Define categorizing factor and groups
  #########################################################################
  # Check that all parameters were filled correctly and are available
  if(1){
    if(is.null(stratification)){
      warning("Please provide the name of the column containing the stratification group that will be compared in the single-cell metadata object under the `stratification` parameter")
      return()
    }
    if(length(which(stratification  %in% colnames(meta) ==F))>0){
      warning("The selected stratification name parameter was not found in col.names of the single-cell object")
      return()
    }
    if(is.null(group_up)){
      warning("Please provide the stratification factor group down `group_up` that will be used for comparison with `group_dn` contained in the `stratification` column")
      return()
    }
    idx<-which(meta[,stratification]==group_up)
    if(length(idx)==0){
      warning(paste0("Stratification factor not found for `group_up`; available ones are: ",paste(unique(meta[,stratification]),collapse=", ")
      ))
      return()
    }
    if(is.null(group_dn)){
      warning("Please provide the stratification factor group down `group_dn` that will be used for comparison with `group_up` contained in the `stratification` column")
      return()
    }
    idx<-which(meta[,stratification]==group_dn)
    if(length(idx)==0){
      warning(paste0("Stratification factor not found for `group_dn`; available ones are: ",paste(unique(meta[,stratification]),collapse=", ")
      ))
      return()
    }
  }

  # Creating Sample annotation file
  Sample_annot <- meta %>%
    group_by(across(sample_colname)) %>%
    summarise(across(all_of(stratification), unique))

  Sample_annot <- as.data.frame(Sample_annot)
  Sample_annot$Stratif <- Sample_annot[,stratification]

  idx_up <- which(Sample_annot$Stratif %in% group_up)
  idx_dn <- which(Sample_annot$Stratif %in% group_dn)
  if(length(intersect(idx_up,idx_dn)) > 0){
    warning("Group.up and Group_dn are overlapping - Please provide non-overlapping categories")
    return()
  }

  Sample_annot$Stratif[idx_up] <- "up"
  Sample_annot$Stratif[idx_dn] <- "dn"
  Sample_annot <- Sample_annot[which(Sample_annot$Stratif %in% c("up","dn")),]
  Sample_annot$Stratif <- factor(Sample_annot$Stratif,levels = c("up","dn"))

  # Define patient levels
  sample_levels <- levels(factor(as.character(Sample_annot[,sample_colname])))
  meta[,sample_colname] <- factor(meta[,sample_colname], levels = sample_levels)
  #########################################################################

  # List and stratification of all cell types
  #########################################################################
  order_strat <- order(as.numeric(gsub("Cell_Type_Strat","",colnames(meta)[grep("Cell_Type_Strat",colnames(meta))])))
  cell_type_strats <- colnames(meta)[grep("Cell_Type_Strat",colnames(meta))][order_strat]

  column_cell_type<-cell_type_strats[length(cell_type_strats)]

  # Vector of granular cell types
  cell_subtypes <- levels(droplevels(meta[,cell_type_strats[length(cell_type_strats)]]))

  # Vector of bulked categories
  cell_types <- c()
  for(strata in rev(cell_type_strats[-length(cell_type_strats)])){
    tmp<-meta[order(match(meta[,column_cell_type],cell_subtypes)),strata]
    for(str in unique(tmp)){
      if(length(unique(meta[,column_cell_type][which(meta[,strata] == str)]))>1){
        cell_types<-c(cell_types,str)
      }
    }
  }
  if(length(which(duplicated(cell_types)))>0){
    cell_types <- cell_types[-which(duplicated(cell_types))]
  }

  # Rename the bulk cell types (the parent cell type that have more than one daughter cell type)
  meta_backup <- meta
  cell_types_bulk <- paste0(cell_types,"_bulk")
  for(strata in rev(cell_type_strats[-length(cell_type_strats)])){
    tmp <- meta[,strata]
    idx <- which(!is.na(match(tmp,cell_types)))
    tmp[idx] <- cell_types_bulk[match(tmp,cell_types)[idx]]
    meta[,strata] <- tmp
  }
  cell_types <- cell_types_bulk

  Collect_Proportion <- data.frame()

  collect_stats_ratio <- as.data.frame(matrix(nrow=length(c(cell_subtypes,cell_types)),ncol=length(c(cell_subtypes,cell_types)),NA))
  rownames(collect_stats_ratio)<-c(cell_subtypes,cell_types)
  colnames(collect_stats_ratio)<-c(cell_subtypes,cell_types)
  #########################################################################

  # Proportions out of total cells
  ########################################################################
  Proportions <- prop.table(table(meta[,which(colnames(meta) %in% c(sample_colname,column_cell_type))]),1)
  Proportions <- reshape2::melt(Proportions)
  Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
  colnames(Proportions)[2] <- "Cell_Type"
  Proportions$Type <- "out_of_:total"
  Collect_Proportion <- rbind(Collect_Proportion,Proportions)

  strats <- cell_type_strats[-length(cell_type_strats)]
  for(strat in strats){
    Proportions <- prop.table(table(meta[,which(colnames(meta) %in% c(sample_colname,strat))]),1)
    Proportions <- reshape2::melt(Proportions)
    Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
    colnames(Proportions)[2] <- "Cell_Type"
    Proportions$Type <- paste0("out_of_:total")
    Collect_Proportion <- rbind(Collect_Proportion,Proportions)
  }
  ########################################################################

  # Proportions out of stratification layers
  ########################################################################
  # Interrogate Fine cell subtypes out of all strats
  strats <- cell_type_strats[-length(cell_type_strats)]
  for(strat in strats){
    outofs <- unique(meta[,strat])
    for(outof in outofs){
      idx_outof <- which(meta[,strat] == outof)

      Proportions <- prop.table(table(meta[idx_outof,which(colnames(meta) %in% c(sample_colname,column_cell_type))]),1)
      Proportions <- reshape2::melt(Proportions)
      Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions)[2] <- "Cell_Type"
      Proportions$Type <- paste0("out_of_:",outof)
      Collect_Proportion <- rbind(Collect_Proportion,Proportions)
    }
  }

  # Interrogate Bulk populations in lower strats
  strats <- cell_type_strats[-length(cell_type_strats)]
  for(s in 1:length(strats)){
    strat <- strats[s]
    for(o in 1:(s-1)){
      outof_strat <- strats[o]

      outofs <- unique(meta[,outof_strat])
      for(outof in outofs){
        idx_outof <- which(meta[,outof_strat] == outof)

        Proportions <- prop.table(table(meta[idx_outof,which(colnames(meta) %in% c(sample_colname,strat))]),1)
        Proportions <- reshape2::melt(Proportions)
        Proportions$Stratif <- Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
        colnames(Proportions)[2] <- "Cell_Type"
        Proportions$Type <- paste0("out_of_:",outof)
        Collect_Proportion <- rbind(Collect_Proportion,Proportions)
      }
    }
  }
  ########################################################################

  # Ratio of cell types out of total cells
  ########################################################################

  # Subset Data For Ratio Computation
  if(!is.null(stratification_for_ratio)){
    if(!is.list(stratification_for_ratio)) {
      warning("The `stratification_for_ratio` parameter should be a list")
      return()
    }
    if(length(which(names(stratification_for_ratio)  %in% colnames(meta) ==F))>0){
      warning("The selected stratification_for_ratio parameters were not found in col.names of the single-cell object")
      return()
    }
    idx_to_keep <- 1:nrow(meta_backup)
    for(l in 1:length(stratification_for_ratio)){
      idx_to_keep <- intersect(idx_to_keep,which(meta_backup[,names(stratification_for_ratio)[l]] %in% stratification_for_ratio[[l]]))
    }
    meta_ratio <- meta[idx_to_keep,]
  } else {
    meta_ratio <- meta
  }

  # Compute Ratio for all Cell Type at all stratification and for all patients
  strats <- cell_type_strats
  for(s in 1:length(strats)){
    strat<-strats[s]
    if(s==1){
      Proportions<-prop.table(table(meta_ratio[,which(colnames(meta_ratio) %in% c(sample_colname,strat))]),1)
      Proportions<-reshape2::melt(Proportions)
      Proportions$Stratif<-Sample_annot$Stratif[match(Proportions[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions)[2]<-"Cell_Type"
    } else {
      Proportions_tmp<-prop.table(table(meta_ratio[,which(colnames(meta_ratio) %in% c(sample_colname,strat))]),1)
      Proportions_tmp<-reshape2::melt(Proportions_tmp)
      Proportions_tmp$Stratif<-Sample_annot$Stratif[match(Proportions_tmp[,sample_colname],Sample_annot[,sample_colname])]
      colnames(Proportions_tmp)[2]<-"Cell_Type"

      idx <- which(unique(Proportions_tmp$Cell_Type) %in% unique(Proportions$Cell_Type))
      if(length(idx)>0){
        Proportions_tmp<-Proportions_tmp[-which(Proportions_tmp$Cell_Type %in% unique(Proportions_tmp$Cell_Type)[idx]),]
      }
      Proportions<-rbind(Proportions,Proportions_tmp)
      rm(Proportions_tmp)
    }
  }

  # Compute and collect statistics for Ratio comparisons
  for(ct_row in rownames(collect_stats_ratio)){
    for(ct_col in colnames(collect_stats_ratio)){

      tmp1 <- Proportions[which(Proportions$Cell_Type==ct_col),]
      tmp2 <- Proportions[which(Proportions$Cell_Type==ct_row),]
      tmp2 <- tmp2[match(tmp1[,sample_colname],tmp2[,sample_colname]),]

      if(nrow(tmp1) > 0){
        tmp1$Ratio<-log2(tmp1$value/tmp2$value)

        tmp1$value<-tmp1$Ratio;tmp1$Ratio<-NULL
        tmp1$Type<-paste0("ratio:",ct_col,":",ct_row)
        Collect_Proportion<-rbind(Collect_Proportion,tmp1)
      }
    }
  }
  ########################################################################

  # Save objects
  Collect_Proportion$value[is.nan(Collect_Proportion$value)] <- NA
  Collect_Proportion[,1] <- as.character(Collect_Proportion[,1])
  Collect_Proportion$Cell_Type <- as.character(Collect_Proportion$Cell_Type)

  # Remove duplicates
  Collect_Proportion$concatenated <- paste0(Collect_Proportion[,1],
                                            Collect_Proportion$Cell_Type,
                                            Collect_Proportion$Type)
  if(length(which(duplicated(Collect_Proportion$concatenated)))>0){
    Collect_Proportion<-Collect_Proportion[-which(duplicated(Collect_Proportion$concatenated)),]
  }
  Collect_Proportion$concatenated<-NULL

  return(Collect_Proportion)
}

#' @export
CellType_Proportion_Boxplot <- function(cell_proportion_object,
                                        subset_data = NULL,
                                        stratification,
                                        groups = list(),
                                        type = c("Proportion","Ratio"),
                                        cell_type = NULL,
                                        out_of = NULL,
                                        cell_type_ratio = NULL,
                                        col_var){

  require(ComplexHeatmap)
  require(limma)
  require(reshape2)
  require(crayon)
  require(circlize)
  require(usethis)
  require(scales)
  require(dplyr)

  if(is.null(cell_proportion_object)){
    warning("Please provide a `cell_proportion_object` by running the `Compute_Proportions_Ratios` function")
    return()
  }
  Collect_Proportion <- cell_proportion_object

  if(is.na(match(type,c("Proportion","Ratio")))){
    warning("The selected type should be either `Proportion` or `Ratio`")
    return()
  }

  ## Subset Data
  #########################################################################
  if(!is.null(subset_data)){
    if (!is.list(subset_data)) {
      warning("The `subset_data` parameter should be a list")
      return()
    }
    if(length(which(names(subset_data)  %in% colnames(Collect_Proportion) ==F))>0){
      warning("The selected subset_data parameters were not found in col.names of the `cell_proportion_object` object")
      return()
    }
    idx_to_keep <- 1:nrow(Collect_Proportion)
    for(l in 1:length(subset_data)){
      idx_to_keep<-intersect(idx_to_keep,which(Collect_Proportion[,names(subset_data)[l]] %in% subset_data[[l]]))
    }
    Collect_Proportion<-Collect_Proportion[idx_to_keep,]
  }

  if(type=="Proportion"){

    idx<-grep("^out_of_:",Collect_Proportion$Type)
    Collect_Proportion<-Collect_Proportion[idx,]

    if(is.na(match(out_of,unique(gsub("^out_of_:","",Collect_Proportion$Type))))){
      warning(paste0("The selected value for the `out_of` parameter should be one of :",paste(unique(gsub("^out_of_:","",Collect_Proportion$Type)),collapse=", ")))
      return()
    }
    idx<-grep(paste0("^out_of_:",out_of),Collect_Proportion$Type)
    Collect_Proportion<-Collect_Proportion[idx,]

    if(is.na(match(cell_type,unique(Collect_Proportion$Cell_Type)))){
      warning(paste0("The selected value for the `cell_type` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type),collapse=", ")))
      return()
    }
    idx<-which(Collect_Proportion$Cell_Type %in% cell_type)
    Collect_Proportion<-Collect_Proportion[idx,]

    y_value<-paste0("Proportion out of ",out_of)
    title<-paste0(cell_type," out of ",out_of)
  }

  if(type=="Ratio"){

    idx<-grep("ratio",Collect_Proportion$Type)
    Collect_Proportion<-Collect_Proportion[idx,]
    Collect_Proportion$Type<-gsub("ratio:","",Collect_Proportion$Type)
    Collect_Proportion$Cell_Type1<-strsplit2(Collect_Proportion$Type,split="[:::]")[,1]
    Collect_Proportion$Cell_Type2<-strsplit2(Collect_Proportion$Type,split="[:::]")[,2]

    if(is.na(match(cell_type,unique(Collect_Proportion$Cell_Type1)))){
      warning(paste0("The selected value for the `cell_type` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type1),collapse=", ")))
      return()
    }

    if(is.na(match(cell_type_ratio,unique(Collect_Proportion$Cell_Type2)))){
      warning(paste0("The selected value for the `cell_type_ratio` parameter should be one of :",paste(unique(Collect_Proportion$Cell_Type2),collapse=", ")))
      return()
    }
    idx1<-which(Collect_Proportion$Cell_Type1 %in% cell_type)
    idx2<-which(Collect_Proportion$Cell_Type2 %in% cell_type_ratio)
    idx<-intersect(idx1,idx2)
    Collect_Proportion<-Collect_Proportion[idx,]

    y_value<-paste0("Ratio ",cell_type," / ",cell_type_ratio)
    title<-paste0("Ratio ",cell_type," / ",cell_type_ratio)
  }

  if(length(which(stratification  %in% colnames(Collect_Proportion) ==F))>0){
    warning("`stratification` factor not found in the col.names of the `cell_proportion_object` object")
    return()
  }

  if(length(which(unlist(groups)  %in% Collect_Proportion[,stratification]==F))>0){
    warning(paste0("The selected groups should be in : ",paste(unique(Collect_Proportion[,stratification]),collapse=", ")))
    return()
  }
  #########################################################################

  ## Collecting values for boxplots
  #########################################################################
  Collect_Proportion$pch<-21
  Collect_Proportion$col<-"lightgray"
  Collect_Proportion<-Collect_Proportion[complete.cases(Collect_Proportion),]
  Collect_Proportion$x<-NA

  list_Imm<-list()
  max<-c()
  min<-c()
  for(j in 1:length(groups)){
    idx<-which(Collect_Proportion[,stratification] %in% groups[[j]])
    tmp<-Collect_Proportion$value[idx]
    tmp <- tmp[!is.na(tmp)][!is.infinite(tmp[!is.na(tmp)])]
    assign(paste0("Imm",j),tmp)
    list_Imm[[j]]<-tmp[!is.na(tmp)]
    max<-c(max,tmp[!is.na(tmp)])
    min<-c(min,tmp[!is.na(tmp)])
    Collect_Proportion$x[idx]<-j
  }
  max <- max(max)
  min <- min(min)
  range <- (max-min)

  # Go out of the function if any of the vector has a zero length
  if (any(sapply(mget(ls(pattern = "^Imm\\d+$"), inherits = TRUE), length) == 0)) {
    warning("One condition had zero values without NAs or infinite")
    return()
  }
  #########################################################################

  ## Boxplot
  #########################################################################
  b<-boxplot(list_Imm,at=1:length(groups),col=col_var,
             names=names(groups),main=title,ylab=y_value, xaxt = "n",whisklty=1,lwd=1.5,
             ylim=c(min,max+(0.3*range)), outline = FALSE,width=rep(0.3,length(groups)))
  axis(side = 1, at = seq_along(b$names), labels = b$names, tick = FALSE)

  for(i in 1:nrow(Collect_Proportion)){
    stripchart(Collect_Proportion$value[i],
               vertical = TRUE, method = "jitter",jitter=0.2,at=Collect_Proportion$x[i],
               pch = Collect_Proportion$pch[i],
               bg=alpha(Collect_Proportion$col[i],1),col = "black",
               add = TRUE,cex=1.4)
  }

  if(length(groups)==2){
    used_test<-"wilcox.test"
    res<-do.call(used_test,list(Imm1,Imm2))
    text(1.2,max+(0.05*range),labels = paste0("p-val: ",signif(res$p.value,digits = 2)),cex=0.8,pos=4)
  }


  if(length(groups)==3){
    adjusted<-F
    pval<-as.data.frame(pairwise.wilcox.test(Collect_Proportion$value, Collect_Proportion[,stratification],p.adjust.method = ifelse(adjusted,"fdr","none"))$p.value)
    pval<-signif(pval,digit=2)

    text(2,max+(0.3*range),colnames(pval)[1],cex=0.8)
    text(3,max+(0.3*range),colnames(pval)[2],cex=0.8)
    text(1,max+(0.23*range),rownames(pval)[1],cex=0.8)
    text(1,max+(0.16*range),rownames(pval)[2],cex=0.8)

    text(2,max+(0.23*range),pval[1,1],cex=0.7)
    text(2,max+(0.16*range),pval[2,1],cex=0.7)
    text(3,max+(0.16*range),pval[2,2],cex=0.7)

  }

  if(length(groups)==4){
    adjusted<-F
    pval<-as.data.frame(pairwise.wilcox.test(Collect_Proportion$value, Collect_Proportion[,stratification],p.adjust.method = ifelse(adjusted,"fdr","none"))$p.value)
    pval<-signif(pval,digit=2)

    text(2,max+(0.3*range),colnames(pval)[1],cex=0.8)
    text(3,max+(0.3*range),colnames(pval)[2],cex=0.8)
    text(4,max+(0.3*range),colnames(pval)[3],cex=0.8)
    text(1,max+(0.23*range),rownames(pval)[1],cex=0.8)
    text(1,max+(0.16*range),rownames(pval)[2],cex=0.8)
    text(1,max+(0.09*range),rownames(pval)[3],cex=0.8)

    text(2,max+(0.23*range),pval[1,1],cex=0.7)
    text(2,max+(0.16*range),pval[2,1],cex=0.7)
    text(2,max+(0.09*range),pval[3,1],cex=0.7)
    text(3,max+(0.16*range),pval[2,2],cex=0.7)
    text(3,max+(0.09*range),pval[3,2],cex=0.7)
    text(4,max+(0.09*range),pval[3,3],cex=0.7)
  }
  #########################################################################

}
