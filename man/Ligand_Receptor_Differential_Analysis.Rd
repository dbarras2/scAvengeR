\name{Ligand_Receptor_Differential_Analysis}
\alias{Ligand_Receptor_Differential_Analysis}
\title{Differential Analysis of Ligand-Receptor Interactions Between Groups}
\description{
  Performs statistical differential analysis of ligand-receptor interaction scores between two groups (e.g., treatment vs control, responder vs non-responder). Creates comprehensive visualizations including circos plots, heatmaps, and bar plots to illustrate significant differences in cell-cell communication patterns.
}
\usage{
Ligand_Receptor_Differential_Analysis(interaction_score_object,
                                      group_up,
                                      label_up = NULL,
                                      group_dn,
                                      label_dn = NULL,
                                      statistics = c("ttest", "wilcoxon", "lm"),
                                      lr_pathways_to_plot = NULL,
                                      color_pathways = NULL,
                                      cell_type_groups)
}
\arguments{
  \item{interaction_score_object}{
    Output from \code{Ligand_Receptor_Interaction_Scores} function containing interaction_score and possible_combination data.frames.
  }
  \item{group_up}{
    Character vector of sample names for the first group (e.g., treatment group). These samples must exist in the columns of the interaction_score matrix. Interactions enriched in this group will be labeled as "up".
  }
  \item{label_up}{
    Optional character string to label the first group in plots. Default is "up".
  }
  \item{group_dn}{
    Character vector of sample names for the second group (e.g., control group). These samples must exist in the columns of the interaction_score matrix. Interactions enriched in this group will be labeled as "down".
  }
  \item{label_dn}{
    Optional character string to label the second group in plots. Default is "down".
  }
  \item{statistics}{
    Statistical test to use for differential analysis. Options are:
    \describe{
      \item{"ttest"}{Two-sample t-test (parametric)}
      \item{"wilcoxon"}{Wilcoxon rank-sum test (non-parametric)}
      \item{"lm"}{Linear model using limma package}
    }
  }
  \item{lr_pathways_to_plot}{
    Optional named list of data.frames, each containing "receptor" and "ligand" columns defining pathway-specific ligand-receptor pairs. If NULL, uses default pathways from scAvengeR package.
  }
  \item{color_pathways}{
    Optional character vector of colors for each pathway in \code{lr_pathways_to_plot}. Must be same length as number of pathways. If NULL, default colors are used.
  }
  \item{cell_type_groups}{
    Named list defining cell type groupings. Each element should be a character vector of cell types belonging to that group. For example: \code{list(Immune = c("T_cells", "B_cells", "Macrophages"), Stromal = c("Fibroblasts", "Endothelial"))}.
  }
}
\details{
  The function performs comprehensive differential analysis:

  \enumerate{
    \item \strong{Statistical Testing}: Compares interaction scores between groups using the specified test
    \item \strong{Multiple Visualizations}:
      \itemize{
        \item Bar plots showing total and significant interactions per group
        \item Circos plots displaying significant interactions by pathway
        \item Heatmaps of interaction counts between cell types
        \item Bar plots of interactions grouped by pathways
      }
    \item \strong{Chi-square Test}: Tests for overall enrichment of interactions between groups
    \item \strong{Pathway Analysis}: Groups interactions by biological pathways for interpretation
  }

  Significance threshold is set at p-value < 0.05. Directionality indicates which group has higher interaction scores.
}
\value{
  A list containing nine elements:
  \item{Differential_Analysis}{
    Data.frame with all interactions and their statistical results, including directionality and p-values.
  }
  \item{Barplot_Total_Interactions}{
    ggplot object showing total and significant interactions per group with chi-square test result.
  }
  \item{circos_group_up}{
    Circos plot showing significant interactions enriched in group_up.
  }
  \item{circos_group_dn}{
    Circos plot showing significant interactions enriched in group_dn.
  }
  \item{Significant_Interaction_Number}{
    Data.frame summarizing the number of significant interactions by cell type pairs and pathways.
  }
  \item{Barplot_Interaction_Pathway_Fine}{
    ggplot showing interactions for fine-grained cell type pairs grouped by pathway.
  }
  \item{Barplot_Interaction_Pathway_Main}{
    ggplot showing interactions for main cell type groups grouped by pathway.
  }
  \item{Heatmap_Interacion_Pathway}{
    Heatmap showing pathway-specific interaction counts between cell types.
  }
  \item{Heatmap_Interacion_All_LR}{
    Heatmap showing all ligand-receptor interaction counts between cell types.
  }
}
\references{
  \url{https://www.science.org/doi/10.1126/sciimmunol.adg7995}
}
\author{
  David Barras
}
\note{
  \itemize{
    \item Requires at least 2 samples per group for statistical testing
    \item Circos plots require the circlize package
    \item Heatmaps require the ComplexHeatmap package
    \item Large numbers of cell types may make circos plots crowded
    \item P-values are not adjusted for multiple testing by default
  }
}
\seealso{
  \code{\link{Ligand_Receptor_Interaction_Scores}} for calculating interaction scores
}
\examples{
\dontrun{
# First calculate interaction scores
interaction_results <- Ligand_Receptor_Interaction_Scores(
  single_cell_metadata = metadata,
  single_cell_gex = expression,
  sample_colname = "patient_id",
  celltype_colname = "cell_type"
)

# Define cell type groups
cell_groups <- list(
  Immune = c("CD4_T", "CD8_T", "B_cells", "NK", "Macrophages"),
  Stromal = c("Fibroblasts", "Endothelial"),
  Epithelial = c("Tumor", "Normal_epithelial")
)

# Perform differential analysis
diff_results <- Ligand_Receptor_Differential_Analysis(
  interaction_score_object = interaction_results,
  group_up = c("Patient1", "Patient2", "Patient3"),
  label_up = "Responders",
  group_dn = c("Patient4", "Patient5", "Patient6"),
  label_dn = "Non-responders",
  statistics = "wilcoxon",
  cell_type_groups = cell_groups
)

# Access specific results
significant_interactions <- diff_results$Differential_Analysis[
  diff_results$Differential_Analysis$p.value < 0.05, ]

# Plot individual visualizations
print(diff_results$Barplot_Total_Interactions)
print(diff_results$circos_group_up)
print(diff_results$Heatmap_Interacion_All_LR)
}
}
\keyword{differential analysis}
\keyword{cell communication}
\keyword{visualization}
