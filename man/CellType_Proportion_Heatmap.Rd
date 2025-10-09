\name{CellType_Proportion_Heatmap}
\alias{CellType_Proportion_Heatmap}
\title{Creation of a heatmap displaying Cell Type proportions and ratios between cell type proportions in a single-cell annotated object}
\description{The function will interrogate the indicated cohort(s) for any potential cell type enrichment in the two groups provided. The cell type proportion are shown according to several levels of resolution ranging from the overall stoechiometry to their proportion within their cell lineage. The reported heatmap also shows ratio between cell type proportions. The function supports cancer type selection within the cohort and exclusion of specific patients.}
\usage{CellType_Proportion_Heatmap <- function(single_cell_data,
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
                                        stratification_for_ratio = NULL)}

\arguments{
  \item{single_cell_data}{a data.frame (typically the meta.data from a Seurat object) where rows are individual cells and columns contain informations about cell types, their stratifications, the sample origin and the sample information needed to perform the comparison and for potential subsetting of the data. The cell type information has to be encoded in columns with the `Cell_Type_Strat` pattern followed by a numeric (e.g. Cell_Type_Strat1, Cell_Type_Strat2, Cell_Type_Strat3 and Cell_Type_Strat4) where the lower stratification (Cell_Type_Strat1) encodes the less granular cell types (e.g. `immune` and `non_immune`) and the highest stratification (Cell_Type_Strat4) encodes for the highest granularity (e.g. `exhausted_T_cell`, `CD4_TRegs` or `Pericytes`). The highest granularity stratification has to be a factor with the levels encoding the order of the cell types as they will appear in the heatmap.}
  \item{sample_colname}{a character indicating which column name of the single_cell_data object corresponds to the column with the sample identifiers.}
  \item{subset_data}{if provided, should be a list indicating which data to subset/use (such as list(`Technology`=c("5prime"),`Patient`=c(`Patient1`,`Patient2`))) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed.}
  \item{stratification}{a character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of the single_cell_data object.}
  \item{group_up}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_dn` group. The cell types enriched in this group will appear in red (unless other color is specific in `group_up_col`) in the heatmap}
  \item{group_dn}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_up` group. The cell types enriched in this group will appear in blue (unless other color is specific in `group_dn_col`) in the heatmap}
  \item{group_up_col}{if provided, a color that will be used to highlight cell types enriched in the `group_up`.}
  \item{group_dn_col}{if provided, a color that will be used to highlight cell types enriched in the `group_dn`.}
  \item{statistics}{which statistical test to use. Should be either "wilcoxon" or "ttest". Default set to wilcoxon.}
  \item{include_ratio}{a boolean indicating if ratio between cell type proportions should be added to the heatmap. If set to FALSE, only cell type proportions will be shown. Default set to TRUE.}
  \item{stratification_names}{a vector that can used to supply the names of the stratification so that they are more meaningful that `Cell_Type_Strat1`, `Cell_Type_Strat2` etc... . It can be something like c(`Immune_NotImmune`,`Tcell_Bcell_Malignant_Stromal`,`CD4_CD8_Malignant_CAFs_Endo`).}
  \item{pval_range}{a numerical setting the maxima of the -log10 p-values in the color code. That will force the range of -log10 p-values (encoded by the color code) to be contain within this range. If set to 3, the color code will range between -3 and 3.}
  \item{stratification_for_ratio}{if provided, should be a list indicating which data to subset/use (such as list(`Cell_Type_Strat1`=c("Immune")) ) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed for comouting the ratio of between cell proportions.}
}
\details{}
\value{a list containing a heatmap (`heatmap`) and the statistics p-value table used to make the heatmap (`statistics`). The output can be plotted directly by calling the function or alternatively saved into a variable and plot using the draw() function. The values in the heatmap are showing directional log10 of the p-values. The red color-code indicates enrichment in the `group_up` group while blue color indicates enrichment in the `group_dn` group. Significance is defined, without adjustment for multiple comparisons, by `*` = p-value < 0.05, `**` = p-value < 0.01 and `***` = p-value < 0.001.}
\references{https://www.science.org/doi/10.1126/sciimmunol.adg7995}
\author{David Barras}
\note{}
\seealso{}
\examples{
data("example_data", package = "CellProportionR")
Proportion_Heatmap <- CellType_Proportion_Heatmap(single_cell_data = single_cell_data,
                            sample_colname = "Patient",
                            subset_data = list("Technology"=c("5prime"),
                                               "Patient"=c("Patient11","Patient27","Patient28","Patient32",
                                                           "Patient31","Patient17","Patient25","Patient30")),
                            stratification = "Category",
                            group_up = c("group1"),
                            group_dn=c("group2"),
                            stratification_names = c("Immune_vs_Not_vs_Doublet",
                                                     "Mal_T_B_Myeloid_Stromal_NK",
                                                     "Mal_CD8_CD4_DN_Tgd_B_Macro_Mono_DC_Endo_CAFs_Pericyte",
                                                     "Fine_Cell_Types"))
Proportion_Heatmap$heatmap
pvalues <- Proportion_Heatmap$statistics}

