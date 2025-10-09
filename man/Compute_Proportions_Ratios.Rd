\name{Compute_Proportions_Ratios}
\alias{Compute_Proportions_Ratios}
\title{Computation of Cell Type proportions at all stratifications and ratios between cell type proportions}
\description{The function will compute the cell type proportions and their ratios at various stratification level (if provided) and per sample. This object can then be used for plotting using the CellType_Proportion_Boxplot function.}
\usage{Compute_Proportions_Ratios <- function(single_cell_data,
                                       sample_colname,
                                       stratification,
                                       group_up,
                                       group_dn,
                                       stratification_for_ratio = NULL)}

\arguments{
  \item{single_cell_data}{a data.frame (typically the meta.data from a Seurat object) where rows are individual cells and columns contain informations about cell types, their stratifications, the sample origin and the sample information needed to perform the comparison and for potential subsetting of the data. The cell type information has to be encoded in columns with the `Cell_Type_Strat` pattern followed by a numeric (e.g. Cell_Type_Strat1, Cell_Type_Strat2, Cell_Type_Strat3 and Cell_Type_Strat4) where the lower stratification (Cell_Type_Strat1) encodes the less granular cell types (e.g. `immune` and `non_immune`) and the highest stratification (Cell_Type_Strat4) encodes for the highest granularity (e.g. `exhausted_T_cell`, `CD4_TRegs` or `Pericytes`).}
  \item{sample_colname}{a character indicating which column name of the single_cell_data object corresponds to the column with the sample identifiers.}
  \item{stratification}{a character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of the single_cell_data object.}
  \item{group_up}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_dn` group. The cell types enriched in this group will appear in red (unless other color is specific in `group_up_col`) in the heatmap}
  \item{group_dn}{a vector or character value indicating which group to use in the selected stratification. This group will be compared against the `group_up` group. The cell types enriched in this group will appear in blue (unless other color is specific in `group_dn_col`) in the heatmap}
  \item{stratification_for_ratio}{if provided, should be a list indicating which data to subset/use (such as list(`Cell_Type_Strat1`=c("Immune")) ) where each element is named with the same name than the column in `single_cell_data` and each elements is a vector containing the identifiers contained in these columns and that allow to select the cells to be analyzed for comouting the ratio of between cell proportions.}
}
\details{}
\value{a data.frame containing the values for all cell proportion and ratios per sample.}
\references{https://www.science.org/doi/10.1126/sciimmunol.adg7995}
\author{David Barras}
\note{}
\seealso{}
\examples{
data("example_data", package = "CellProportionR")
Compute_Proportions_Ratios(single_cell_data = single_cell_data,
                           sample_colname = "Patient",
                           stratification = "Category",
                           group_up = "group1",
                           group_dn = "group2")
}

