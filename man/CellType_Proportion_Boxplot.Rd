\name{CellType_Proportion_Boxplot}
\alias{CellType_Proportion_Boxplot}
\title{Creation of a boxplot displaying Cell Type proportions or ratio in a single-cell annotated object}
\description{The function will plot the cell type frequencies or ratio by splitting for the chosen stratification factor. The cell type proportion can be shown according to several levels of resolution ranging from the overall stoechiometry (total cells) to their proportion within their cell lineage. The ratio between different cell type proportion (out of total cells) can also be reported. The function supports subsetting of the data and exclusion of specific samples.}
\usage{CellType_Proportion_Boxplot(cell_proportion_object,
                                   subset_data = NULL,
                                   stratification,
                                   groups = list(),
                                   type = c("Proportion","Ratio"),
                                   cell_type = NULL,
                                   out_of = NULL,
                                   cell_type_ratio = NULL,
                                   col_var)}

\arguments{
  \item{cell_proportion_object}{a data.frame as output by the `Compute_Proportions_Ratios` function.}
  \item{subset_data}{if provided, should be a list indicating which data to subset/use (such as list(`Technology`=c("5prime"),`Patient`=c(`Patient1`,`Patient2`))) where each element is named with the same name than the column in the `cell_proportion_object` object. Important, don't forget to annotate the `cell_proportion_object` object if the subset information is not already present.}
  \item{stratification}{a character mentioning which stratification factor to use for splitting the data into groups. This character has to be one of the column name of `cell_proportion_object`.}
  \item{groups}{a list indicating all the groups that will be used to make the boxplots. Each group will represented by a box. Each element of the list can be a simple character or a vector indicating groups to form (e.g "Responders"=c("CR","PR")). The name of each element will be used as the name of the category for plotting.}
  \item{type}{a character indicating whether to plot `Proportion` or `Ratio`}
  \item{cell_type}{a character indicating which cell type to report. If type="Proportion", then the cell_type indicates the proportion of this cell type will be computed relative to the indicated cell population (see `out_of` below). If type="Ratio", then the cell_type will be the cell type that will be divided by the second cell type indicated by the `cell_type_ratio` parameter (see below). Of note, bulk cell types (parental cell types haveing more than one daughter cell type) were being added the "_bulk" extension to their name to avoid any misinterpretation with more granular/daughter cell type bearing the same name.}
  \item{out_of}{If type="Proportion", `out_of` should indicated to which cell population the frequency of cell_type should be computed. If the proportion of CD8_Pex out of T cells has to be computed, then cell_type="CD8_Pex" and out_of="T_cells". Any values that corresponds to any levels of stratification is allowed as long as it was computable. Of note, bulk cell types (parental cell types haveing more than one daughter cell type) were being added the "_bulk" extension to their name to avoid any misinterpretation with more granular/daughter cell type bearing the same name.}
  \item{cell_type_ratio}{If type="Ratio", `cell_type_ratio` indicates the cell type to normalize to. So the ratio would be the ratio of cell_type divided by cell_type_ratio. Importantly, the log2 ratio is reported, so positive values indicates that the proportion of cell_type is higher than the proportion of cell_type_ratio. If the ratio between CD8_Pex and CD4_CXCL13 has to be computed, then cell_type="CD8_Pex" and cell_type_ratio="CD4_CXCL13"}
  \item{col_var}{a vector indicating which colors should be used for the boxes. Should be the same length as the length of the `groups` parameter}
}
\details{}
\value{Running the function will plot the boxplot in the graphic device. p-values in the plot are computed using non-adjusted Wilcoxon rank-sum test}
\references{https://www.science.org/doi/10.1126/sciimmunol.adg7995}
\author{David Barras}
\note{}
\seealso{}
\examples{
data("example_data", package = "CellProportionR")
Proportion_object <- Compute_Proportions_Ratios(single_cell_data = single_cell_data,
                                                sample_colname = "Patient",
                                                stratification = "Category",
                                                group_up = "group1",
                                                group_dn = "group2")

Proportion_object$excluded<-"no"
Proportion_object$excluded[which(Proportion_object$Patient=="Patient1")]<-"yes"

CellType_Proportion_Boxplot(cell_proportion_object = Proportion_object,
                            subset_data = list("excluded"="no"),
                            stratification="Stratif",
                            groups = list("up"="up","dn"="dn"),
                            type = "Ratio",
                            cell_type = "CD8_EM-like",
                            out_of = NULL,
                            cell_type_ratio = "Neutrophils",
                            col_var = c("red","blue"))

CellType_Proportion_Boxplot(cell_proportion_object = Proportion_object,
                            subset_data = list("excluded"="no"),
                            stratification="Stratif",
                            groups = list("up"="up","dn"="dn"),
                            type = "Proportion",
                            cell_type = "CD8_Tex",
                            out_of = "CD8_T_cells",
                            cell_type_ratio = NULL,
                            col_var = c("red","blue"))}

