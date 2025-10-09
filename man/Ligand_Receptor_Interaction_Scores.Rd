\name{Ligand_Receptor_Interaction_Scores}
\alias{Ligand_Receptor_Interaction_Scores}
\title{Calculate Ligand-Receptor Interaction Scores Between Cell Types}
\description{
  Computes interaction scores between all pairs of cell types based on ligand-receptor gene expression and cell type proportions. The interaction score for each ligand-receptor pair between two cell types is calculated as the product of: (1) proportion of cell type 1, (2) average expression of gene 1 in cell type 1, (3) proportion of cell type 2, and (4) average expression of gene 2 in cell type 2.
}
\usage{
Ligand_Receptor_Interaction_Scores(single_cell_metadata,
                                   single_cell_gex,
                                   sample_colname,
                                   celltype_colname,
                                   subset_data = NULL,
                                   lr_database = NULL,
                                   celltype_to_exclude = NULL)
}
\arguments{
  \item{single_cell_metadata}{
    A data.frame (typically the meta.data from a Seurat object) where rows are individual cells and columns contain information about cell types, sample origin, and other metadata. Row names must match column names of \code{single_cell_gex}.
  }
  \item{single_cell_gex}{
    A matrix or sparse matrix (dgCMatrix) of single-cell gene expression with genes in rows and cells in columns. Column names must match row names of \code{single_cell_metadata}. Gene names should include those present in the ligand-receptor database.
  }
  \item{sample_colname}{
    Character string specifying the column name in \code{single_cell_metadata} that contains sample or patient identifiers.
  }
  \item{celltype_colname}{
    Character string specifying the column name in \code{single_cell_metadata} that contains cell type annotations.
  }
  \item{subset_data}{
    Optional list for subsetting data. Each element should be named according to a column in \code{single_cell_metadata} and contain a vector of values to filter by. For example: \code{list(Treatment = c("Control", "Drug"), Tissue = "Tumor")}.
  }
  \item{lr_database}{
    Optional data.frame with two columns named "receptor" and "ligand" containing ligand-receptor pairs. If NULL, uses the default database from the scAvengeR package.
  }
  \item{celltype_to_exclude}{
    Optional character vector of cell type names to exclude from the analysis.
  }
}
\details{
  The function performs the following steps:
  \enumerate{
    \item Validates input data consistency and checks for matching cells between metadata and expression matrix
    \item Filters genes to those present in both the expression matrix and ligand-receptor database
    \item Computes cell type proportions per sample
    \item Creates pseudobulk expression profiles by averaging gene expression per cell type per sample
    \item Calculates interaction scores for all possible cell type pairs and ligand-receptor combinations
    \item Sets scores to 0 when cell types are absent in a sample
  }

  The interaction score formula is:
  \deqn{Score = Prop_{CT1} \times Expr_{Gene1,CT1} \times Prop_{CT2} \times Expr_{Gene2,CT2}}

  where Prop is the cell type proportion and Expr is the average gene expression.
}
\value{
  A list containing three elements:
  \item{interaction_score}{
    Data.frame with rows representing unique cell type-gene pairs and columns representing samples. Each value is the computed interaction score.
  }
  \item{possible_combination}{
    Data.frame detailing all possible combinations with columns: celltype1, gene1, celltype2, gene2, and unique_id.
  }
  \item{proportion_cell_type}{
    Matrix showing the proportion of each cell type in each sample.
  }
}
\references{
  \url{https://www.science.org/doi/10.1126/sciimmunol.adg7995}
}
\author{
  David Barras
}
\note{
  Large datasets may require significant memory. The function will issue warnings if:
  \itemize{
    \item Cell names don't match between metadata and expression matrix
    \item Required columns are missing from metadata
    \item Few ligand-receptor genes are found in the expression data
  }
}
\seealso{
  \code{\link{Ligand_Receptor_Differential_Analysis}} for differential analysis of interaction scores
}
\examples{
\dontrun{
# Load example single-cell data
library(Seurat)
library(scAvengeR)

# Prepare input from Seurat object
metadata <- seurat_obj@meta.data
expression <- GetAssayData(seurat_obj, layer = "data")

# Calculate interaction scores
interaction_results <- Ligand_Receptor_Interaction_Scores(
  single_cell_metadata = metadata,
  single_cell_gex = expression,
  sample_colname = "patient_id",
  celltype_colname = "cell_type",
  subset_data = list(tissue = "tumor"),
  celltype_to_exclude = c("Doublets", "Unknown")
)

# Access results
interaction_scores <- interaction_results$interaction_score
cell_proportions <- interaction_results$proportion_cell_type
}
}
\keyword{cell communication}
\keyword{ligand-receptor}
\keyword{single-cell}
