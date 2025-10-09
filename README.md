# scAvengeR

## Overview

**scAvengeR** is a comprehensive R package for single-cell RNA-seq analysis, providing a collection of utilities for investigating cell type proportions and cell-cell communication through ligand-receptor interactions. This package combines and extends functionality for both statistical analysis and visualization of single-cell datasets.

The package includes two main analysis modules:
- **Cell Proportion Analysis**: Statistical comparison and visualization of cell type proportions between conditions
- **Ligand-Receptor Interaction Analysis**: Quantification and differential analysis of cell-cell communication

## Installation

To install **scAvengeR**, use the following command in R:

```r
# Install from GitHub
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("dbarras2/scAvengeR")
```

## Module 1: Cell Proportion Analysis

These functions enable comprehensive analysis of cell type proportions across different conditions, with hierarchical stratification support and multiple visualization options.

### 1.1 CellType_Proportion_Heatmap

**Description:** Generates a heatmap displaying statistical comparisons of cell type proportions and ratios between two groups in single-cell data.

#### Example
<img width="730" alt="Heatmap_example" src="https://github.com/user-attachments/assets/b4f3a359-10c9-4cf4-9ca3-351fe030547a" />

#### Usage
```r
data("example_data", package = "scAvengeR")
Proportion_Heatmap <- CellType_Proportion_Heatmap(
  single_cell_data = single_cell_data,
  sample_colname = "Patient",
  subset_data = list(
    Technology = c("5prime"),
    Patient = c("Patient11", "Patient27", "Patient28", "Patient32",
                "Patient31", "Patient17", "Patient25", "Patient30")
  ),
  stratification = "Category",
  group_up = c("group1"),
  group_dn = c("group2"),
  stratification_names = c(
    "Immune_vs_Not_vs_Doublet",
    "Mal_T_B_Myeloid_Stromal_NK",
    "Mal_CD8_CD4_DN_Tgd_B_Macro_Mono_DC_Endo_CAFs_Pericyte",
    "Fine_Cell_Types"
  )
)
Proportion_Heatmap$heatmap
pvalues <- Proportion_Heatmap$statistics
```

### 1.2 Compute_Proportions_Ratios

**Description:** Computes proportions and ratios of different cell types within a single-cell dataset for downstream analysis.

#### Usage
```r
Proportion_object <- Compute_Proportions_Ratios(
  single_cell_data = single_cell_data,
  sample_colname = "Patient",
  stratification = "Category",
  group_up = "group1",
  group_dn = "group2"
)
```

### 1.3 CellType_Proportion_Boxplot

**Description:** Creates boxplots comparing cell type proportions or ratios across groups with statistical testing.

#### Usage
```r
# Proportion boxplot
CellType_Proportion_Boxplot(
  cell_proportion_object = Proportion_object,
  subset_data = list(excluded = "no"),
  stratification = "Stratif",
  groups = list(up = "up", dn = "dn"),
  type = "Proportion",
  cell_type = "CD8_Tex",
  out_of = "CD8_T_cells",
  col_var = c("red", "blue")
)

# Ratio boxplot
CellType_Proportion_Boxplot(
  cell_proportion_object = Proportion_object,
  subset_data = list(excluded = "no"),
  stratification = "Stratif",
  groups = list(up = "up", dn = "dn"),
  type = "Ratio",
  cell_type = "CD8_EM-like",
  cell_type_ratio = "Neutrophils",
  col_var = c("red", "blue")
)
```

## Module 2: Ligand-Receptor Interaction Analysis

These functions enable comprehensive analysis of cell-cell communication through ligand-receptor interactions, with support for differential analysis and pathway-level investigation.

### 2.1 Ligand_Receptor_Interaction_Scores

**Description:** Calculates interaction scores between all pairs of cell types based on ligand-receptor gene expression and cell type proportions. The score represents the potential for cell-cell communication based on the product of cell proportions and gene expression levels.

#### Usage
```r
# Calculate interaction scores from Seurat object
library(Seurat)
metadata <- seurat_obj@meta.data
expression <- GetAssayData(seurat_obj, layer = "data")

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
```

### 2.2 Ligand_Receptor_Differential_Analysis

**Description:** Performs statistical differential analysis of ligand-receptor interactions between conditions and creates comprehensive visualizations including circos plots, heatmaps, and bar plots.

#### Example Visualizations
[INSERT YOUR LIGAND-RECEPTOR FIGURE HERE]

#### Usage
```r
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

# Access visualizations
print(diff_results$Barplot_Total_Interactions)
print(diff_results$circos_group_up)
print(diff_results$circos_group_dn)
print(diff_results$Heatmap_Interacion_All_LR)

# Get significant interactions
sig_interactions <- diff_results$Differential_Analysis[
  diff_results$Differential_Analysis$p.value < 0.05, 
]
```

## Key Features

### Cell Proportion Analysis
- Hierarchical cell type stratification support (from broad to fine-grained annotations)
- Multiple statistical tests (Wilcoxon, t-test)
- Proportion and ratio calculations
- Customizable heatmaps with significance indicators
- Group-wise boxplot comparisons

### Ligand-Receptor Interaction Analysis
- Interaction score calculation based on cell proportions and gene expression
- Three statistical test options (t-test, Wilcoxon, linear model)
- Pathway-level interaction grouping
- Multiple visualization types:
  - Circos plots for interaction networks
  - Heatmaps for interaction patterns
  - Bar plots for quantitative comparisons
- Support for custom ligand-receptor databases

## Dependencies

### Core Dependencies
- `limma`
- `reshape2`
- `dplyr`
- `Matrix`

### Visualization Dependencies
- `ComplexHeatmap`
- `ggplot2`
- `patchwork`
- `circlize`
- `RColorBrewer`
- `colorspace`
- `scales`

### Additional Dependencies for Ligand-Receptor Analysis
- `pbapply`
- `usethis`

## Input Data Format

### For Cell Proportion Analysis
- Single-cell metadata with hierarchical cell type annotations (Cell_Type_Strat1, Cell_Type_Strat2, etc.)
- Sample identifiers
- Grouping variables for comparison

### For Ligand-Receptor Analysis
- Single-cell gene expression matrix (genes × cells)
- Cell type annotations
- Sample identifiers
- Optional: Custom ligand-receptor database

## Citations

If you use **scAvengeR** in your research, please cite:

```
Barras, D. et al. (2024). scAvengeR: A comprehensive toolkit for single-cell 
analysis of cell proportions and cell-cell communication. 
https://github.com/dbarras2/scAvengeR

Reference publication:
https://www.science.org/doi/10.1126/sciimmunol.adg7995
```

## Author

Developed by **David Barras**. 

## Contributing

Contributions are welcome! Please feel free to submit issues, feature requests, or pull requests on the [GitHub repository](https://github.com/dbarras2/scAvengeR).

## License

This package is licensed under the MIT License.

## Support

For questions, bug reports, or feature requests, please open an issue on the [GitHub repository](https://github.com/dbarras2/scAvengeR) or contact the author directly.