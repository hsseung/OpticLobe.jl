# Parts List Paper Analysis Scripts

This directory contains Julia scripts for generating plots and analyses for the Nature papers on the Drosophila optic lobe parts list.

## Setup

This directory has its own Julia environment with all required dependencies. To set up:

```bash
cd PartsListPaper
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Running Scripts

Run any script using:

```bash
julia --project=. "script_name.jl"
```

For example:
```bash
julia --project=. "Fig 1 cell numbers.jl"
```

## Scripts Overview

### Main Figures
- **`Fig 1 cell numbers.jl`** - Generates cell count visualizations for main Figure 1
- **`Fig 2a featurevectors.jl`** - Creates feature vector plots for main Figure 2a

### Extended Data Figures
- **`Fig S1 histogram type cardinalities.jl`** - Histogram of cell type cardinalities (Extended Data Figure 1)
- **`Fig S3ab TypeRadii.jl`** - Type radius analysis plots (Extended Data Figure 3a-b)  
- **`Fig S4 Type-to-type connectivity as a matrix.jl`** - Connectivity matrix visualization (Extended Data Figure 4)

### Supplementary Data
- **`Data S3 Discriminators2D.jl`** - 2D discriminator analysis for Supplementary Data S3
- **`Data S4 input output heatmaps.jl`** - Input/output heatmap generation for Supplementary Data S4

### Analysis Tools
- **`Dendrogram.jl`** - Hierarchical clustering and dendrogram generation
- **`perplexity.jl`** - Perplexity calculations for dimensionality analysis

## Dependencies

The environment includes all necessary packages:
- **OpticLobe.jl** - Main analysis package (local development version)
- **Plotting**: Plots.jl, StatsPlots.jl, PlotlyKaleido.jl, Measures.jl
- **Data Processing**: CSV.jl, DataFrames.jl, NamedArrays.jl, SparseArrays.jl
- **Analysis**: Clustering.jl, Distances.jl, StatsBase.jl, FreqTables.jl
- **Specialized**: NewickTree.jl, Phylo.jl, MissingsAsFalse.jl

## Interactive Development

For interactive analysis:

```bash
julia --project=.
```

Then:
```julia
using OpticLobe
# All plotting and analysis packages are available
```

## Data Dependencies

Scripts automatically download required FlyWire v783 data via OpticLobe.jl's DataDeps system on first run. Data includes:
- Cell statistics and morphology
- Synaptic connectivity (Princeton synapses by default)
- Cell type classifications
- Spatial coordinates

## Output

Scripts generate plots and save them to the current directory. File formats vary by script (PNG, PDF, SVG).