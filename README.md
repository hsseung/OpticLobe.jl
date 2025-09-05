# OpticLobe.jl

A Julia package for exploring the neurons and connections of the *Drosophila* optic lobe, based on the FlyWire connectome (v783 proofreading).

[![Build Status](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hsseung/OpticLobe.jl/actions/workflows/CI.yml?query=branch%3Amain)

## About

This package accompanies two visual system papers in the [FlyWire paper package](https://www.nature.com/collections/hgcfafejia):

* **Matsliah, Yu, et al.**, [Neuronal parts list and wiring diagram for a visual system](https://doi.org/10.1038/s41586-024-07981-1), *Nature* 634:166-180 (2024)
* **Seung**, [Predicting visual function by interpreting a neuronal wiring diagram](https://doi.org/10.1038/s41586-024-07953-5), *Nature* 634:113-123 (2024)

The package provides easy access to:
- **740 visual cell types** (230 intrinsic + 510 boundary types)
- **Synaptic connectivity matrices** at cell and type levels  
- **Cell morphology and spatial coordinates**
- **Analysis tools and visualization utilities**

All data is automatically downloaded from archived sources and processed into convenient Julia data structures.

## Installation

```julia
julia> import Pkg
julia> Pkg.add(url = "https://github.com/hsseung/OpticLobe.jl")
julia> using OpticLobe
```

**Note**: The first time you load the package, it will automatically download data files. You may see "waiting for IO to finish" messages, which can be ignored.

## Examples

IDs of five `Tm1` cells.
``` julia
julia> first(type2ids("Tm1"), 5)
5-element Vector{Int64}:
 720575940599333574
 720575940603884512
 720575940604009062
 720575940604094240
 720575940604151520
```

Top ten target cells of `720575940599333574`, with number of synapses.

``` julia
julia> first(sort(W[Name(720575940599333574), :], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
cellid             │
───────────────────┼───
720575940620875399 │ 57
720575940629884688 │ 39
720575940616048757 │ 38
720575940615102003 │ 29
720575940611830245 │ 27
720575940622106996 │ 27
720575940619661648 │ 25
720575940629129692 │ 22
720575940639008319 │ 21
720575940626667928 │ 20
```

Cell type of `720575940620875399` (top cell in the list above).
``` julia
julia> ind2type[id2ind(720575940620875399)]
"Pm05"
```

Top ten target cell types of `720575940599333574`, with number of synapses.
``` julia
julia> first(sort(Wct[Name(720575940599333574), :], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
celltype  │
──────────┼────
Pm03      │ 131
Pm08      │  86
Pm02      │  72
Pm05      │  57
T2a       │  41
LMa5      │  29
T5a       │  29
T3        │  26
Tm21      │  26
Tm4       │  22
```

Number of synapses from cell type `Tm1` to cell type `TmY4`:
```julia
julia> Wtt["Tm1", "TmY4"]
5068
```
Top ten intrinsic cell types targeted by cell type `Tm1`:
``` julia
julia> first(sort(Wtt["Tm1", intrinsictypes], rev=true), 10)
10-element Named SparseArrays.SparseVector{Int32, Int32}
celltype  │
──────────┼──────
Pm03      │ 78694
Pm02      │ 70095
Pm08      │ 49431
T3        │ 35806
T2a       │ 35162
Pm05      │ 29965
Tm21      │ 21059
Y3        │ 20896
LMa5      │ 17871
Pm06      │ 17143
```

## Data Structures

### Cell Identification
- `ind2id`, `id2ind` - Convert between cell indices and FlyWire root IDs
- `ind2type` - Map cell indices to cell type names  
- `id2pq` - Spatial coordinates (hexagonal p,q system)

### Cell Type Classifications  
- `intrinsictypes` (230) - Neurons intrinsic to the optic lobe
- `boundarytypes` (510) - Visual projection and centrifugal neurons
- `othertypes` (7805) - All other cell types (central brain, etc.)
- `visualtypes` - Combined intrinsic + boundary types (740 total)
- `alltypes` - All cell types in the dataset

### Connectivity Matrices
- `W` - Cell-to-cell synaptic weight matrix (sparse, ~130K × 130K cells)
- `Wtt` - Type-to-type connectivity (aggregated synapses)
- `Wct` - Cell-to-type connectivity 
- `Wtc` - Type-to-cell connectivity
- `infraction`, `outfraction` - Normalized connectivity fractions
- `inmean`, `outmean` - Mean synapses per cell by type

### Cell Type Assignment
- `A` - Boolean matrix: cells × all types  
- `Ai` - Boolean matrix: cells × intrinsic types only

### Type Hierarchies
- `class2families` - Visual classes to type families mapping
- `family2types` - Type families to individual types mapping

## Useful Functions

### Analysis Functions
- `toppre(celltype, n)` - Top n presynaptic partners of a cell type
- `toppost(cellid, n)` - Top n postsynaptic partners of a cell  
- `type2ids(typename)` - Get all cell IDs belonging to a cell type
- `showall(vector)` - Display all elements of a named vector

### Visualization & Integration
- `codex_open(cellid)` - Open cell in FlyWire Codex browser
- `ng_open(cellids)` - Open cells in Neuroglancer 3D viewer
- `ng_hyper(cellids)` - Create Neuroglancer hyperlink
- `strings2ticks(strings)` - Convert strings to plot tick marks

### Spatial Analysis
- `pq2column(p, q)` - Convert hex coordinates to column ID
- Hexagonal visualization tools in `hexgraphics` module

## Analysis Examples

The `PartsListPaper/` directory contains complete analysis scripts that reproduce the figures from the Nature papers. Each script has its own environment with all required dependencies.

```bash
cd PartsListPaper
julia --project=. "Fig 1 cell numbers.jl"    # Reproduce Figure 1
julia --project=. "dendrograms.jl"            # Hierarchical clustering
```

See `PartsListPaper/README.md` for detailed usage instructions.

## Data Sources

All data is automatically downloaded via DataDeps.jl from archived sources:
- **FlyWire Codex data** (v783 proofreading) archived on Zenodo

Data files are cached locally after first download.

