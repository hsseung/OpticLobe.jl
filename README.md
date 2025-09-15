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

All data is automatically downloaded from the FlyWire Codex and processed into convenient Julia data structures.

## Installation

```julia
julia> import Pkg
julia> Pkg.add(url = "https://github.com/hsseung/OpticLobe.jl")
julia> using OpticLobe
```

**Note**: The first time you load the package, it will automatically download ~350MB of data files. You may see "waiting for IO to finish" messages, which can be ignored.

## Synapse Versions

The package supports two synapse detection methods:
- **Princeton** (default): Latest synapse predictions with improved accuracy
- **Buhmann**: Original predictions from Buhmann et al. (2021)

```julia
# Switch to Buhmann synapses
set_default_synapses("Buhmann")    # Restart Julia to take effect

# Load both versions simultaneously (uses more memory)  
enable_both_synapses(true)         # Restart Julia to take effect
```

After enabling both versions, access them as follows:
```julia
W_Princeton[cellid1, cellid2]  # Princeton synapses
W_Buhmann[cellid1, cellid2]    # Buhmann synapses
W[cellid1, cellid2]             # Default version (Princeton unless changed)
```

### Checking Current Configuration
You can check which synapse version is currently active:
```julia
# Using internal variables (quick check)
OpticLobe.default_synapses  # "Princeton" or "Buhmann"
OpticLobe.load_both         # true if both versions loaded, false otherwise

# Using Preferences.jl (standard approach)
using Preferences
load_preference(OpticLobe, "default_synapses")    # "Princeton" or "Buhmann" 
load_preference(OpticLobe, "load_both_synapses")  # true or false
```

## Examples

**Note**: The examples below show approximate output. Numbers may vary due to updated synapse predictions and cell type annotations in newer versions.

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
OpticLobe uses two different cell identifiers:

- **Cell ID**: FlyWire root IDs (64-bit integers like `720575940599333574`) that uniquely identify cells in the FlyWire dataset
- **Cell Index**: Sequential integers (1, 2, 3, ...) used internally for array indexing and matrix operations

Conversion functions:
- `ind2id`, `id2ind` - Convert between cell indices and FlyWire root IDs
- `ind2type` - Map cell indices to cell type names  
- `id2pq` - Map cell IDs to spatial coordinates (hexagonal p,q system)

### Cell Type Classifications  
These are vectors of strings containing cell type names:

- `intrinsictypes` (230) - Neurons intrinsic to the optic lobe
- `boundarytypes` (510) - Visual projection and centrifugal neurons
- `othertypes` (7805) - All other cell types (central brain, etc.)
- `visualtypes` - Combined intrinsic + boundary types (740 total)
- `alltypes` - All cell types in the dataset

### Connectivity Matrices

All connectivity matrices are `NamedArray` objects that can be indexed by either position or name using `Name()`. For example, `W[Name(720575940599333574), Name(720575940620875399)]` or `Wtt["Tm1", "Dm3v"]`. The core connectivity data is stored as synapse counts (`Int32` values) representing the number of synaptic connections between neurons.

**Cell-to-cell connectivity:**
- `W` - Full synaptic weight matrix (sparse, ~130K × 130K cells). `W[i, j]` is the number of synapses from neuron `i` to neuron `j`.

**Cell-to-type and type-to-cell connectivity:**
- `Wct` - Cell-to-type connectivity: `Wct[c, t]` gives total synapses from cell `c` to type `t`
- `Wtc` - Type-to-cell connectivity: `Wtc[t, c]` gives total synapses from type `t` to cell `c`

**Type-to-type connectivity:**
Matrices that summarize connections between cell types:
- `Wtt` - Type-to-type connectivity: `Wtt[pretype, posttype]` is raw synapse counts.
- `infraction`, `outfraction` - Normalized versions of `Wtt` (0-1 scale).
- `inmean`, `outmean` - Alternative normalization of `Wtt` giving mean synapses per cell of a type.
- `importance` - Connection importance: `max(infraction, outfraction)` for each connection.
- `inrank`, `outrank` - Rank ordering of connections by strength (intrinsic types only).

### Cell Type Assignment

Boolean matrices that encode which cells belong to which types. These sparse matrices enable efficient filtering and selection of cells by type membership:

- `A` - Boolean matrix assigning cells to types. where `A[c, t]` is `true` if cell `c` belongs to type `t`.
- `Ai` - Submatrix of `A` containing only intrinsic types (cells × intrinsic types).

### Type Hierarchies

The visual system types are organized in a three-level hierarchy: **classes** → **families** → **individual types**. This hierarchical organization reflects functional and anatomical relationships:

- `class2families` - Maps 5 visual classes ("receptor", "columnar", "interneuron", etc.) to their constituent type families
- `family2types` - Maps type families ("Tm", "Dm", "R7-8", etc.) to individual cell types within each family

## Useful Functions

### Analysis Functions

Basic functions for exploring connectivity patterns and cell properties:
- `toppre(celltype; nresults=15, sort=:in, values=:fraction)` - Top presynaptic partners of a cell type
- `toppost(celltype; nresults=15, sort=:out, values=:fraction)` - Top postsynaptic partners of a cell type  
- `type2ids(typename; side="right")` - Get all cell IDs belonging to a cell type
- `showall(vector)` - Display all elements of a named vector
- `inmaps(paths)` - Convert connectivity matrix to spatial input maps (receptive fields)
- `outmaps(paths)` - Convert connectivity matrix to spatial output maps (projective fields)

### Pathway Analysis (with Side Filtering)

Functions for analyzing multi-step synaptic pathways with hemisphere-specific filtering. All functions support a `side` parameter (default: "right") to restrict analysis to cells in a specific hemisphere, which is essential for studying lateralized visual processing:

- `tracetypes(celltypes; side="right")` - Multi-step connectivity through cell type sequence
- `tracebacktypes(celltypes; side="right")` - Like tracetypes but with normalized connectivity
- `preimage(pretype, posttype; side="right")` - Spatial map of presynaptic inputs by type  
- `prepreimage(prepretype, pretype, posttype; side="right")` - Two-step presynaptic maps

### Visualization & Integration

Functions for viewing cells in external tools and preparing data for plotting. The FlyWire integration functions automatically detect your environment (Jupyter vs REPL) and adapt their behavior accordingly:

- `codex_open(cellids; version=783)` - Open cells in FlyWire Codex browser (works in both Jupyter and REPL)
- `ng_open(cellids; version=783)` - Open cells in Neuroglancer 3D viewer (works in both Jupyter and REPL)
- `ng_hyper(cellids; anchor, version=783)` - Create clickable Neuroglancer hyperlinks (Jupyter) or print URLs (REPL)
- `strings2ticks(strings)` - Convert string vector to (positions, labels) for plot tick marks

### Spatial Analysis

Tools for working with the hexagonal column coordinate system of the optic lobe and creating spatial visualizations. The optic lobe uses a hexagonal lattice where each column has (p,q) coordinates:

- `id2pq[cellid]` - Map cell IDs to hexagonal (p,q) column coordinates 
- `pq2column[p, q]` - Map hex coordinates to column IDs (right hemisphere only)
- `crop(image, center, radius; sat)` - Extract square region from image centered at coordinates
- `square2hex(square)` - Convert square array to hexagonal shape by masking corners
- `rect2hex(rect; hexelsize, pointannotation)` - Display rectangular array on hexagonal lattice
- `montage(images; hexelsize, labels, ellipses, ...)` - Create grid layout of hexagonal eye maps
- `triad(celltype)` - Generate spatial triad visualization for cell type
- Hexagonal visualization tools in `hexgraphics` module

### Cell Properties

Morphological measurements for individual cells. These are NamedArrays indexed directly by FlyWire cell IDs, with measurements converted from nanometers to micrometers for convenience:

- `cell_length[cellid]` - Cable length of cell in micrometers
- `cell_area[cellid]` - Surface area of cell in square micrometers  
- `cell_volume[cellid]` - Volume of cell in cubic micrometers

### Cell Classification & Mapping

Functions for mapping between cell indices (not cell IDs) and various classification schemes. These are vectors indexed by sequential cell indices, often returning `missing` for cells without annotations:

- `ind2type[cellindex]` - Map cell indices to primary type names (vector of strings)
- `ind2category[cellindex]` - Map cell indices to category (intrinsic/boundary/missing)
- `ind2side[cellindex]` - Map cell indices to hemisphere (left/right/missing) 
- `ind2superclass[cellindex]`, `ind2class[cellindex]`, `ind2subclass[cellindex]` - Hierarchical classifications
- `ind2nt[cellindex]`, `type2nt(typename)` - Map to neurotransmitter types

### Advanced Spatial Functions

Specialized functions for creating complex spatial visualizations and analyzing spatial patterns in the hexagonal coordinate system:

- `eyetriad(cellid)` - Generate eye-centered triad for specific cell
- `typetriad(celltype)` - Generate type-centered triad visualization  
- `celltriad(cellid)` - Generate cell-centered triad visualization
- `eyehot(image; sat)` - Convert array to heatmap using hot colormap with eye masking
- `ellipsesummary(image)`, `drawellipse(ellipse)` - Ellipse fitting and visualization
- `hexproject`, `drawpqaxes`, `hexannulus` - Hexagonal projection utilities
- `HexagonEye(p, q, hexelsize)` - Create hexagon at eye coordinates

### Utility Functions

Miscellaneous helper functions for data processing and analysis:

- `convert2arrows(matrix)` - Convert connectivity matrix to arrow format
- `findcenter(coordinates)` - Find center of coordinate cluster
- `convcluster(data)` - Convex clustering analysis
- `seven` - Utility constant/function
- `Name` - NamedArrays indexing utility

### Configuration Functions
- `set_default_synapses(version)` - Set default synapse version (Princeton/Buhmann)
- `enable_both_synapses(enabled)` - Enable simultaneous loading of both synapse versions

## Analysis Examples

The repository contains complete analysis scripts that reproduce the figures from both Nature papers:

### Parts List Paper Analysis
The `PartsListPaper/` directory contains scripts for **Matsliah, Yu, et al.**, [Neuronal parts list and wiring diagram for a visual system](https://doi.org/10.1038/s41586-024-07981-1), *Nature* 634:166-180 (2024).

```bash
cd PartsListPaper
julia --project=. "Fig 1 cell numbers.jl"    # Reproduce Figure 1
julia --project=. "dendrograms.jl"            # Hierarchical clustering
```

See `PartsListPaper/README.md` for detailed usage instructions.

### Form Vision Paper Analysis  
The `FormVisionPaper/` directory contains scripts for **Seung**, [Predicting visual function by interpreting a neuronal wiring diagram](https://doi.org/10.1038/s41586-024-07953-5), *Nature* 634:113-123 (2024).

```bash
cd FormVisionPaper
julia --project=. FiguresSpatial.jl          # Spatial analysis and montages
julia --project=. FiguresNonspatial.jl       # Non-spatial connectivity analysis
```

See `FormVisionPaper/README.md` for Jupytext setup and detailed instructions.

## Data Sources

All data is automatically downloaded via DataDeps.jl from the latest versions on the FlyWire Codex.
- **Cell morphology** and spatial coordinates  
- **Synaptic connectivity** (Princeton and Buhmann predictions)
- **Cell type annotations** and hierarchical classifications

Data files are cached locally after first download (~350MB total).

