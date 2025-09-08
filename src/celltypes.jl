# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Cell Types including visual and central brain
#
# read Codex downloads
# - `consolidated_cell_types.csv.gz` for all cell type annotations
# - `visual_neuron_types.csv.gz` for intrinsic vs boundary categories, and side
#
# compute sorted lists of type names, all with "natural" sort
# - `intrinsictypes` - intrinsic to optic lobe
# - `boundarytypes` - VPN and VCN, bilateral?

# - `othertypes` - other cell types
#
# compute assignments of cells to types
# - `ind2type` - each element of this vector is a cell type or `missing`
# - `A[cell, type] = true` if `cell` is assigned to `type`

# %%
using DataDeps, CSV, DataFrames

# %%
using NaturalSort

# %%
using MissingsAsFalse
import Missings.levels

# %%
# for the assignment matrix A.
# A is *not* a NamedArray. This conversion is made later, for the convenience of human user
using SparseArrays


# %%
import StatsBase.countmap

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    include("cellids.jl")
end

# %%
visual = CSV.read(datadep"Codex visual neuron types/visual_neuron_types.csv.gz", DataFrame)
consolidated = CSV.read(datadep"Codex cell types/consolidated_cell_types.csv.gz", DataFrame)

# %%
# note that `side` annotation in this file is "corrected" relative to Jefferis annotation based on soma side
countmap(visual.side)

# %% [markdown]
# ## special characters

# %%
#visual.type[visual.type .== "TmY9q__perp"] .= "TmY9q⊥"
consolidated.primary_type[consolidated.primary_type .== "TmY9q__perp"] .= "TmY9q⊥"

# %% [markdown]
# ## create `ind2type`
# This includes intrinsic and boundary types as defined in Codex.\
# The Codex download currently includes only type annotations for cells on the right side.

# %%
#ind2type = Vector{Union{Missing, String}}(missing, length(ind2id))
#ind2type[id2ind.(visual.root_id)] .= visual.type

# %%
"""
    ind2category::Vector{Union{Missing, String}}

Maps cell indices to category classification ("intrinsic", "boundary", or `missing`).

Categories classify visual system neurons based on their anatomical location:
- "intrinsic": Neurons with processes confined to the optic lobe
- "boundary": Visual projection neurons (VPN) and visual centrifugal neurons (VCN)
- `missing`: Non-visual neurons (central brain, etc.) that lack category annotation

# Examples
```julia
ind2category[1]                    # Category of first cell
ind2category[id2ind[720575940599333574]]  # Category of specific cell
findall(ind2category .== "intrinsic")    # All intrinsic cell indices
```

# Notes
- Only visual system neurons have category annotations
- Based on FlyWire Codex visual neuron classifications
- Used to define `intrinsictypes` and `boundarytypes` vectors
"""
ind2category = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2category[id2ind.(visual.root_id)] .= visual.category

# %%
"""
    ind2side::Vector{Union{Missing, String}}

Maps cell indices to hemisphere side ("left", "right", or `missing`).

Side classification for visual system neurons based on their spatial location:
- "left": Neurons primarily in the left optic lobe
- "right": Neurons primarily in the right optic lobe  
- `missing`: Non-visual neurons or cells without clear laterality

# Examples
```julia
ind2side[1]                     # Side of first cell
ind2side[id2ind[720575940599333574]]  # Side of specific cell
findall(ind2side .== "right")   # All right hemisphere cell indices
```

# Notes
- Only visual system neurons have side annotations
- Side filtering uses this information to restrict analyses to single hemispheres
- Essential for `tracetypes`, `tracebacktypes`, and spatial analysis functions
"""
ind2side = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2side[id2ind.(visual.root_id)] .= visual.side

# %%
"""
    ind2type::Vector{Union{Missing, String}}

Maps cell indices to primary cell type names (e.g., "Tm1", "Dm3v", or `missing`).

The primary cell type classification for all neurons in the dataset:
- Visual types: "Tm1", "Dm3v", "LC10", etc. (740 types total)
- Central brain types: Thousands of additional types
- `missing`: Unclassified or untyped cells

# Examples
```julia
ind2type[1]                     # Type of first cell
ind2type[id2ind[720575940599333574]]  # Type of specific cell
findall(ind2type .== "Tm1")     # All Tm1 cell indices
```

# Notes
- Based on FlyWire Codex consolidated cell type annotations
- Foundation for all type-based analyses and connectivity computations
- Used to construct assignment matrix `A` and type vectors (`alltypes`, etc.)
"""
ind2type = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2type[id2ind.(consolidated.root_id)] .= consolidated.primary_type

# %% [markdown]
# ## sorted lists of intrinsic and boundary types

# %%
#intrinsictypes = unique(visual.type[visual.category .== "intrinsic"])
intrinsictypes = @mfalse unique(skipmissing(ind2type[ind2category .== "intrinsic"]))
intrinsictypes = sort(String.(intrinsictypes), lt = natural)

# %%
#boundarytypes = unique(visual.type[visual.category .== "boundary"])
boundarytypes = @mfalse unique(skipmissing(ind2type[ind2category .== "boundary"]))
boundarytypes = sort(String.(boundarytypes), lt = natural)

# %%
othertypes = setdiff(levels(ind2type), union(intrinsictypes, boundarytypes))
othertypes = sort(othertypes, lt=natural)
println(length(intrinsictypes), " intrinsic, ", length(boundarytypes), " boundary, and ", length(othertypes), " other types")

# %% [markdown]
# ## weirdos
# developmental accidents
# reset type to `missing`

# %%
#weirdos = ["LM102a_L234-M89", "Li124_56", "SPm101_78-78"]
weirdos = ["LM102a_L234-M89", "Li124_56"]

# %%
println("eliminating weirdo types from `intrinsictypes`")
filter!(e -> e∉weirdos, intrinsictypes)  # remove weirdo types

println("eliminating weirdo types from `ind2type`")
for w in weirdos
    @mfalse ind2type[ind2type .== w] .= missing
end

# %% [markdown]
# ## central brain types

# %%
#classification = sort(CSV.read(datadep"Codex classification/classification.csv.gz", DataFrame))

# Verify that rows of `classification` have the canonical order, the same order as `neurons`
#@assert classification.root_id == ind2id

# %% [markdown]
# ## assignment matrix assigns cell to type

# %%
visualtypes = vcat(intrinsictypes, boundarytypes)
alltypes = vcat(visualtypes, othertypes)

println("computing assignment matrix: all")
"""
    A::SparseMatrixCSC{Bool, Int32}

The assignment matrix `A[c, t]` assigns cell `c` to type `t`.

# Indexing
- `c` can be either a cell index (1, 2, 3, ...) or cell ID (FlyWire root ID)
- `t` is one of the strings in `alltypes`, or the corresponding type index

# Returns
- `true` if cell `c` is assigned to type `t`, `false` otherwise

# Examples
```julia
A[1, 1]                    # Check if first cell belongs to first type
A[Name(720575940599333574), Name("Tm1")]  # Check if specific cell is Tm1
A[:, Name("Tm1")]          # All cells of type Tm1 (boolean vector)
```

# Notes
- Sparse matrix for memory efficiency (~130K cells × ~8K types)
- Use `NamedArrays.Name()` for string-based indexing
- Dimensions: `length(ind2id) × length(alltypes)`
"""
A = zeros(Bool, length(ind2type), length(alltypes))
@mfalse for (i, celltype) in enumerate(alltypes)
    A[ind2type .== celltype, i]  .= 1
end
A = SparseMatrixCSC{Bool, Int32}(A);

# %%
# for typing the right optic lobe, we eliminate cells in the left optic lobe
# i.e. A_right includes right optic lobe, bilateral optic lobes, central brain, etc.
# println("computing assignment matrix: left, right")
# @mfalse A_right = sparse(ind2side .!= "left") .* A;
# @mfalse A_left = sparse(ind2side .!= "right") .* A;


"""
    type2ids(celltype::String; side::String="right") -> Vector{Int64}

Get cell IDs for all cells of a specified type, with optional side filtering for visual types.

For visual cell types (those in `visualtypes`), applies side filtering to return only cells
from the specified hemisphere. For non-visual cell types (central brain, etc.), returns all
cells regardless of side since ind2side doesn't yet contain information about non-visual types.

# Arguments
- `celltype::String`: Name of the cell type
- `side::String`: Side filter for visual types ("left", "right"). Default "right".
  Ignored for non-visual cell types.

# Returns
- `Vector{Int64}`: Cell IDs of the specified type and side (for visual types)

# Examples
```julia
# Visual cell type - applies side filtering
tm1_right = type2ids("Tm1")                    # right side (default)
tm1_left = type2ids("Tm1", side="left")        # left side

# Non-visual cell type - ignores side parameter
central_cells = type2ids("SomeCentralType")     # returns all cells
```

# Notes
- Uses global variables `visualtypes`, `ind2id`, `ind2type`, `ind2side`
- Visual vs non-visual determination based on membership in `visualtypes`
- TODO: Extend side filtering to work for non-visual types that have side information
"""
function type2ids(celltype::String; side::String="right")
    if celltype in visualtypes
        @mfalse ind2id[(ind2type .== celltype) .& (ind2side .== side)]
    else
        @mfalse ind2id[ind2type .== celltype]
    end
end
