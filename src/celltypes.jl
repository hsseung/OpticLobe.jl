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
# if you want to run as a standalone script, first run these scripts:
#   include("codexdependencies.jl")
#   include("cellids.jl")
#
# and define this function
#   (id2ind::Dict)(k) = id2ind[k]
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
using Missings # for `levels`

# %%
# for the assignment matrix A.
# A is *not* a NamedArray. This conversion is made later, for the convenience of human user
using SparseArrays

# %%
import StatsBase.countmap

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
ind2category = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2category[id2ind.(visual.root_id)] .= visual.category

# %%
ind2side = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2side[id2ind.(visual.root_id)] .= visual.side

# %%
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
A = zeros(Bool, length(ind2type), length(alltypes))
@mfalse for (i, celltype) in enumerate(alltypes)
    A[ind2type .== celltype, i]  .= 1
end 
A = SparseMatrixCSC{Bool, Int32}(A);
