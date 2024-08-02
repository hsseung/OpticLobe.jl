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
# read `cellids.jld2`
#
# read Codex downloads
# - `consolidated_cell_types.csv.gz` for optic lobe types (formerly `visual_neuron_types.csv.gz`)
# - `classification.csv.gz` for central brain types
#
# compute sorted lists of type names, all with "natural" sort
# - `intrinsictypes` - intrinsic to optic lobe
# - `boundarytypes` - VPN and VCN

# - `centraltypes` - central brain`
#
# compute assignments of cells to types
# - `ind2type` - each element of this vector is a celltype or `missing`
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
#ind2id, id2ind = load("cellids.jld2", "ind2id", "id2ind")
#(id2ind::Dict)(k) = id2ind[k]

# %%
visual = CSV.read(datadep"Codex visual neuron types/visual_neuron_types.csv.gz", DataFrame)

# %%
# note that `side` annotation in this file is "corrected" relative to Jefferis annotation based on soma side
println(sum(visual.side .!= "right"), " cells not on right side")

# %% [markdown]
# ## special characters

# %%
visual.type[visual.type .== "TmY9q__perp"] .= "TmY9q⊥"

# %% [markdown]
# ## create `ind2type`
# This includes intrinsic and boundary types as defined in Codex.\
# The Codex download currently includes only type annotations for cells on the right side.

# %%
ind2type = Vector{Union{Missing, String}}(missing, length(ind2id))
ind2type[id2ind.(visual.root_id)] .= visual.type

# %% [markdown]
# ## sorted lists of intrinsic and boundary types

# %%
intrinsictypes = unique(visual.type[visual.category .== "OL intrinsic"])
intrinsictypes = sort(String.(intrinsictypes), lt = natural)

# %%
boundarytypes = unique(visual.type[visual.category .== "OL boundary"])
boundarytypes = sort(String.(boundarytypes), lt = natural)

# %%
println(length(intrinsictypes), " intrinsic and ", length(boundarytypes), " boundary types")

# %% [markdown]
# ## weirdos
# developmental accidents
# reset type to `missing`

# %%
weirdos = ["LM102a_L234-M89", "Li124_56", "SPm101_78-78"]

# %%
println("eliminating weirdo types from `intrinsictypes`")
filter!(e -> e∉weirdos, intrinsictypes)  # remove weirdo types

println("eliminating weirdo types from `ind2type`")
for w in weirdos
    @mfalse ind2type[ind2type .== w] .== missing
end

# %% [markdown]
# ## central brain types

# %%
classification = sort(CSV.read(datadep"Codex classification/classification.csv.gz", DataFrame))
@assert classification.root_id == ind2id

# %%
## old version
# classification.hemibrain_type[.~ismissing.(ind2type)] .= missing
# tochange = .~ismissing.(classification.hemibrain_type) .&& classification.super_class .== "central"

# %%
# We assume that rows of `classification` have the canonical order, the same order as `neurons`
tochange = ismissing.(ind2type) .&& classification.super_class .== "central"  # this could include some cells that are missing type annotations altogether, but that's OK
ind2type[tochange] = classification.hemibrain_type[tochange]

# %%
centraltypes = setdiff(levels(ind2type), union(intrinsictypes, boundarytypes))
centraltypes = sort(centraltypes, lt=natural)
println(length(intrinsictypes), " intrinsic, ", length(boundarytypes), " boundary, and ", length(centraltypes), " central brain types")

# %% [markdown]
# ## assignment matrix assigns cell to type

# %%
visualtypes = vcat(intrinsictypes, boundarytypes)
alltypes = vcat(visualtypes, centraltypes)
A = zeros(Bool, length(ind2type), length(alltypes))
@mfalse for (i, celltype) in enumerate(alltypes)
    A[ind2type .== celltype, i]  .= 1
end 
A = SparseMatrixCSC{Bool, Int32}(A);
