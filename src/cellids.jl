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
# Reads cell IDs from codex download and creates:
#
# - `ind2id` - `Int64` cell IDs in canonical order corresponding to indices of weight matrix
# - `id2ind` - dictionary for looking up `Int32` index corresponding to a cell ID.
#
# These variables are written to `jldfilename`, which should be specified prior to running this script
# If not specified, this defaults to `cellids.jld2` in the working directory.

# %%
using CSV, DataFrames

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    println("Running in standalone mode")

    include("codexdependencies.jl")
end

# %% [markdown]
# ### root IDs from Codex download
# `Int64` is compatible with the largest cell ID, even though root IDs are technically `UInt64`
# This is the official list of proofread cells in v783, and includes cells that have no connections.

# %%
neurons = CSV.read(datadep"Codex neuron IDs/neurons.csv.gz", DataFrame)

"""
    ind2id::Vector{Int64}

Vector of FlyWire root IDs in canonical order corresponding to matrix indices.

Maps cell indices (1, 2, 3, ...) to FlyWire root IDs (64-bit integers like 720575940599333574).
This vector defines the canonical ordering used throughout OpticLobe for matrix operations.

# Examples
```julia
ind2id[1]           # Get FlyWire root ID for first cell
ind2id[1:5]         # Get first 5 cell IDs
length(ind2id)      # Total number of cells (~130K)
```

# Notes
- Length corresponds to total number of proofread cells in FlyWire v783
- Used as row/column names for connectivity matrices when converted to NamedArrays
- Canonical ordering ensures consistency across all OpticLobe data structures
"""
ind2id = neurons.root_id

# %%
# ### dictionary for reverse lookup
"""
    id2ind::Dict{Int64, Int32}

Dictionary for reverse lookup: maps FlyWire root IDs to cell indices.

Converts FlyWire root IDs (64-bit integers) back to sequential indices (1, 2, 3, ...) 
used for array operations and matrix indexing.

# Examples  
```julia
id2ind[720575940599333574]     # Get index for specific cell ID
id2ind.(some_cell_ids)         # Broadcast to convert vector of IDs to indices
id2ind(missing_id)             # Returns missing for unknown IDs (see function method)
```

# Notes
- Returns `Int32` indices for memory efficiency
- Function method `id2ind(key)` returns `missing` for unknown keys (safer than direct dict access)
- Inverse of `ind2id`: `id2ind[ind2id[i]] == i` for valid indices
"""
id2ind = Dict(ind2id .=> Int32.(1:length(ind2id)))

# %%
# use dictionary as function
# useful when broadcasting e.g. id2ind
function (d::Dict)(key)
    get(d, key, missing)
end
