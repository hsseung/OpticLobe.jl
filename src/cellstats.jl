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
# # FlyWire cell stats
# Reads cell stats from codex download, convert from nm to um
# Create arrays of same size as `ind2id`:
#
# - `cell_length` -
# - `cell_area` -
# - `cell_volume` -

# Values are `missing` for nine cells.

# %%
using CSV, DataFrames
using NamedArrays

# %%
# For standalone execution, include dependencies
if (@__MODULE__) == Main
    include("cellids.jl")
end

# %% [markdown]
# ### root IDs from Codex download
# `Int64` is compatible with the largest cell ID, even though root IDs are technically `UInt64`
# This is the official list of proofread cells in v783, and includes cells that have no connections.

# %%
cell_stats_df = CSV.read(datadep"Codex cell stats/cell_stats.csv.gz", DataFrame)

# %%
"""
    cell_length::NamedArray{Union{Missing, Float32}, 1}

Cable length of each cell in micrometers (μm).

Maps cell IDs to the total cable length of dendrites and axon for each neuron.
Converted from nanometers in the original FlyWire data.

# Examples
```julia
cell_length[1]                              # Length of first cell in μm
cell_length[Name(720575940599333574)]       # Length of specific cell by ID
mean(skipmissing(cell_length))              # Average cell length across all cells
```

# Notes  
- Values are `missing` for ~9 cells that lack morphology data
- Units: micrometers (μm), converted from nanometers (nm) in source data
- NamedArray with FlyWire root IDs as names for direct indexing
"""
cell_length_vec = missings(Float32, size(ind2id))
cell_length_vec[id2ind.(cell_stats_df.root_id)] = cell_stats_df.length_nm/1e3
cell_length = NamedArray(cell_length_vec, ind2id)

# %%
"""
    cell_area::NamedArray{Union{Missing, Float32}, 1}

Surface area of each cell in square micrometers (μm²).

Maps cell IDs to the total surface area of dendrites and axon for each neuron.
Converted from square nanometers in the original FlyWire data.

# Examples
```julia
cell_area[1]                            # Surface area of first cell in μm²
cell_area[Name(720575940599333574)]     # Surface area of specific cell by ID
mean(skipmissing(cell_area))            # Average cell surface area
```

# Notes
- Values are `missing` for ~9 cells that lack morphology data  
- Units: square micrometers (μm²), converted from nm² in source data
- NamedArray with FlyWire root IDs as names for direct indexing
"""
cell_area_vec = missings(Float32, size(ind2id))
cell_area_vec[id2ind.(cell_stats_df.root_id)] = cell_stats_df.area_nm/1e6
cell_area = NamedArray(cell_area_vec, ind2id)

# %%
"""
    cell_volume::NamedArray{Union{Missing, Float32}, 1}

Volume of each cell in cubic micrometers (μm³).

Maps cell IDs to the total volume of dendrites and axon for each neuron.
Converted from cubic nanometers in the original FlyWire data.

# Examples
```julia
cell_volume[1]                          # Volume of first cell in μm³
cell_volume[Name(720575940599333574)]   # Volume of specific cell by ID
mean(skipmissing(cell_volume))          # Average cell volume
```

# Notes
- Values are `missing` for ~9 cells that lack morphology data
- Units: cubic micrometers (μm³), converted from nm³ in source data  
- NamedArray with FlyWire root IDs as names for direct indexing
"""
cell_volume_vec = missings(Float32, size(ind2id))
cell_volume_vec[id2ind.(cell_stats_df.root_id)] = cell_stats_df.size_nm/1e9
cell_volume = NamedArray(cell_volume_vec, ind2id)
