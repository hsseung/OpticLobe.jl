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

# %% [markdown]
# ### root IDs from Codex download
# `Int64` is compatible with the largest cell ID, even though root IDs are technically `UInt64`
# This is the official list of proofread cells in v783, and includes cells that have no connections.

# %%
cell_stats_df = CSV.read(datadep"Codex cell stats/cell_stats.csv.gz", DataFrame)

# %%
cell_length = missings(Float32, size(ind2id))
cell_length[id2ind.(cell_stats_df.root_id)] = cell_stats_df.length_nm/1e3

# %%
cell_area = missings(Float32, size(ind2id))
cell_area[id2ind.(cell_stats_df.root_id)] = cell_stats_df.area_nm/1e6

# %%
cell_volume = missings(Float32, size(ind2id))
cell_volume[id2ind.(cell_stats_df.root_id)] = cell_stats_df.size_nm/1e9
