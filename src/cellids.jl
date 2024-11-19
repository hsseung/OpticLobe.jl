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
# # Export FlyWire cell IDs to JLD2 file
# Reads cell IDs from codex download and creates:
#
# - `ind2id` - `Int64` cell IDs in canonical order corresponding to indices of weight matrix
# - `id2ind` - dictionary for looking up `Int32` index corresponding to a cell ID.
#
# These variables are written to `jldfilename`, which should be specified prior to running this script
# If not specified, this defaults to `cellids.jld2` in the working directory.

# %%
using CSV, DataFrames

# %% [markdown]
# ### root IDs from Codex download
# `Int64` is compatible with the largest cell ID, even though root IDs are technically `UInt64`
# This is the official list of proofread cells in v783, and includes cells that have no connections.

# %%
neurons = CSV.read(datadep"Codex neuron IDs/neurons.csv.gz", DataFrame)

ind2id = neurons.root_id

# %%
# ### dictionary for reverse lookup
id2ind = Dict(ind2id .=> Int32.(1:length(ind2id)))
