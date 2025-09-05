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

# %%
visualcolumns = CSV.read(datadep"Codex visual neuron columns/column_assignment.csv.gz", DataFrame)

# id2pq is defined for both eyes
# The offset is backward compatible with right eye, and makes 1 the minimum value of p and q for the right eye
# But 0 is the minimum value for the left eye. This should be fixed later.
pq = [visualcolumns.p .+ 19 visualcolumns.q .+ 17]
id2pq = Dict(visualcolumns.root_id .=> map(collect, eachrow(pq)))

# Note that pq2column is only defined at present for the right eye
# Create OffsetArray for right hemisphere columns indexed by p and q
using OffsetArrays
right_rows = visualcolumns[visualcolumns.hemisphere .== "right", :]
p_range = minimum(right_rows.p):maximum(right_rows.p)
q_range = minimum(right_rows.q):maximum(right_rows.q)
pq2column = OffsetArray(Matrix{Union{Missing, Int}}(missing, length(p_range), length(q_range)), p_range, q_range)

for row in eachrow(right_rows)
    pq2column[row.p, row.q] = row.column_id
end

pq2column = parent(pq2column)
