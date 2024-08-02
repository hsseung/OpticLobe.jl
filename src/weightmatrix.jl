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
#     display_name: Julia 1.10.4
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Export FlyWire weight matrix to JLD2 file
# Downloads synapses from Codex and creates:
#
# - `W[pre, post]` sparse matrix, number of synapses (`Int32`) from `pre` to `post` (`Int32`)
#
# Depends on `id2ind` from `cellids.jld2`.
#
# Note that `W` does not depend on cell types, except that outgoing synapses from T1 cells are zeroed out.
#
# `W` is *not* yet a `NamedArray`. This conversion is made later, as it is for the convenience of human user.

# %%
using DataDeps, CSV, DataFrames
using SparseArrays, NamedArrays
using MissingsAsFalse

# %% [markdown]
# ## download synapse table from Codex into DataFrame

# %%
println("reading synapse table")
df = CSV.read(datadep"Codex connections no threshold/connections_no_threshold.csv.gz", DataFrame)

# %% [markdown]
# ## use grouped data frame to combine synapses that share the same pre- and post-synaptic cells

# %%
println("counting synapses for each pre-post pair of neurons")
gdf = groupby(df, [:pre_root_id, :post_root_id])
weighttable = combine(gdf, :syn_count => sum)

# %% [markdown]
# ## convert data frame to sparse matrix

# %%
println("converting list of connections to sparse matrix")

# id2ind, ind2id = load("cellids.jld2", "id2ind", "ind2id")   # if NamedArray is desired
# id2ind = load("cellids.jld2", "id2ind")

II = Int32[]     # Int32 sufficient here as there are only 10^5 neurons
JJ = Int32[]
WW = Int32[]
for (idpre, idpost, count) in zip(weighttable.pre_root_id, weighttable.post_root_id, weighttable.syn_count_sum)
    if haskey(id2ind, idpre) && haskey(id2ind, idpost)
        push!(II, id2ind[idpre])
        push!(JJ, id2ind[idpost])
        push!(WW, count)
    end
end

ncell = length(id2ind)
W = sparse(II, JJ, WW, ncell, ncell)

# %%
# this number was 54519562 in Oct 2023 file
# now updated for Feb 2024 file
println("total number of synapses = ", sum(W))
println("This number should match 54492922 for v783 Feb 2024")

# %% [markdown]
# ## use prior info to eliminate T1 outgoing synapses
# beware as this can lead to NaNs in normalized connection strengths

# %%
println("zeroing outgoing synapses of T1")
@mfalse W[ind2type .== "T1", :] .= 0  # manual adjustment as T1 is supposed to lack outgoing synapses
