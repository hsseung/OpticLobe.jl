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
# Downloads synapses from Codex and creates:
#
# - `W[pre, post]` sparse matrix, number of synapses (`Int32`) from `pre` to `post` (`Int32`)
#
# Depends on `id2ind` defined in `cellids.jl`.
#
# Note that `W` does not depend on cell types, except that outgoing synapses from T1 cells are zeroed out.
#
# `W` is *not* yet a `NamedArray`. This conversion is made later, as it is for the convenience of human user.

# %%
using Preferences

# For standalone execution, include dependencies and set defaults
if (@__MODULE__) == Main
    include("celltypes.jl")
    
    # Standalone defaults
    const load_both = false
    const default_synapses = "Buhmann"
else
    # Load preferences from package
    const load_both = @load_preference("load_both_synapses", false)
    const default_synapses = @load_preference("default_synapses", "Princeton")
end

# %%
using DataDeps, CSV, DataFrames
using SparseArrays, NamedArrays
using MissingsAsFalse
using LinearAlgebra

# %% [markdown]
# ## download synapse table from Codex into DataFrame

# %%
function load_synapses(synapse_type::String)
    println("reading synapse table ($(synapse_type) version)")
    if synapse_type == "Buhmann"
        df = CSV.read(datadep"Codex Buhmann connections no threshold/connections_buhmann_no_threshold.csv.gz", DataFrame)
    elseif synapse_type == "Princeton"
        df = CSV.read(datadep"Codex Princeton connections no threshold/connections_princeton_no_threshold.csv.gz", DataFrame)
    else
        error("Unknown synapses version: $(synapse_type)")
    end
    return df
end

# %% [markdown]
# ## function to build weight matrix from synapse data

# %%
function build_weight_matrix(synapse_type::String)
    df = load_synapses(synapse_type)
    
    # use grouped data frame to combine synapses that share the same pre- and post-synaptic cells
    println("counting synapses for each pre-post pair of neurons")
    gdf = groupby(df, [:pre_root_id, :post_root_id])
    weighttable = combine(gdf, :syn_count => sum)
    
    # convert data frame to sparse matrix
    println("converting list of connections to sparse matrix")
    
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
    
    # eliminate autapses (self-connections)
    W[diagind(W)] .= 0
    
    println("total number of synapses = ", sum(W))
    if synapse_type == "Buhmann"
        println("This number should match 54492922 for Buhmann v783 Feb 2024")
    elseif synapse_type == "Princeton"
        println("This number should match 76944499 for Princeton synapses")
    end
    
    return W
end

# Load both versions if requested, otherwise just the default
if load_both
    println("Loading both synapse versions...")
    W_Buhmann = build_weight_matrix("Buhmann")
    W_Princeton = build_weight_matrix("Princeton")
    W = (default_synapses == "Buhmann") ? W_Buhmann : W_Princeton
else
    W = build_weight_matrix(default_synapses)
end

# %% [markdown]
# ## use prior info to eliminate T1 outgoing synapses
# beware as this can lead to NaNs in normalized connection strengths

# %%
println("zeroing outgoing synapses of T1")
@mfalse W[ind2type .== "T1", :] .= 0  # manual adjustment as T1 is supposed to lack outgoing synapses
