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
#     display_name: Julia 1.10.6
#     language: julia
#     name: julia-1.10
# ---

# %% [markdown]
# # Visualizations of the mapping of cells to hexagonal lattice
# Visualize cells of a given type with constant `p` or constant `q`.
#
# Creates clickable neuroglancer links.
#
# Used to generate visualizations in Supplementary Data 2, which demonstrate mapping of types to hexagonal lattice.

# %%
using OpticLobe

# %%
using Missings, MissingsAsFalse

# %%
pmax, qmax = size(pq2column)

# %% [markdown]
# ## visualize L1

# %%
celltype = "L1"  # substitute your desired cell type here

# %%
ids = type2ids(celltype)
coords = id2pq.(ids)

# %%
@mfalse for p = 1:pmax
    ng_hyper(ids[passmissing(getindex).(coords, 1) .== p], anchor = "p=$p") |> display
end

# %%
@mfalse for q = 1:qmax
    ng_hyper(ids[passmissing(getindex).(coords, 2) .== q], anchor = "q=$q") |> display
end

# %% [markdown]
# overall, Hungarian assignment looks really good for finding correct locations
#
# big gap in p=15
#
# 720575940605423532 in two columns? p=26 not right. should be p=25
#
# anterior three cells of p=30 questionable
#
# p=31 anterior cell also seems out of place
