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
# # Figure S4. Type-to-type connectivity matrix

# %%
using OpticLobe

include("config.jl")

# %% jupyter={"outputs_hidden": false}
using NamedArrays, SparseArrays
using Plots
using OrderedCollections
using Measures

# %% [markdown]
# ## visualize matrix of synapse number

# %%
# reorder types by families
orderedtypes = vcat(values(family2types)...)

# %%
I, J, vals = findnz(Wtt[orderedtypes, orderedtypes].array)

# %%
dotsizes = min.(sqrt.(vals), 60)/15  # nonlinear scaling to make weaker connections visible. strong connections saturate at a maximum dot size of four

# %%
# decided that this wasn't worth it, because some dots could become visible if reader zooms in
# this threshold retains roughly half the dots, compressing the SVG file by about half
# change `I[visible], J[visible], and dotsizes[visible]` in scatter plot below
#visible = dotsizes .> 0.3
#println("fraction of dots retained = ", sum(visible)/length(vals))

# %%
default(
    fontfamily = "Helvetica",
    label = "",
)

cmap = palette(:tol_bright) # a Tol color scheme that is color blind safe

# %%
ntypes = length(orderedtypes)
odd = 1:2:ntypes
even = 2:2:ntypes

# note J and I are transposed due to matrix layout
#p = scatter(J[visible], I[visible], markersize = dotsizes[visible],
p = scatter(J, I, markersize = dotsizes,
    size = (1500, 1500), yflip = true, msc = :auto,
    yticks = (1:2:length(orderedtypes), orderedtypes[1:2:end]), tickfontcolor = cmap[2],
    xticks = (1:2:length(orderedtypes), orderedtypes[1:2:end]), 
    xrotation = 90,
    xlim = (0, ntypes) .+ 0.5,
    ylim = (0, ntypes) .+ 0.5,
    grid = false,
    ylabel = "Presynaptic",
    xlabel = "Postsynaptic",
    left_margin = 5mm,
    guidefontsize = 14
)
scatter!(twinx(), [], [], yflip=true, ylim = (0, ntypes) .+ 0.5,
    yticks = (2:2:length(orderedtypes), orderedtypes[2:2:end]), tickfontcolor = cmap[3],
)
scatter!(twiny(), [], [], xlim = (0, ntypes) .+ 0.5,
    xticks = (2:2:length(orderedtypes), orderedtypes[2:2:end]), tickfontcolor = cmap[3],
    xrotation = 90
)
hline!(p, odd, linealpha = 0.3, linecolor=cmap[2])
hline!(p, even, linealpha = 0.3, linecolor=cmap[3])
vline!(p, odd, linealpha = 0.3, linecolor=cmap[2])
vline!(p, even, linealpha = 0.3, linecolor=cmap[3])

# %%
savefig(joinpath(TARGETDIR, "Fig S4.pdf"))
