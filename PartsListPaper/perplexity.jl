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
# # graph perplexity

# %%
using OpticLobe

# %%
using SparseArrays, NamedArrays

# %%
using StatsBase

# %%
using Plots, Measures

# %%
include("config.jl")

# %% [markdown]
# ## entropy (intrinsic)

# %%
outentropy = [entropy(outfraction[celltype, :]) for celltype in intrinsictypes]
outentropy = NamedArray(outentropy, intrinsictypes)
#cellnumbers["R1-6"] = round(cellnumbers["R1-6"]/6)
outentropy = outentropy[intrinsictypes .!= "R1-6"]   # leave out, because it has 6x the number of one type

# %%
inentropy = [entropy(infraction[:, celltype]) for celltype in intrinsictypes]
inentropy = NamedArray(inentropy, intrinsictypes)
#cellnumbers["R1-6"] = round(cellnumbers["R1-6"]/6)
inentropy = inentropy[intrinsictypes .!= "R1-6"]   # leave out, because it has 6x the number of one type

# %%
# r = [1:57, 58:114, 115:171, 172:228]
r = [1:46, 47:92, 93:138, 139:184, 185:229]
ncol = length(r)
default(; # Plots defaults
    fontfamily="Helvetica",
    label="" # only explicit legend entries
    )

# %% [markdown]
# ## sort by average of in- and out-entropy

# %%
perm = sortperm(inentropy + outentropy)
pavg = [plot(
        [exp.(inentropy[perm[r[i]]].array), exp.(outentropy[perm[r[i]]].array)],
        permute = (:y, :x),
        xlim = (0.5, length(r[i])+0.5), 
        xticks=((1:length(r[i])), names(inentropy)[1][perm[r[i]]]),
        ylabel = "perplexity",
        seriestype = :stepmid,
        yrotation = 45,
        labels = (i == 1 ? ["in" "out"] : false),
        legend = :bottomright
    ) for i = ncol:-1:1]

plot(pavg..., 
    size = (200*ncol, 800*4/ncol), 
    layout = (1, ncol), 
    legendfontsize = 11, 
    bottom_margin = 8mm
)

# %%
#savefig(joinpath(TARGETDIR, "PerplexityAvgSorted.svg"))

# %% [markdown]
# ## sort by in-entropy

# %%
perm = sortperm(inentropy)
pin = [plot(
        [exp.(inentropy[perm[r[i]]].array), exp.(outentropy[perm[r[i]]].array)],
        permute = (:y, :x),
        xlim = (0.5, length(r[i])+0.5), 
        xticks=((1:length(r[i])), names(inentropy)[1][perm[r[i]]]),
        ylabel = "perplexity",
        seriestype = :stepmid,
        yrotation = 45,
        labels = (i == 1 ? ["in" "out"] : false),
        legend = :bottomright
    ) for i = ncol:-1:1]

plot(pin..., 
    size = (200*ncol, 800*4/ncol), 
    layout = (1, ncol), 
    legendfontsize = 11, 
    bottom_margin = 8mm
)

# %%
#savefig(joinpath(TARGETDIR, "PerplexityInSorted.svg"))

# %% [markdown]
# ## sort by out-entropy

# %%
perm = sortperm(outentropy)
pout = [plot(
        [exp.(inentropy[perm[r[i]]].array), exp.(outentropy[perm[r[i]]].array)],
        permute = (:y, :x),
        xlim = (0.5, length(r[i])+0.5), 
        xticks=((1:length(r[i])), names(inentropy)[1][perm[r[i]]]),
        ylabel = "perplexity",
        seriestype = :stepmid,
        yrotation = 45,
        labels = (i == 1 ? ["in" "out"] : false),
        legend = :bottomright
    ) for i = ncol:-1:1]

plot(pout..., 
    size = (200*ncol, 800*4/ncol), 
    layout = (1, ncol), 
    legendfontsize = 11, 
    bottom_margin = 8mm
)

# %%
#savefig(joinpath(TARGETDIR, "PerplexityOutSorted.svg"))

# %% [markdown]
# ## out vs in perplexity

# %%
outvsin = scatter(exp.(inentropy), exp.(outentropy), 
    xlabel = "in perplexity", ylabel = "out perplexity", 
#    hover = names(inentropy)[1], 
    legend = :none,
    tickfontsize = 12,
    guidefontsize = 14,
#    xlim = (0, 80), ylim = (0, 120)
)
Plots.abline!(1, 0, line=:dash)

# %%
#savefig(joinpath(TARGETDIR, "OutVsInEntropy.svg"))

# %%
plot(pavg..., outvsin,
    size = (200*ncol, 800*5/ncol + 200), 
#    layout = (1, ncol),
    layout = @layout([grid(1, ncol){0.6h}; [_ b{0.7w} _]]),
    left_margin = 5mm
)

# %%
savefig(joinpath(TARGETDIR, "Fig S7 raw.svg"))

# %% [markdown]
# ## difference: out- minus in-entropy

# %%
perm = sortperm(outentropy - inentropy)
pdiff = [plot(
        [outentropy[perm[r[i]]].array - inentropy[perm[r[i]]].array],
        permute = (:y, :x),
        xlim = (0.5, length(r[i])+0.5),
#        hover=optictypes[p2], 
        xticks=((1:length(r[i])), names(inentropy)[1][perm[r[i]]]),
        ylabel = "Δentropy",
#        ytickfontsize = 6,
        seriestype = :stepmid,
        yrotation = 45,
    ) for i = ncol:-1:1]

plot(pdiff..., size = (200*ncol, 800*5/ncol), layout = (1, ncol), bottom_margin = 8mm)

# %%
savefig(joinpath(TARGETDIR, "Fig S8 OutMinusInEntropy.svg"))
