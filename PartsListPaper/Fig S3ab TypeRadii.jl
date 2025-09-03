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
# # Fig S3ab Type Radii

# %%
using OpticLobe

include("config.jl")

# %%
using SparseArrays, NamedArrays, LinearAlgebra, Distances
using StatsBase
using Plots, Measures
using MissingsAsFalse
using ProgressMeter

# %% [markdown]
# ## compute centers by heuristic minimization of average Jaccard distance

# %%
"""
find center of columns in Z by minimizing mean of jaccard distances
"""
function jaccardcenter(Z)
    c = median(Z, dims = 2)[:]
    initialcost = mean(colwise(jaccard, Z, c))
    oldcost = initialcost
    newcost = 0
    while true
        for i = 1:size(Z, 1)
            vals = sort(unique(Z[i, :]))
            C = repeat(c, 1, length(vals))
            C[i, :] = vals
            newcost, winner = findmin(mean(pairwise(jaccard, Z, C), dims=1)[:])
            c[i] = vals[winner]
        end
        if newcost == oldcost
            break
        else
            oldcost = newcost
        end
    end
    return c, newcost, initialcost
end

# %% [markdown]
# ## projection vectors (intrinsic types only)

# %%
featurenames = vcat("in-".*intrinsictypes, "out-".*intrinsictypes)
Xi = NamedArray(
    vcat(Wtc[intrinsictypes, :].array, Wct[:, intrinsictypes].array'),    # combine inputs and outputs
    names=(featurenames, ind2id)
    )

# %% [markdown]
# ## compute cluster centers: intrinsic feature vector, intrinsic types only

# %%
println("computing cluster centers")
centers = zeros(length(featurenames), length(intrinsictypes))
@showprogress for (i, celltype) in enumerate(intrinsictypes)
    @mfalse centers[:, i] = jaccardcenter(collect(Xi[:, ind2type .== celltype]))[1]
end

# %%
centers = NamedArray(centers, names = (featurenames, intrinsictypes))

# %% [markdown]
# ### type radius
# average distance from cells in a cluster to cluster centers

# %%
@time celltocenter_dist = pairwise(jaccard, float(collect(Xi)), centers.array)

# %%
celltocenter_dist_avg = Ai'*celltocenter_dist./sum(Ai, dims=1)'

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries,
    )

# %%
radii = diag(celltocenter_dist_avg.array)

# r = [1:57, 58:114, 115:171, 172:228]
r = [1:46, 47:92, 93:138, 139:184, 185:230]
ncol = length(r)

plt = [plot(
        radii[r[i]], 
        permute = (:y, :x), 
        xticks = strings2ticks(intrinsictypes[r[i]]),
        yticks = [0, 0.2, 0.4, 0.6],
        xlim = (0.5, length(r[i]) + 0.5),
        ylim = (0, 0.7), 
        xflip = true, 
        ylabel = "radius", 
        yrotation = 45,
        seriestype = :stepmid,
        ytickfontsize = 7,
        xtickfontsize = 7,
        guidefontsize = 10,
        ) for i = 1:ncol]

h = histogram(radii, 
    xlabel = "type radius", 
    ylabel = "# types",
    tickfontsize = 12,
    guidefontsize = 14,
#    bottom_margin = 10mm
)

plot(plt..., h,
    size = (150*ncol, 600*5/ncol + 200), 
#    layout = (1, ncol),
    layout = @layout([grid(1, ncol){0.8h}; [b _]]),
    left_margin = 5mm
)

# %%
# Now using shared config.jl

# %%
savefig(joinpath(TARGETDIR, "Fig S3ab TypeRadii.pdf"))

# %%
#savefig("/Users/sseung/sseung@princeton.edu/OpticLobeCellTypesPaper/panels/RadiiTypes.svg")

# %%
#savefig("/Users/sseung/sseung@princeton.edu/OpticLobeCellTypesPaper/panels/TypeRadiusHistogram.svg")
