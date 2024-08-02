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
# # Figure 2a - feature vectors for cells

# %% jupyter={"outputs_hidden": false}
using OpticLobe
using NamedArrays
using Distances
using Plots, Measures, Printf

# %% [markdown]
# ## projection vectors (intrinsic only)

# %% jupyter={"outputs_hidden": false}
Xi = vcat(Ai'*W, (W*Ai)')
featurenames = vcat("in-" .* intrinsictypes, "out-" .* intrinsictypes)
Xi = NamedArray(Xi, names = (featurenames, ind2id))

# %% [markdown]
# ## visualize projection vectors for three example cells
# Note that the ylabels are now "in-Am1" instead of "in-C2" as the intrinsic types are now sorted alphabetically

# %%
ids = [720575940606527618, 720575940623179177, 720575940607643058]

# %% jupyter={"outputs_hidden": false}
titles = [@sprintf("Cell %d\n%s", i, "($celltype)") for (i, celltype) in enumerate(ind2type[id2ind.(ids)])]

# %% jupyter={"outputs_hidden": false}
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries
    xtickfontsize = 9, 
    ytickfontsize = 9,
    guidefontsize = 10,
    titlefontsize = 9
    )

# %% jupyter={"outputs_hidden": false}
skipping = [1, length(intrinsictypes)+1]
plt = [
    plot(
        Xi[:, Name.(ids[i])]/100,
        permute = (:y, :x),
        seriestype = :stepmid,
        xflip = true,
        ylabel = i == 2 ? "# synapses/100" : "",
        xticks = i == 1 ? (skipping, featurenames[skipping]) : false,
#    yformatter=Returns(""),
#    left_margin = 10mm,
        xlim = (0.5, 0.5 + length(featurenames)),
        yticks = [[0, 5], [0, 10], [0, 3]][i],
        title = titles[i],
        margin = 0mm,
        xrotation = 45,
    ) for i in 1:3]
for i = 1:3
    annotate!(plt[i], :bottomright, Plots.text(ids[i], "Helvetica", 8, rotation=270))
end

# %%
p = plot(
    plt...,
    size = (250, 800),
    layout = (1, 3),
    left_margin = 7mm
    )
for i=1:3
    hspan!(p[i], skipping .+ 0.5; alpha = 0.1)
    hspan!(p[i], skipping .+ 0.5 .+ length(intrinsictypes); alpha = 0.1)
end
display(p)

# %% jupyter={"outputs_hidden": false}
# savefig("/Users/sseung/sseung@princeton.edu/OpticLobeCellTypesPaper/panels/ThreeFeatureVectors.svg")
savefig("Fig 2a.svg")

# %% [markdown]
# ## Jaccard distances between three example cells

# %% jupyter={"outputs_hidden": false}
jaccard(Xi[:, Name.(ids[1])], Xi[:, Name.(ids[2])])

# %% jupyter={"outputs_hidden": false}
jaccard(Xi[:, Name.(ids[1])], Xi[:, Name.(ids[3])])

# %% jupyter={"outputs_hidden": false}
jaccard(Xi[:, Name.(ids[2])], Xi[:, Name.(ids[3])])
