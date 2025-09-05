# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:percent
#     text_representation:
#       extension: .jl
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.17.3
#   kernelspec:
#     display_name: Julia 1.11.6
#     language: julia
#     name: julia-1.11
# ---

# %% [markdown]
# # Form vision paper (nonspatial figures)
# For reproducing figures in Seung, [Predicting visual function by interpreting a neuronal wiring diagram](https://doi.org/10.1038/s41586-024-07953-5), Nature 634:113-123 (2024).
#
# This notebook is restricted to figures not involving space.

# %%
using OpticLobe

include("config.jl")

# %%
using NamedArrays, SparseArrays

# %%
using Plots, Measures

# %%
using MissingsAsFalse, Missings

# %%
using StatsBase

# %%
using StatsPlots

# %%
using Distances, Clustering

# %%
# because labels requires a row vector
Base.transpose(s::String) = s

# %% [markdown]
# ## canonical ordering of types

# %%
Dm3types = ["Dm3v", "Dm3p", "Dm3q"]
TmYtypes = ["TmY4", "TmY9q", "TmY9q⊥"]
formtypes = vcat(Dm3types, TmYtypes)

# %% [markdown]
# ## cell numbers
# - Dm3 is more numerous than TmY
# - vertical is least numerous for Dm3
# - horizontal is most numerous for TmY

# %%
[Dm3types length.(type2ids.(Dm3types))]

# %%
[TmYtypes length.(type2ids.(TmYtypes))]

# %%
length(type2ids("Tm1"))

# %%
length(type2ids("T2a"))

# %% [markdown]
# ## Figure 1b-d. Visualizations of all Dm3v cells
# These show all cells of each type, rather than the subsets in the paper.

# %%
for Dm3type in Dm3types
    ng_hyper(Dm3type) |> display
end

# %% [markdown]
# ## Figure 2a-c. Visualizations of all TmY4 and TmY9 cells
# These show all cells of each type, rather than the subsets in the paper.

# %%
for TmYtype in TmYtypes
    ng_hyper(TmYtype) |> display
end

# %% [markdown]
# ## population-to-population connectivity
# two populations:
# - Dm3 includes trio of Dm3 types
# - TmY includes trio of TmY types

# %%
Wpp = [sum(Wtt[a, b]) for a in [Dm3types, TmYtypes], b in [Dm3types, TmYtypes]]  # 2x2 matrix of synapse numbers
Wpout = [sum(Wtc[a, :]) for a in [Dm3types, TmYtypes]]  # outgoing synapses from two populations
Wpin = [sum(Wct[:, a]) for a in [Dm3types, TmYtypes]]   # incoming synapses to two populations

# %%
plot(heatmap(Wpp, title = "number"), 
    heatmap(Wpp./Wpin', title = "in fraction"), 
    heatmap(Wpp./Wpout, title = "out fraction"),
    size = (1200, 250),
    layout = (1, 3),
    ticks = strings2ticks(["Dm3", "TmY"]),
    yflip = true
)

# %% [markdown]
# ## Figure 3a, type-to-type connectivity

# %%
# include LC in postsynaptic types
formtypespost = vcat(Dm3types, TmYtypes, ["LC15", "LC10e"])

# %%
NamedArray(collect(Wtt[formtypes, formtypespost]), names = (formtypes, formtypespost))

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 12,
    linewidth = 2,
    left_margin = 10mm,
    bottom_margin = 0mm,
    c = :hot, 
    xticks = strings2ticks(formtypes), yticks = strings2ticks(formtypes), 
    yflip = true
#    xrotation = 60,
    )

# %%
# version with cell types as tick labels
heatmap(infraction[formtypes, formtypespost], 
    yticks = strings2ticks(formtypes),
    xticks = strings2ticks(formtypespost), 
    xrot = 30, legend = :none
)

# %%
# version without tick labels for inclusion in figure file
heatmap(infraction[formtypes, formtypespost], xticks = strings2ticks(formtypespost), axis = ([], false), legend = :none)

# %%
savefig(joinpath(TARGETDIR, "Dm3TmYPopulationInteractions.pdf"))

# %% [markdown]
# ## Figure 3b. identify connections included in cartoon wiring diagram

# %%
NamedArray(collect(infraction[formtypes, formtypespost]), names = (formtypes, formtypespost)) |> showall

# %%
# connections included in cartoon wiring diagram Figure 3b
collect(infraction[formtypes, formtypespost] .> 0.03)

# %%
collect(infraction[formtypes, formtypespost] .> 0.03 .&& infraction[formtypes, formtypespost] .< 0.04)

# %% [markdown]
# ## functions for tracing pathways backwards and forwards

# %%
"""
find cells of type `celltype` presynaptic to cell `id`.
connections containing more than `thres` synapses
"""
function backwardtrace(id, celltype, thres = 4)
    strongpre = findall(W[:, id2ind(id)] .> thres)
    return @mfalse ind2id[strongpre[ind2type[strongpre] .== celltype]]
end

"""
find cells of type `celltype` postsynaptic to cell `id`.
connections containing more than `thres` synapses
"""
function forwardtrace(id, celltype, thres = 4)
    strongpost = findall(W[id2ind(id), :] .> thres)
    return @mfalse ind2id[strongpost[ind2type[strongpost] .== celltype]]
end

# %% [markdown]
# ## Figure 3c. Dm3q-Dm3p-TmY9q visualization

# %%
first(type2ids("Dm3p"), 10)

# %%
id = 720575940605138860 # Dm3p cell

# %%
cellids = vcat(backwardtrace(id, "Dm3q", 5), id, forwardtrace(id, "TmY9q", 5))

# %%
ind2type[id2ind.(cellids)]

# %%
ng_hyper(cellids)

# %%
W[Name.(cellids), Name.(cellids)] |> collect

# %%
[backwardtrace(id, "Tm1") for id in cellids]  # presynaptic Tm1 cells

# %% [markdown]
# ## Figure 3d. TmY4 visualization
# visualize TmY4-TmY4 pair reciprocally connected by more than thres synapses

# %%
ids = type2ids("TmY4")
WintraTmY4 = W[Name.(ids), Name.(ids)]
thres = 4
strongreciprocalpairs = findall(WintraTmY4 .> thres .&& WintraTmY4' .> thres)
println(length(strongreciprocalpairs), " reciprocal pairs connected by >", thres, " synapses")

# %%
# first pair in list
i, j = collect(Tuple(strongreciprocalpairs[1]))
println("$(ids[i])→$(ids[j]): ", WintraTmY4[i, j], " synapses")
println("$(ids[j])→$(ids[i]): ", WintraTmY4[j, i], " synapses")

# %%
ng_hyper(ids[[i, j]])

# %% [markdown]
# ## Figure 6a. LC15 visualization

# %%
id = 720575940604229152   # LC15 cell

# %%
ng_hyper(vcat([id], [backwardtrace(id, TmYtype) for TmYtype in TmYtypes]...))

# %% [markdown]
# ## Figure 6b. LC10e visualization

# %%
id = 720575940606274592   # LC10e cell

# %%
ng_hyper(vcat([id], [backwardtrace(id, TmYtype, 3) for TmYtype in TmYtypes]...))

# %% [markdown]
# ## Dm3p-Dm3q-TmY9q⊥ visualization

# %%
@mfalse ids = type2ids("Dm3q")

# %%
id = 720575940600623020  # Dm3q cell

# %%
backwardtrace(id, "Dm3p", 3)

# %%
forwardtrace(id, "TmY9q⊥", 5)

# %%
ng_hyper(vcat([id], backwardtrace(id, "Dm3p", 3), forwardtrace(id, "TmY9q⊥", 5)))

# %%
cellids = [720575940605123701, 720575940600623020, 720575940617250233]

# %%
ind2type[id2ind.(cellids)]

# %%
W[Name.(cellids), Name.(cellids)] |> collect

# %% [markdown]
# ## Figure S7. disynaptic pathways

# %%
# defined as numerous and receive L input + R7-8
hexel = ["L1", "L2", "L3", "L4", "L5", "Tm1", "Tm2", "Tm9", "Mi1", "Mi4", "Mi9"]
length.(type2ids.(hexel))

# %%
nothexel = setdiff(intrinsictypes, hexel)

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 10,
    legendfontsize = 10,
    linewidth = 2,
    left_margin = 0mm,
    bottom_margin = 0mm,
#    xrotation = 60,
    )

# %%
h = []
n = 12

for target in vcat(Dm3types, TmYtypes, "LC15", "LC10e")
    scores = [scorepath([c1, c2, target]) for c1 in hexel, c2 in nothexel]
    scores = NamedArray(scores, names = (hexel, nothexel))
#    println(target)
#    showall(sum(scores, dims = 2))
#    println("\n")
    scoresums = NamedArray(sum(scores.array, dims = 1)[:], nothexel)
    scoremaxs = NamedArray(maximum(scores.array, dims = 1)[:], nothexel)
    p = sortperm(scoresums, rev = true)

    push!(h, plot(100*[scoresums[p[1:n]], 
        scoremaxs[p[1:n]], scores["Tm1", p[1:n]]], 
        xticks = strings2ticks(nothexel[p[1:n]]), 
        label = ["sum" "max" "Tm1"],
        annotations = (:E, Plots.text(target, "Helvetica", :right)),
        permute = (:y, :x), 
        xlabel = startswith(target, "Dm3") ? "intermediary type" : "",
        leftmargin = startswith(target, "Dm3") ? 5mm : 0mm,
        xflip = true)
    )
end

# %%
plot(h[[1 4 7 2 5 8 3 6]]...,
    layout = grid(3, 3),
    size = (800, 1000),
    legend = :bottomright,
)

# %%
savefig(joinpath(TARGETDIR, "DisynapticInFractions.pdf"))

# %%
type2nt[["TmY10", "TmY11", "Tm20", "Tm25", "Tm27", "Tm8a", "Tm16", "Tm7"]]

# %% [markdown]
# ## list top partners (collapse)
# These are modified versions of `toppre` and `toppost` utility functions defined in `OpticLobe`.
#
# number of synapses per cell

# %%
Av = A[:, visualtypes]

# %%
# function toppost2(ids::Vector{Int64}, nresults = 15, reduce=false)
#     if reduce
#         pre = sort(NamedArray(mean(W[id2ind.(ids), :]*Av*collapse, dims=1)[:], collapsetypes), rev=true)
#     else    
#         pre = sort(NamedArray(mean(W[id2ind.(ids), :]*Av, dims=1)[:], visualtypes), rev=true)
#     end
#     (names(pre)[1][1:nresults], round.(Int64, pre.array[1:nresults]))
# end

# %%
function toppost2(ids::Vector{Int64}, nresults = 15, reduce=false)
    if reduce
        pre = sort(NamedArray(mean(W[id2ind.(ids), :]*Av*collapse, dims=1)[:], collapsetypes), rev=true)
    else    
        pre = sort(NamedArray(mean(W[id2ind.(ids), :]*Av, dims=1)[:], visualtypes), rev=true)
    end
    (disallowmissing(names(pre, 1))[1:nresults], round.(Int64, pre.array[1:nresults]))
end

# %%
function toppre2(ids::Vector{Int64}, nresults = 15, reduce=false)
    if reduce
        post = sort(NamedArray(mean(collapse'*Av'*W[:, id2ind.(ids)], dims=2)[:], collapsetypes), rev=true)
    else
        post = sort(NamedArray(mean(Av'*W[:, id2ind.(ids)], dims=2)[:], visualtypes), rev=true)
    end
    (disallowmissing(names(post, 1))[1:nresults], round.(Int64, post.array[1:nresults]))
end

# %% [markdown]
# ## create matrix that collapses all Dm3 types to "Dm3", and all TmY4/TmY9 types to "TmY"

# %%
using LinearAlgebra

# %%
ntype = length(visualtypes)

# %%
mask = NamedArray(trues(ntype), visualtypes)

mask["TmY9q"] = false
mask["TmY9q⊥"] = false
mask["Dm3p"] = false
mask["Dm3q"] = false

# %%
collapse = NamedArray(Matrix(I, ntype, ntype), names = (visualtypes, visualtypes))

collapse = collapse[:, mask]

# %%
collapse["Dm3p", "Dm3v"] = 1 
collapse["Dm3q", "Dm3v"] = 1 
collapse["TmY9q", "TmY4"] = 1 
collapse["TmY9q⊥", "TmY4"] = 1 

# %%
collapsetypes = visualtypes[mask]

# %%
collapsetypes[collapsetypes .== "TmY4"] .= "TmY4/9"

# %%
collapsetypes[collapsetypes .== "Dm3v"] .= "Dm3"

# %%
collapse = NamedArray(collapse.array, names = (visualtypes, collapsetypes))

# %%
visualtypes[collapse[:, "TmY4/9"] .> 0]

# %% [markdown]
# ## Figure S1. Dm3 inputs and outputs

# %%
@mfalse Dm3ids = vcat([ind2id[ind2type .== celltype] for celltype in ["Dm3v", "Dm3p", "Dm3q"]]...)

# top presynaptic and postsynaptic partners of Dm3
npartners = 20
@mfalse Dm3prestrings = toppre2(Dm3ids, npartners)[1]
@mfalse Dm3poststrings = toppost2(Dm3ids, npartners)[1]

Dm3strings = ["Dm3v", "Dm3p", "Dm3q"]
Dm3inputs = collect(A[:, Dm3prestrings]'*W*A[:, Dm3strings])
Dm3inputs = Dm3inputs./collect(sum(A[:, Dm3strings], dims=1))

Dm3outputs = collect(A[:, Dm3strings]'*W*A[:, Dm3poststrings])'
Dm3outputs = Dm3outputs./collect(sum(A[:, Dm3strings], dims=1))

# %%
sum(Dm3outputs, dims = 2)

# %%
@mfalse Dm3prestringsreduced = toppre2(Dm3ids, npartners, true)[1]
@mfalse Dm3poststringsreduced = toppost2(Dm3ids, npartners, true)[1]

Dm3inputsreduced = collect((Av*collapse)[:, Dm3prestringsreduced]'*W*Av[:, Dm3strings])
Dm3inputsreduced = Dm3inputsreduced./collect(sum(Av[:, Dm3strings], dims=1))

Dm3outputsreduced = collect(Av[:, Dm3strings]'*W*(Av*collapse)[:, Dm3poststringsreduced])'
Dm3outputsreduced = Dm3outputsreduced./collect(sum(Av[:, Dm3strings], dims=1))

# %% [markdown]
# ### Plots.jl defaults

# %%
Plots.reset_defaults()

# %%
default(; # Plots defaults
    fontfamily="Helvetica",
    label="", # only explicit legend entries
    guidefontsize = 12,
    tickfontsize = 12,
    legendfontsize = 12,
    linewidth = 2,
    left_margin = 10mm,
    bottom_margin = 0mm,
#    xrotation = 60,
    )

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(Dm3strings), 
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true 
    )
for (data, partnertypes, direction) in zip(
        [Dm3inputs, Dm3outputs, Dm3inputsreduced, Dm3outputsreduced], 
        [Dm3prestrings, Dm3poststrings, Dm3prestringsreduced, Dm3poststringsreduced], 
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "Dm3InOut.pdf"))

# %% [markdown]
# ## dendrogram for Dm3 and TmY types

# %%
embeddings = vcat(infraction[visualtypes, formtypes], outfraction[formtypes, visualtypes]')

# %%
dist = pairwise(jaccard, embeddings)
c = hclust(dist, linkage = :average)
plot(c, xticks = strings2ticks(formtypes[c.order]), size = (1200, 300), xrot = 30, bottom_margin = 7mm)

# %%
clist = []
for embeddings in [vcat(infraction[celltypes, formtypes], outfraction[formtypes, celltypes]') for celltypes in [intrinsictypes, visualtypes, alltypes]]
    dist = pairwise(jaccard, embeddings)
    push!(clist, hclust(dist, linkage = :average))
end

# %%
h = [plot(c, xticks = strings2ticks(formtypes[c.order]), size = (1200, 300), xrot = 30, bottom_margin = 7mm) for c in clist]

# %%
plot(h..., layout = (1, 3))

# %% [markdown]
# ## dendrogram for Dm3, TmY, and other types connected with them

# %%
involved = vcat(formtypes,[ "L3", "L5", "Tm1", "Tm2", "Tm9", "Mi2", "Mi4", "Mi9", "T2a", "T3", "Tm20", "Tm5f", "Dm15", "Dm18", "Tm21", "Tm25", "Tm7", "Tm8a", "Tm16", "TmY10", "TmY11", "Y3", "T5a", "T4a", "T5c", "T5d", "T4c", "T4d"])

# %%
embeddings = vcat(infraction[visualtypes, involved], outfraction[involved, visualtypes]')

# %%
dist = pairwise(jaccard, embeddings)
c = hclust(dist, linkage = :average)
plot(c, xticks = strings2ticks(involved[c.order]), size = (1200, 300), xrot = 30, bottom_margin = 7mm)

# %% [markdown]
# ## Figure S2. Tm1 inputs and outputs

# %%
# top presynaptic and postsynaptic partners of Tm1
#npartners = 25
npartners = 26
@mfalse Tm1prestrings = toppre2(ind2id[ind2type .== "Tm1"], npartners)[1]
@mfalse Tm1poststrings = toppost2(ind2id[ind2type .== "Tm1"], npartners)[1]

Tm1inputs = collect(A[:, Tm1prestrings]'*W*A[:, ["Tm1"]])
Tm1inputs = Tm1inputs./collect(sum(A[:, ["Tm1"]], dims=1))
Tm1outputs = collect(A[:, ["Tm1"]]'*W*A[:, Tm1poststrings])'
Tm1outputs = Tm1outputs./collect(sum(A[:, ["Tm1"]], dims=1))

# %%
# top presynaptic and postsynaptic partners of Tm1
@mfalse Tm1prestringsreduced = toppre2(ind2id[ind2type .== "Tm1"], npartners, true)[1]
@mfalse Tm1poststringsreduced = toppost2(ind2id[ind2type .== "Tm1"], npartners, true)[1]

Tm1inputsreduced = collect((Av*collapse)[:, Tm1prestringsreduced]'*W*Av[:, ["Tm1"]])
Tm1inputsreduced = Tm1inputsreduced./collect(sum(Av[:, ["Tm1"]], dims=1))
Tm1outputsreduced = collect(Av[:, ["Tm1"]]'*W*Av*collapse[:, Tm1poststringsreduced])[:]
Tm1outputsreduced = Tm1outputsreduced./collect(sum(Av[:, ["Tm1"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = "Tm1", 
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, direction) in zip(
        [Tm1inputs, Tm1outputs, Tm1inputsreduced, Tm1outputsreduced], 
        [Tm1prestrings, Tm1poststrings, Tm1prestringsreduced, Tm1poststringsreduced], 
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "Tm1InOut.pdf"))

# %% [markdown]
# ## Figure S4. TmY4 and TmY9 four panels

# %% [markdown]
# ### TmY9

# %%
TmY9strings = ["TmY9q", "TmY9q⊥"]
@mfalse TmY9ids = vcat([ind2id[ind2type .== celltype] for celltype in TmY9strings]...)

# top presynaptic and postsynaptic partners of TmY
npartners = 30
@mfalse TmY9prestrings = toppre2(TmY9ids, npartners)[1]

@mfalse TmY9poststrings = toppost2(TmY9ids, npartners)[1]

TmY9inputs = collect(A[:, TmY9prestrings]'*W*A[:, TmY9strings])
TmY9inputs = TmY9inputs./collect(sum(A[:, TmY9strings], dims=1))
TmY9outputs = collect(A[:, TmY9strings]'*W*A[:, TmY9poststrings])'
TmY9outputs = TmY9outputs./collect(sum(A[:, TmY9strings], dims=1))

# %% [markdown]
# ### TmY4

# %%
# top presynaptic and postsynaptic partners of TmY4
npartners = 20
@mfalse TmY4prestrings = toppre2(ind2id[ind2type .== "TmY4"], npartners)[1]
@mfalse TmY4poststrings = toppost2(ind2id[ind2type .== "TmY4"], npartners)[1]

TmY4inputs = collect(Av[:, TmY4prestrings]'*W*Av[:, ["TmY4"]])
TmY4inputs = TmY4inputs./collect(sum(Av[:, ["TmY4"]], dims=1))
TmY4outputs = collect(Av[:, ["TmY4"]]'*W*Av[:, TmY4poststrings])'
TmY4outputs = TmY4outputs./collect(sum(Av[:, ["TmY4"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [TmY4inputs, TmY4outputs, TmY9inputs, TmY9outputs], 
        [TmY4prestrings, TmY4poststrings, TmY9prestrings, TmY9poststrings], 
        [ "TmY4", "TmY4", TmY9strings, TmY9strings],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.4, 0.6]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "TmY4TmY9InOut.pdf"))

# %% [markdown]
# ## Figure S5. TmY4 and TmY9 four panels, reduced

# %% [markdown]
# ### TmY9

# %%
TmY9strings = ["TmY9q", "TmY9q⊥"]
@mfalse TmY9ids = vcat([ind2id[ind2type .== celltype] for celltype in TmY9strings]...)

# top presynaptic and postsynaptic partners of TmY
npartners = 30
@mfalse TmY9prestringsreduced = toppre2(TmY9ids, npartners, true)[1]

@mfalse TmY9poststringsreduced = toppost2(TmY9ids, npartners, true)[1]

TmY9inputsreduced = collect((Av*collapse)[:, TmY9prestringsreduced]'*W*Av[:, TmY9strings])
TmY9inputsreduced = TmY9inputsreduced./collect(sum(Av[:, TmY9strings], dims=1))
TmY9outputsreduced = collect(Av[:, TmY9strings]'*W*(Av*collapse)[:, TmY9poststringsreduced])'
TmY9outputsreduced = TmY9outputsreduced./collect(sum(Av[:, TmY9strings], dims=1))

# %% [markdown]
# ### TmY4

# %%
# top presynaptic and postsynaptic partners of TmY4
npartners = 20
@mfalse TmY4prestringsreduced = toppre2(ind2id[ind2type .== "TmY4"], npartners, true)[1]
@mfalse TmY4poststringsreduced  = toppost2(ind2id[ind2type .== "TmY4"], npartners, true)[1]

TmY4inputsreduced = collect((Av*collapse)[:, TmY4prestringsreduced]'*W*Av[:, ["TmY4"]])
TmY4inputsreduced = TmY4inputsreduced./collect(sum(Av[:, ["TmY4"]], dims=1))
TmY4outputsreduced = collect(Av[:, ["TmY4"]]'*W*(Av*collapse)[:, TmY4poststringsreduced])'
TmY4outputsreduced = TmY4outputsreduced./collect(sum(Av[:, ["TmY4"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [TmY4inputsreduced, TmY4outputsreduced, TmY9inputsreduced, TmY9outputsreduced], 
        [TmY4prestringsreduced, TmY4poststringsreduced, TmY9prestringsreduced, TmY9poststringsreduced], 
        [ "TmY4", "TmY4", TmY9strings, TmY9strings],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.4, 0.6]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "TmY4TmY9InOutReduced.pdf"))

# %% [markdown]
# ## Y3 and TmY5a

# %%
# top presynaptic and postsynaptic partners of Y3
npartners = 25
@mfalse Y3prestrings = toppre2(ind2id[ind2type .== "Y3"], npartners)[1]
@mfalse Y3poststrings = toppost2(ind2id[ind2type .== "Y3"], npartners)[1]

Y3inputs = collect(Av[:, Y3prestrings]'*W*Av[:, ["Y3"]])
Y3inputs = Y3inputs./collect(sum(Av[:, ["Y3"]], dims=1))
Y3outputs = collect(Av[:, ["Y3"]]'*W*Av[:, Y3poststrings])'
Y3outputs = Y3outputs./collect(sum(Av[:, ["Y3"]], dims=1))

# %%
# top presynaptic and postsynaptic partners of TmY5a
npartners = 25
@mfalse TmY5aprestrings = toppre2(ind2id[ind2type .== "TmY5a"], npartners)[1]
@mfalse TmY5apoststrings = toppost2(ind2id[ind2type .== "TmY5a"], npartners)[1]

TmY5ainputs = collect(Av[:, TmY5aprestrings]'*W*Av[:, ["TmY5a"]])
TmY5ainputs = TmY5ainputs./collect(sum(Av[:, ["TmY5a"]], dims=1))
TmY5aoutputs = collect(Av[:, ["TmY5a"]]'*W*Av[:, TmY5apoststrings])'
TmY5aoutputs = TmY5aoutputs./collect(sum(Av[:, ["TmY5a"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [Y3inputs, Y3outputs, TmY5ainputs, TmY5aoutputs], 
        [Y3prestrings, Y3poststrings, TmY5aprestrings, TmY5apoststrings], 
        [ "Y3", "Y3", "TmY5a", "TmY5a"],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "Y3TmY5aInOut.pdf"))

# %% [markdown]
# ## T2a

# %%
T2astrings = ["T2a"]
@mfalse T2aids = vcat([ind2id[ind2type .== celltype] for celltype in T2astrings]...)

# top presynaptic and postsynaptic partners of TmY
npartners = 20
@mfalse T2aprestrings = toppre2(T2aids, npartners)[1]

@mfalse T2apoststrings = toppost2(T2aids, npartners)[1]

T2ainputs = collect(Av[:, T2aprestrings]'*W*Av[:, T2astrings])
T2ainputs = T2ainputs./collect(sum(Av[:, T2astrings], dims=1))
T2aoutputs = collect(Av[:, T2astrings]'*W*Av[:, T2apoststrings])'
T2aoutputs = T2aoutputs./collect(sum(Av[:, T2astrings], dims=1))

# %%
plot(
    plot(T2ainputs, 
    xticks = (1:length(T2aprestrings), T2aprestrings), 
    labels = "T2a",
    xlabel = "input cell type",
    ylabel = "# of input synapses",
    permute = (:x, :y), xflip = true,
#    left_margin = 10mm,
        ),
plot(T2aoutputs, 
    xticks = (1:length(T2apoststrings), T2apoststrings), 
    labels = "T2a",
    xlabel = "output cell type",
    ylabel = "# of output synapses",
    permute = (:x, :y), xflip = true,
#    left_margin = 5mm,
), 
    layout = grid(1, 2),
    bottom_margin = 5mm,
    size = (800, 500),
    legend = :bottomright
    )

# %%
savefig(joinpath(TARGETDIR, "T2aInOut.pdf"))

# %% [markdown]
# ## Figure S8. Y3 and T2a

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [Y3inputs, Y3outputs, T2ainputs, T2aoutputs], 
        [Y3prestrings, Y3poststrings, T2aprestrings, T2apoststrings], 
        [ "Y3", "Y3", "T2a", "T2a"],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "Y3T2aInOut.pdf"))

# %% [markdown]
# ## Figure S9. TmY10 TmY11

# %%
# top presynaptic and postsynaptic partners of TmY10
npartners = 20
@mfalse TmY10prestrings = toppre2(ind2id[ind2type .== "TmY10"], npartners)[1]
@mfalse TmY10poststrings = toppost2(ind2id[ind2type .== "TmY10"], npartners)[1]

TmY10inputs = collect(Av[:, TmY10prestrings]'*W*Av[:, ["TmY10"]])
TmY10inputs = TmY10inputs./collect(sum(Av[:, ["TmY10"]], dims=1))
TmY10outputs = collect(Av[:, ["TmY10"]]'*W*Av[:, TmY10poststrings])'
TmY10outputs = TmY10outputs./collect(sum(Av[:, ["TmY10"]], dims=1))

# %%
# top presynaptic and postsynaptic partners of TmY10
npartners = 20
@mfalse TmY11prestrings = toppre2(ind2id[ind2type .== "TmY11"], npartners)[1]
@mfalse TmY11poststrings = toppost2(ind2id[ind2type .== "TmY11"], npartners)[1]

TmY11inputs = collect(Av[:, TmY11prestrings]'*W*Av[:, ["TmY11"]])
TmY11inputs = TmY11inputs./collect(sum(Av[:, ["TmY11"]], dims=1))
TmY11outputs = collect(Av[:, ["TmY11"]]'*W*Av[:, TmY11poststrings])'
TmY11outputs = TmY11outputs./collect(sum(Av[:, ["TmY11"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [TmY10inputs, TmY10outputs, TmY11inputs, TmY11outputs], 
        [TmY10prestrings, TmY10poststrings, TmY11prestrings, TmY11poststrings], 
        [ "TmY10", "TmY10", "TmY11", "TmY11"],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (750, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "TmY10TmY11InOut.pdf"))

# %% [markdown]
# ## Figure S10. LC10e and LC15

# %%
### new clustering
LC10ed = [720575940603823136, 720575940604935089, 720575940610075029, 720575940614496735, 720575940614746178, 720575940616179050, 720575940618360741, 720575940620416176, 720575940621175979, 720575940622744414, 720575940622805854, 720575940623979815, 720575940624787504, 720575940625143396, 720575940625981080, 720575940628023296, 720575940628359783, 720575940628615752, 720575940629385943, 720575940629398492, 720575940629958764, 720575940631576546, 720575940633027615, 720575940637401321, 720575940637706574, 720575940637731955, 720575940638701375, 720575940644393070]
LC10ev = [720575940604125920, 720575940606274592, 720575940608977988, 720575940609733387, 720575940609890402, 720575940612305250, 720575940612332401, 720575940614897662, 720575940621233157, 720575940621602047, 720575940622486057, 720575940622602992, 720575940623601193, 720575940624050407, 720575940624595893, 720575940628632064, 720575940629076042, 720575940629756556, 720575940630438651, 720575940631305467, 720575940632034385, 720575940634088602, 720575940640974555]

# %%
# top presynaptic and postsynaptic partners of LC10e
npartners = 25
@mfalse LC10eprestrings = toppre2(vcat(LC10ed, LC10ev), npartners)[1]
@mfalse LC10epoststrings = toppost2(vcat(LC10ed, LC10ev), npartners)[1]

LC10evinputs = mean(Av[:, LC10eprestrings]'*W[:, id2ind.(LC10ev)], dims=2)[:]
LC10evoutputs = mean(W[id2ind.(LC10ev), :]*Av[:, LC10epoststrings], dims=1)[:]

LC10edinputs = mean(Av[:, LC10eprestrings]'*W[:, id2ind.(LC10ed)], dims=2)[:]
LC10edoutputs = mean(W[id2ind.(LC10ed), :]*Av[:, LC10epoststrings], dims=1)[:]

# %%
# top presynaptic and postsynaptic partners of LC15
npartners = 25
@mfalse LC15prestrings = toppre2(ind2id[ind2type .== "LC15"], npartners)[1]
@mfalse LC15poststrings = toppost2(ind2id[ind2type .== "LC15"], npartners)[1]

LC15inputs = collect(Av[:, LC15prestrings]'*W*Av[:, ["LC15"]])
LC15inputs = LC15inputs./collect(sum(Av[:, ["LC15"]], dims=1))
LC15outputs = collect(Av[:, ["LC15"]]'*W*Av[:, LC15poststrings])'
LC15outputs = LC15outputs./collect(sum(Av[:, ["LC15"]], dims=1))

# %%
h = [plot(data, xticks = strings2ticks(partnertypes), 
        labels = transpose(visualtypes),
        xlabel = "$direction cell type", 
        ylabel = "# $direction synapses", 
        permute = (:x, :y), xflip = true,
        ylim = (0, Inf), widen = true) 
for (data, partnertypes, visualtypes, direction) in zip(
        [[LC10edinputs, LC10evinputs], [LC10edoutputs, LC10evoutputs], LC15inputs, LC15outputs], 
        [LC10eprestrings, LC10epoststrings, LC15prestrings, LC15poststrings], 
        [ ["LC10e dorsal", "LC10e ventral"], ["LC10e dorsal", "LC10e ventral"], "LC15", "LC15"],
        ["input", "output", "input", "output"]
        )
    ]
plot(h...,
    layout = grid(2, 2, heights=[0.5, 0.5]),
    size = (800, 1000),
    legend = :bottomright
)

# %%
savefig(joinpath(TARGETDIR, "LC10eLC15InOut.svg"))

# %% [markdown]
# ## search for VPNs with high infraction from TmY

# %%
thres = 0.045

# %%
abovethresholdtypes = unique(visualtypes[sum(infraction[TmYtypes, visualtypes], dims = 1)[:] .> 0.08])

# %%
upstreamTmY = NamedArray(sum(infraction[TmYtypes, visualtypes], dims = 1)[:].array, visualtypes)

# %%
vpntypes = vcat(filter(startswith("LC"), visualtypes), filter(startswith("LT"), visualtypes))

# %%
rankedlist = sort(upstreamTmY[vpntypes], rev=true)

# %%
for celltype in names(rankedlist)[1][1:20]
    @mfalse println(sum(ind2type .== celltype))
    toppre(celltype)
end

# %%
toppre("LTe57")

# %%
toppre("LTe37")

# %%
toppre("LC45")

# %%
toppre("LT84")

# %%
toppre("LT86")

# %%
toppre("LC46")

# %%
toppre("LTe35")

# %%
first(sort(infraction["TmY9q", :], rev = true), 20)

# %%
toppre("LTe63")

# %%
toppre("TmY9q", nresults = 30)

# %%
toppre("LTe08")

# %%
# toppre("CB0472") # doesn't exist in latest annotations

# %%
toppre("LTe28")

# %%
toppre("LTe19")

# %%
toppre("LTe45")

# %%
toppre("LC41")

# %%
toppre("LC10e")

# %%
toppre("LT80")

# %%
toppre("LT72")

# %%
first(sort(infraction["TmY9q⊥", :], rev = true), 20)

# %%
toppre("LTe08")

# %%
for id in type2ids("LCe08")
    toppre([id])
end
